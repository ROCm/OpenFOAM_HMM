/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2018 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sampledIsoSurfaceTopo.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "isoSurfaceTopo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledIsoSurfaceTopo, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledIsoSurfaceTopo,
        word,
        isoSurfaceTopo
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledIsoSurfaceTopo::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // No update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    prevTimeIndex_ = fvm.time().timeIndex();

    // Clear derived data
    sampledSurface::clearGeom();

    // Use field from database, or try to read it in

    const auto* cellFldPtr = fvm.findObject<volScalarField>(isoField_);

    if (debug)
    {
        if (cellFldPtr)
        {
            InfoInFunction << "Lookup " << isoField_ << endl;
        }
        else
        {
            InfoInFunction
                << "Reading " << isoField_
                << " from time " << fvm.time().timeName()
                << endl;
        }
    }

    // For holding the volScalarField read in.
    autoPtr<volScalarField> fieldReadPtr;

    if (!cellFldPtr)
    {
        // Bit of a hack. Read field and store.

        fieldReadPtr = autoPtr<volScalarField>::New
        (
            IOobject
            (
                isoField_,
                fvm.time().timeName(),
                fvm,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvm
        );
    }

    const volScalarField& cellFld =
        (fieldReadPtr.valid() ? *fieldReadPtr : *cellFldPtr);

    auto tpointFld = volPointInterpolation::New(fvm).interpolate(cellFld);

    //- Direct from cell field and point field. Gives bad continuity.
    isoSurfaceTopo surf
    (
        fvm,
        cellFld.primitiveField(),
        tpointFld().primitiveField(),
        isoVal_,
        (regularise_ ? isoSurfaceTopo::DIAGCELL : isoSurfaceTopo::NONE)
    );

    MeshedSurface<face>& mySurface = const_cast<sampledIsoSurfaceTopo&>(*this);

    mySurface.transfer(static_cast<meshedSurface&>(surf));
    meshCells_ = std::move(surf.meshCells());

    // triangulate uses remapFaces()
    // - this is somewhat less efficient since it recopies the faces
    // that we just created, but we probably don't want to do this
    // too often anyhow.
    if (triangulate_)
    {
        labelList faceMap;
        mySurface.triangulate(faceMap);
        meshCells_ = UIndirectList<label>(meshCells_, faceMap)();
    }

    if (debug)
    {
        Pout<< "sampledIsoSurfaceTopo::updateGeometry() : constructed iso:"
            << nl
            << "    regularise     : " << regularise_ << nl
            << "    triangulate    : " << triangulate_ << nl
            << "    isoField       : " << isoField_ << nl
            << "    isoValue       : " << isoVal_ << nl
            << "    points         : " << points().size() << nl
            << "    faces          : " << MeshStorage::size() << nl
            << "    cut cells      : " << meshCells_.size() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledIsoSurfaceTopo::sampledIsoSurfaceTopo
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    MeshStorage(),
    isoField_(dict.get<word>("isoField")),
    isoVal_(dict.get<scalar>("isoValue")),
    regularise_(dict.lookupOrDefault("regularise", true)),
    triangulate_(dict.lookupOrDefault("triangulate", false)),
    prevTimeIndex_(-1),
    meshCells_()
{
    if (triangulate_ && !regularise_)
    {
        FatalIOErrorInFunction(dict) << "Cannot both use regularise"
            << " and triangulate" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledIsoSurfaceTopo::~sampledIsoSurfaceTopo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledIsoSurfaceTopo::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledIsoSurfaceTopo::expire()
{
    // Clear derived data
    sampledSurface::clearGeom();

    // Already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // Force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledIsoSurfaceTopo::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurfaceTopo::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurfaceTopo::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurfaceTopo::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurfaceTopo::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurfaceTopo::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledIsoSurfaceTopo::print(Ostream& os) const
{
    os  << "sampledIsoSurfaceTopo: " << name() << " :"
        << "  field:" << isoField_
        << "  value:" << isoVal_;
        //<< "  faces:" << faces().size()   // possibly no geom yet
        //<< "  points:" << points().size();
}


// ************************************************************************* //
