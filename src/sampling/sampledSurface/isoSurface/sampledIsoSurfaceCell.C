/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "sampledIsoSurfaceCell.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "isoSurfaceCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledIsoSurfaceCell, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledIsoSurfaceCell,
        word,
        isoSurfaceCell
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledIsoSurfaceCell::updateGeometry() const
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
        (fieldReadPtr ? *fieldReadPtr : *cellFldPtr);

    auto tpointFld = volPointInterpolation::New(fvm).interpolate(cellFld);

    // Non-averaged? Use reference
    tmp<scalarField> tcellValues(cellFld.primitiveField());

    if (average_)
    {
        // From point field and interpolated cell.
        tcellValues = tmp<scalarField>::New(fvm.nCells(), Zero);
        auto& cellAvg = tcellValues.ref();

        labelField nPointCells(fvm.nCells(), Zero);

        for (label pointi = 0; pointi < fvm.nPoints(); ++pointi)
        {
            const scalar& val = tpointFld().primitiveField()[pointi];
            const labelList& pCells = fvm.pointCells(pointi);

            for (const label celli : pCells)
            {
                cellAvg[celli] += val;
                ++nPointCells[celli];
            }
        }
        forAll(cellAvg, celli)
        {
            cellAvg[celli] /= nPointCells[celli];
        }
    }


    meshedSurface& mySurface = const_cast<sampledIsoSurfaceCell&>(*this);

    isoSurfaceCell surf
    (
        fvm,
        tcellValues(),
        tpointFld().primitiveField(),
        isoVal_,
        isoParams_
    );

    mySurface.transfer(static_cast<meshedSurface&>(surf));
    meshCells_.transfer(surf.meshCells());

    if (debug)
    {
        Pout<< "isoSurfaceCell::updateGeometry() : constructed iso:"
            << nl
            << "    isoField       : " << isoField_ << nl
            << "    isoValue       : " << isoVal_ << nl
            << "    average        : " << Switch(average_) << nl
            << "    filter         : "
            << Switch(bool(isoParams_.filter())) << nl
            << "    bounds         : " << isoParams_.getClipBounds() << nl
            << "    points         : " << points().size() << nl
            << "    faces          : " << Mesh::size() << nl
            << "    cut cells      : " << meshCells_.size() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledIsoSurfaceCell::sampledIsoSurfaceCell
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    Mesh(),
    isoField_(dict.get<word>("isoField")),
    isoVal_(dict.get<scalar>("isoValue")),
    isoParams_(dict),
    average_(dict.getOrDefault("average", true)),
    prevTimeIndex_(-1),
    meshCells_()
{
    isoParams_.algorithm(isoSurfaceParams::ALGO_CELL);  // Force
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledIsoSurfaceCell::~sampledIsoSurfaceCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledIsoSurfaceCell::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledIsoSurfaceCell::expire()
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


bool Foam::sampledIsoSurfaceCell::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurfaceCell::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurfaceCell::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurfaceCell::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurfaceCell::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurfaceCell::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledIsoSurfaceCell::print(Ostream& os) const
{
    os  << "isoSurfaceCell: " << name() << " :"
        << "  field:" << isoField_
        << "  value:" << isoVal_;
        //<< "  faces:" << faces().size()   // possibly no geom yet
        //<< "  points:" << points().size();
}


// ************************************************************************* //
