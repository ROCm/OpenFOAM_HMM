/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
#include "isoSurfaceTopo.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

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

    // Clear any previously stored topologies
    surface_.clear();
    meshCells_.clear();
    isoSurfacePtr_.reset(nullptr);

    // Clear derived data
    sampledSurface::clearGeom();


    // Handle cell zones as inverse (blocked) selection
    if (!ignoreCellsPtr_)
    {
        ignoreCellsPtr_.reset(new bitSet);

        if (-1 != mesh().cellZones().findIndex(zoneNames_))
        {
            bitSet select(mesh().cellZones().selection(zoneNames_));

            if (select.any() && !select.all())
            {
                // From selection to blocking
                select.flip();

                *ignoreCellsPtr_ = std::move(select);
            }
        }
    }


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

    // Field reference (assuming non-averaged)
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

    {
        isoSurfaceTopo surf
        (
            fvm,
            cellFld.primitiveField(),
            tpointFld().primitiveField(),
            isoVal_,
            isoParams_,
            *ignoreCellsPtr_
        );

        surface_.transfer(static_cast<meshedSurface&>(surf));
        meshCells_.transfer(surf.meshCells());
    }

    // if (subMeshPtr_ && meshCells_.size())
    // {
    //     // With the correct addressing into the full mesh
    //     meshCells_ =
    //         UIndirectList<label>(subMeshPtr_->cellMap(), meshCells_);
    // }


    // triangulate uses remapFaces()
    // - this is somewhat less efficient since it recopies the faces
    // that we just created, but we probably don't want to do this
    // too often anyhow.
    if (triangulate_ && surface_.size())
    {
        labelList faceMap;
        surface_.triangulate(faceMap);
        meshCells_ = UIndirectList<label>(meshCells_, faceMap)();
    }

    if (debug)
    {
        Pout<< "isoSurfaceTopo::updateGeometry() : constructed iso:" << nl
            << "    isoField       : " << isoField_ << nl
            << "    isoValue       : " << isoVal_ << nl
            << "    average        : " << Switch(average_) << nl
            << "    filter         : "
            << isoSurfaceParams::filterNames[isoParams_.filter()] << nl
            << "    triangulate    : " << Switch(triangulate_) << nl
            << "    bounds         : " << isoParams_.getClipBounds() << nl
            << "    points         : " << points().size() << nl
            << "    faces          : " << surface().size() << nl
            << "    cut cells      : " << meshCells().size() << endl;
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
    isoField_(dict.get<word>("isoField")),
    isoVal_(dict.get<scalar>("isoValue")),
    isoParams_(dict),
    average_(dict.getOrDefault("average", false)),
    triangulate_(dict.getOrDefault("triangulate", false)),
    zoneNames_(),
    prevTimeIndex_(-1),
    surface_(),
    meshCells_(),
    isoSurfacePtr_(nullptr)
{
    isoParams_.algorithm(isoSurfaceParams::ALGO_TOPO);  // Force

    if
    (
        triangulate_
     && (isoParams_.filter() == isoSurfaceParams::filterType::NONE)
    )
    {
        FatalIOErrorInFunction(dict)
            << "Cannot triangulate without a regularise filter" << nl
            << exit(FatalIOError);
    }

    if (!dict.readIfPresent("zones", zoneNames_) && dict.found("zone"))
    {
        zoneNames_.resize(1);
        dict.readEntry("zone", zoneNames_.first());
    }

    if (-1 != mesh.cellZones().findIndex(zoneNames_))
    {
        DebugInfo
            << "Restricting to cellZone(s) " << flatOutput(zoneNames_) << endl;
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
    surface_.clear();
    meshCells_.clear();
    isoSurfacePtr_.reset(nullptr);

    // Clear derived data
    sampledSurface::clearGeom();

    ignoreCellsPtr_.reset(nullptr);

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
    os  << "isoSurfaceTopo: " << name() << " :"
        << "  field:" << isoField_
        << "  value:" << isoVal_;
        //<< "  faces:" << faces().size()   // possibly no geom yet
        //<< "  points:" << points().size();
}


// ************************************************************************* //
