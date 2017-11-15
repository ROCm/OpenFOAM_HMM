/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "surfMeshPlaneSampler.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"
#include "coordinateSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshPlaneSampler, 0);
    addNamedToRunTimeSelectionTable
    (
        surfMeshSampler,
        surfMeshPlaneSampler,
        word,
        plane
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMeshPlaneSampler::transferContent()
{
    SurfaceSource& src = static_cast<SurfaceSource&>(*this);
    surfMesh& dst = getOrCreateSurfMesh();

    dst.transfer(src);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshPlaneSampler::surfMeshPlaneSampler
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const keyType& zoneKey,
    const bool triangulate
)
:
    surfMeshSampler(name, mesh),
    SurfaceSource(planeDesc),
    zoneKey_(zoneKey),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) == -1)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
}


Foam::surfMeshPlaneSampler::surfMeshPlaneSampler
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    surfMeshSampler(name, mesh, dict),
    SurfaceSource(plane(dict)),
    zoneKey_(keyType::null),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    // Make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found("coordinateSystem"))
    {
        coordinateSystem cs(mesh, dict.subDict("coordinateSystem"));

        point  base = cs.globalPosition(planeDesc().refPoint());
        vector norm = cs.globalVector(planeDesc().normal());

        // Assign the plane description
        static_cast<plane&>(*this) = plane(base, norm);
    }

    dict.readIfPresent("zone", zoneKey_);

    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) == -1)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMeshPlaneSampler::~surfMeshPlaneSampler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfMeshPlaneSampler::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::surfMeshPlaneSampler::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::surfMeshPlaneSampler::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    labelList selectedCells = mesh().cellZones().findMatching(zoneKey_).used();
    if (selectedCells.empty())
    {
        reCut(mesh(), triangulate_);
    }
    else
    {
        Foam::sort(selectedCells);
        reCut(mesh(), triangulate_, selectedCells);
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    transferContent();

    needsUpdate_ = false;
    return true;
}


bool Foam::surfMeshPlaneSampler::sample
(
    const word& fieldName
) const
{
    return
    (
        sampleType<scalar>(fieldName)
     || sampleType<vector>(fieldName)
     || sampleType<sphericalTensor>(fieldName)
     || sampleType<symmTensor>(fieldName)
     || sampleType<tensor>(fieldName)
    );
}


void Foam::surfMeshPlaneSampler::print(Ostream& os) const
{
    os  << "surfMeshPlaneSampler: " << name() << " :"
        << "  base:" << cuttingPlane::refPoint()
        << "  normal:" << cuttingPlane::normal()
        << "  triangulate:" << triangulate_
        << "  faces:"  << SurfaceSource::surfFaces().size()
        << "  points:" << SurfaceSource::points().size();
}


// ************************************************************************* //
