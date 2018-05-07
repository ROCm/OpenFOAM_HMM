/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "surfMeshSamplePlane.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"
#include "coordinateSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshSamplePlane, 0);
    addNamedToRunTimeSelectionTable
    (
        surfMeshSample,
        surfMeshSamplePlane,
        word,
        plane
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSamplePlane::surfMeshSamplePlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const keyType& zoneKey,
    const bool triangulate
)
:
    surfMeshSample(name, mesh),
    SurfaceSource(planeDesc),
    zoneKey_(zoneKey),
    bounds_(),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) == -1)
    {
        Info<< "cellZone(s) " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
}


Foam::surfMeshSamplePlane::surfMeshSamplePlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    surfMeshSample(name, mesh, dict),
    SurfaceSource(plane(dict)),
    zoneKey_(dict.lookupOrDefault<keyType>("zone", keyType::null)),
    bounds_(dict.lookupOrDefault("bounds", boundBox::invertedBox)),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    // Make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found("coordinateSystem"))
    {
        coordinateSystem cs(mesh, dict.subDict("coordinateSystem"));

        const point  base = cs.globalPosition(planeDesc().refPoint());
        const vector norm = cs.globalVector(planeDesc().normal());

        // Assign the plane description
        static_cast<plane&>(*this) = plane(base, norm);
    }

    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) == -1)
    {
        Info<< "cellZone(s) " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfMeshSamplePlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::surfMeshSamplePlane::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::surfMeshSamplePlane::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    const plane& pln = static_cast<const plane&>(*this);

    // Verify specified bounding box
    if (!bounds_.empty())
    {
        // Bounding box does not overlap with (global) mesh!
        if (!bounds_.overlaps(mesh().bounds()))
        {
            WarningInFunction
                << nl
                << name() << " : "
                << "Bounds " << bounds_
                << " do not overlap the mesh bounding box " << mesh().bounds()
                << nl << endl;
        }

        // Plane does not intersect the bounding box
        if (!bounds_.intersects(pln))
        {
            WarningInFunction
                << nl
                << name() << " : "
                << "Plane "<< pln << " does not intersect the bounds "
                << bounds_
                << nl << endl;
        }
    }

    // Plane does not intersect the (global) mesh!
    if (!mesh().bounds().intersects(pln))
    {
        WarningInFunction
            << nl
            << name() << " : "
            << "Plane "<< pln << " does not intersect the mesh bounds "
            << mesh().bounds()
            << nl << endl;
    }

    labelList selectedCells
    (
        mesh().cellZones().findMatching(zoneKey_).sortedToc()
    );

    bool fullMesh = returnReduce(selectedCells.empty(), andOp<bool>());
    if (!bounds_.empty())
    {
        const auto& cellCentres = static_cast<const fvMesh&>(mesh()).C();

        if (fullMesh)
        {
            const label len = mesh().nCells();

            selectedCells.setSize(len);

            label count = 0;
            for (label celli=0; celli < len; ++celli)
            {
                if (bounds_.contains(cellCentres[celli]))
                {
                    selectedCells[count++] = celli;
                }
            }

            selectedCells.setSize(count);
        }
        else
        {
            label count = 0;
            for (const label celli : selectedCells)
            {
                if (bounds_.contains(cellCentres[celli]))
                {
                    selectedCells[count++] = celli;
                }
            }

            selectedCells.setSize(count);
        }

        fullMesh = false;
    }

    if (fullMesh)
    {
        reCut(mesh(), triangulate_);
    }
    else
    {
        reCut(mesh(), triangulate_, selectedCells);
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    // Transfer content
    getOrCreateSurfMesh().transfer
    (
        static_cast<SurfaceSource&>(*this)
    );

    needsUpdate_ = false;
    return true;
}


bool Foam::surfMeshSamplePlane::sample
(
    const word& fieldName,
    const word& sampleScheme
) const
{
    return
    (
        sampleType<scalar>(fieldName, sampleScheme)
     || sampleType<vector>(fieldName, sampleScheme)
     || sampleType<sphericalTensor>(fieldName, sampleScheme)
     || sampleType<symmTensor>(fieldName, sampleScheme)
     || sampleType<tensor>(fieldName, sampleScheme)
    );
}


void Foam::surfMeshSamplePlane::print(Ostream& os) const
{
    os  << "surfMeshSamplePlane: " << name() << " :"
        << "  base:" << cuttingPlane::refPoint()
        << "  normal:" << cuttingPlane::normal()
        << "  triangulate:" << triangulate_
        << "  faces:"  << SurfaceSource::surfFaces().size()
        << "  points:" << SurfaceSource::points().size();
}


// ************************************************************************* //
