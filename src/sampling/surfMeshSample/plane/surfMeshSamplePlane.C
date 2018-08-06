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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMeshSamplePlane::checkBoundsIntersection
(
    const plane& pln,
    const boundBox& meshBb
) const
{
    // Verify specified bounding box
    if (!bounds_.empty())
    {
        // Bounding box does not overlap with (global) mesh!
        if (!bounds_.overlaps(meshBb))
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
    if (!meshBb.intersects(pln))
    {
        WarningInFunction
            << nl
            << name() << " : "
            << "Plane "<< pln << " does not intersect the mesh bounds "
            << meshBb
            << nl << endl;
    }
}


Foam::bitSet Foam::surfMeshSamplePlane::cellSelection
(
    const bool warnIntersect
) const
{
    // Zones requested and in use?
    const bool hasZones =
        returnReduce
        (
            (-1 != mesh().cellZones().findIndex(zoneNames_)),
            andOp<bool>()
        );


    bitSet cellsToSelect;

    if (hasZones)
    {
        cellsToSelect = mesh().cellZones().selection(zoneNames_);
    }

    // Subset the zoned cells with the bounds_.
    // For a full mesh, use the bounds to define the cell selection.

    // If there are zones cells, use them to build the effective mesh
    // bound box.
    // Note that for convenience we use cell centres here instead of the
    // cell points, since it will only be used for checking.

    boundBox meshBb;

    if (bounds_.empty())
    {
        // No bounds restriction, but will need the effective mesh boundBox
        // for checking intersections

        if (hasZones && warnIntersect)
        {
            const auto& cellCentres = static_cast<const fvMesh&>(mesh()).C();

            for (const label celli : cellsToSelect)
            {
                const point& cc = cellCentres[celli];

                meshBb.add(cc);
            }

            meshBb.reduce();
        }
        else
        {
            meshBb = mesh().bounds();  // use the regular mesh bounding box
        }
    }
    else
    {
        const auto& cellCentres = static_cast<const fvMesh&>(mesh()).C();

        // Only examine cells already set
        if (hasZones)
        {
            for (const label celli : cellsToSelect)
            {
                const point& cc = cellCentres[celli];

                meshBb.add(cc);

                if (!bounds_.contains(cc))
                {
                    cellsToSelect.unset(celli);
                }
            }

            meshBb.reduce();
        }
        else
        {
            const label len = mesh().nCells();

            cellsToSelect.resize(len);

            for (label celli=0; celli < len; ++celli)
            {
                const point& cc = cellCentres[celli];

                if (bounds_.contains(cc))
                {
                    cellsToSelect.set(celli);
                }
            }

            meshBb = mesh().bounds();  // use the regular mesh bounding box
        }
    }

    if (warnIntersect)
    {
        checkBoundsIntersection(*this, meshBb);
    }

    return cellsToSelect;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSamplePlane::surfMeshSamplePlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const wordRes& zones,
    const bool triangulate
)
:
    surfMeshSample(name, mesh),
    SurfaceSource(planeDesc),
    zoneNames_(zones),
    bounds_(),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug)
    {
        if (!zoneNames_.empty())
        {
            Info<< "cellZones " << flatOutput(zoneNames_);

            if (-1 == mesh.cellZones().findIndex(zoneNames_))
            {
                Info<< " not found!";
            }
            Info<< endl;
        }
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
    zoneNames_(),
    bounds_(dict.lookupOrDefault("bounds", boundBox::invertedBox)),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    if (!dict.readIfPresent("zones", zoneNames_) && dict.found("zone"))
    {
        zoneNames_.resize(1);
        dict.readEntry("zone", zoneNames_.first());
    }


    // Make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found("coordinateSystem"))
    {
        coordinateSystem cs(mesh, dict.subDict("coordinateSystem"));

        const point  orig = cs.globalPosition(planeDesc().origin());
        const vector norm = cs.globalVector(planeDesc().normal());

        if (debug)
        {
            Info<< "plane " << name << " :"
                << " origin:" << origin()
                << " normal:" << normal()
                << " defined within a local coordinateSystem" << endl;
        }

        // Assign the plane description
        static_cast<plane&>(*this) = plane(orig, norm);
    }

    if (debug)
    {
        Info<< "plane " << name << " :"
            << " origin:" << origin()
            << " normal:" << normal();

        if (!bounds_.empty())
        {
            Info<< " bounds:" << bounds_;
        }

        if (!zoneNames_.empty())
        {
            Info<< " cellZones " << flatOutput(zoneNames_);

            if (-1 == mesh.cellZones().findIndex(zoneNames_))
            {
                Info<< " not found!";
            }
        }
        Info<< endl;
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

    performCut(mesh(), triangulate_, this->cellSelection(true));

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
        << " base:" << plane::origin()
        << " normal:" << plane::normal()
        << " triangulate:" << triangulate_
        << " faces:"  << SurfaceSource::surfFaces().size()
        << " points:" << SurfaceSource::points().size();
}


// ************************************************************************* //
