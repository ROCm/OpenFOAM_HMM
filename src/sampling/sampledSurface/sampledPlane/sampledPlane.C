/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "sampledPlane.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"
#include "coordinateSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPlane, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledPlane,
        word,
        plane
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledPlane::checkBoundsIntersection
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
                << " do not overlap the mesh bounding box " << meshBb
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


Foam::bitSet Foam::sampledPlane::cellSelection
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

Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const wordRes& zones,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    cuttingPlane(planeDesc),
    zoneNames_(zones),
    bounds_(),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug)
    {
        if (!zoneNames_.empty())
        {
            Info<< " cellZones " << flatOutput(zoneNames_);

            if (-1 == mesh.cellZones().findIndex(zoneNames_))
            {
                Info<< " not found!";
            }
            Info<< endl;
        }
    }
}


Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    cuttingPlane(plane(dict)),
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

bool Foam::sampledPlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledPlane::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledPlane::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    performCut(mesh(), triangulate_, this->cellSelection(true));

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledPlane::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledPlane::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledPlane::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPlane::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledPlane::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledPlane::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledPlane::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledPlane::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPlane::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledPlane::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledPlane::print(Ostream& os) const
{
    os  << "sampledPlane: " << name() << " :"
        << " origin:" << plane::origin()
        << " normal:" << plane::normal()
        << " triangulate:" << triangulate_
        << " faces:" << faces().size()
        << " points:" << points().size();
}


// ************************************************************************* //
