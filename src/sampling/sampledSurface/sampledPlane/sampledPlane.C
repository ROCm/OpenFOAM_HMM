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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const keyType& zoneKey,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    cuttingPlane(planeDesc),
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


Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    cuttingPlane(plane(dict)),
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
        mesh().cellZones().selection(zoneKey_).sortedToc()
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
        << "  base:" << refPoint()
        << "  normal:" << normal()
        << "  triangulate:" << triangulate_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
