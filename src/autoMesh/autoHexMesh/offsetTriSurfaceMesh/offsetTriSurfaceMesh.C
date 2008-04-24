/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "offsetTriSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(offsetTriSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, offsetTriSurfaceMesh, dict);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::offsetTriSurfaceMesh::offsetTriSurfaceMesh
(
    const IOobject& io,
    const triSurface& s,
    const scalar offset)
:
    triSurfaceMesh(io, s),
    offset_(offset)
{}


Foam::offsetTriSurfaceMesh::offsetTriSurfaceMesh
(
    const IOobject& io,
    const scalar offset
)
:
    triSurfaceMesh(io),
    offset_(offset)
{}


Foam::offsetTriSurfaceMesh::offsetTriSurfaceMesh
(
    const word& name,
    const objectRegistry& obj,
    const dictionary& dict
)
:
    triSurfaceMesh(name, obj, dict),
    offset_(readScalar(dict.lookup("offset")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::offsetTriSurfaceMesh::~offsetTriSurfaceMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointIndexHit Foam::offsetTriSurfaceMesh::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    // Find nearest (add offset to search span)
    pointIndexHit surfNearest = triSurfaceMesh::findNearest
    (
        sample,
        nearestDistSqr + Foam::sqr(offset_)
    );

    // Shift back onto surface
    if (surfNearest.hit())
    {
        vector n(sample-surfNearest.hitPoint());
        n /= mag(n)+VSMALL;
        surfNearest.setPoint(surfNearest.hitPoint() + offset_*n);
    }
    return surfNearest;
}


Foam::pointIndexHit Foam::offsetTriSurfaceMesh::findNearestOnEdge
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    // Find nearest (add offset to search span)
    pointIndexHit surfNearest = triSurfaceMesh::findNearestOnEdge
    (
        sample,
        nearestDistSqr + Foam::sqr(offset_)
    );

    // Shift back onto surface
    if (surfNearest.hit())
    {
        vector n = sample-surfNearest.hitPoint();
        n /= mag(n)+VSMALL;
        surfNearest.setPoint(surfNearest.hitPoint() + offset_*n);
    }
    return surfNearest;
}


Foam::searchableSurface::volumeType Foam::offsetTriSurfaceMesh::getVolumeType
(
    const point& sample
) const
{
    // Find the nearest point on background surface
    pointIndexHit surfNearest = triSurfaceMesh::findNearest
    (
        sample,
        Foam::sqr(GREAT)
    );

    if (!surfNearest.hit())
    {
        FatalErrorIn("offsetTriSurfaceMesh::getVolumeType(const point&)")
            << "treeBb:" << tree().bb()
            << " sample:" << sample
            << " surfNearest:" << surfNearest
            << abort(FatalError);
    }

    // Offset sample to the point.
    vector n(surfNearest.hitPoint()-sample);
    n /= mag(n)+VSMALL;

    triSurfaceTools::sideType t = triSurfaceTools::surfaceSide
    (
        *this,
        sample+offset_*n,
        surfNearest.index(),
        surfNearest.hitPoint()
    );

    if (t == triSurfaceTools::UNKNOWN)
    {
        return searchableSurface::UNKNOWN;
    }
    else if (t == triSurfaceTools::INSIDE)
    {
        return searchableSurface::INSIDE;
    }
    else if (t == triSurfaceTools::OUTSIDE)
    {
        return searchableSurface::OUTSIDE;
    }
    else
    {
        FatalErrorIn("offsetTriSurfaceMesh::getVolumeType(const point&)")
            << "problem" << abort(FatalError);
        return searchableSurface::UNKNOWN;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
