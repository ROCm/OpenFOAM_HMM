/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "searchableBox.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableBox, 0);
addToRunTimeSelectionTable(searchableSurface, searchableBox, dict);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableBox::searchableBox
(
    const word& name,
    const treeBoundBox& bb
)
:
    searchableSurface(name),
    treeBoundBox(bb)
{}


Foam::searchableBox::searchableBox
(
    const word& name,
    const objectRegistry& obj,
    const dictionary& dict
)
:
    searchableSurface(name),
    treeBoundBox(dict.lookup("min"), dict.lookup("max"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableBox::~searchableBox()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableBox::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    // Point can be inside or outside. For every component direction can be
    // left of min, right of max or inbetween.
    // - outside points: project first one x plane (either min().x()
    // or max().x()), then onto y plane and finally z. You should be left
    // with intersection point
    // - inside point: find nearest side (compare to mid point). Pick any
    // one of three points.

    const point bbMid(mid());

    // Outside point projected onto cube
    point interPt(sample);
    bool outside = false;

    // (for internal points) Per direction what nearest cube side is
    point near;

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (interPt[dir] < min()[dir])
        {
            interPt[dir] = min()[dir];
            outside = true;
        }
        else if (interPt[dir] > max()[dir])
        {
            interPt[dir] = max()[dir];
            outside = true;
        }
        else if (interPt[dir] > bbMid[dir])
        {
            near[dir] = max()[dir];
        }
        else
        {
            near[dir] = min()[dir];
        }
    }


    // For outside points the interPt will be correct now. Handle inside points
    // using the three near distances. Project onto the nearest plane.
    if (!outside)
    {
        vector dist(cmptMag(interPt - near));

        if (dist.x() < dist.y())
        {
            if (dist.x() < dist.z())
            {
                interPt.x() = near.x();
            }
            else
            {
                interPt.z() = near.z();
            }
        }
        else
        {
            if (dist.y() < dist.z())
            {
                interPt.y() = near.y();
            }
            else
            {
                interPt.z() = near.z();
            }
        }
    }

    return pointIndexHit(true, interPt, 0);
}


Foam::pointIndexHit Foam::searchableBox::findNearestOnEdge
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    const point bbMid(mid());

    point interPt(sample);

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        // Project onto left or right depending on mid
        if (interPt[dir] > bbMid[dir])
        {
            interPt[dir] = max()[dir];
        }
        else
        {
            interPt[dir] = min()[dir];
        }
    }

    return pointIndexHit(true, interPt, 0);
}


Foam::pointIndexHit Foam::searchableBox::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint
) const
{
    notImplemented
    (
        "searchableBox::findNearest"
        "(const linePointRef&, treeBoundBox&, point&)"
    );
    return pointIndexHit();
}


Foam::pointIndexHit Foam::searchableBox::findLine
(
    const point& start,
    const point& end
) const
{
    point intPt;
    bool foundInter = intersects(start, end, intPt);

    return pointIndexHit(foundInter, intPt, 0);
}


Foam::pointIndexHit Foam::searchableBox::findLineAny
(
    const point& start,
    const point& end
) const
{
    return findLine(start, end);
}


Foam::searchableSurface::volumeType Foam::searchableBox::getVolumeType
(
    const point& pt
) const
{
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (pt[dir] < min()[dir] || pt[dir] > max()[dir])
        {
            return OUTSIDE;
        }
    }

    return INSIDE;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
