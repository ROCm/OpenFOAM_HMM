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

#include "error.H"
#include "treeBoundBox.H"
#include "point.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null setting points to zero
inline treeBoundBox::treeBoundBox()
:
    boundBox()
{}


// Construct from components
inline treeBoundBox::treeBoundBox(const point& min, const point& max)
:
    boundBox(min, max)
{}


// Construct from components
inline treeBoundBox::treeBoundBox(const boundBox& bb)
:
    boundBox(bb)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Returns octant in which sample resides. Reverse of subBbox.
inline direction treeBoundBox::subOctant(const point& sample) const
{
    scalar midx = 0.5*(max().x() + min().x());
    scalar midy = 0.5*(max().y() + min().y());
    scalar midz = 0.5*(max().z() + min().z());

    direction octant = 0;

    if (sample.x() > midx)
    {
        octant |= treeBoundBox::RIGHTHALF;
    }

    if (sample.y() > midy)
    {
        octant |= treeBoundBox::TOPHALF;
    }

    if (sample.z() > midz)
    {
        octant |= treeBoundBox::FRONTHALF;
    }

    return octant;
}


// Returns octant in which sample resides. Reverse of subBbox. Precalculated
// midpoint
inline direction treeBoundBox::subOctant
(
    const point& mid,
    const point& sample
)
{
    direction octant = 0;

    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }

    return octant;
}


// Returns octant in which sample resides. Reverse of subBbox. Flags sample
// exactly on edge.
inline direction treeBoundBox::subOctant(const point& sample, bool& onEdge)
 const
{
    point mid = 0.5*(max() + min());

    direction octant = 0;

    onEdge = false;
    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }
    else if (sample.x() == mid.x())
    {
        onEdge = true;
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }
    else if (sample.y() == mid.y())
    {
        onEdge = true;
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }
    else if (sample.z() == mid.z())
    {
        onEdge = true;
    }

    return octant;
}


// Returns octant in which sample resides. Reverse of subBbox. Precalculated
// midpoint
inline direction treeBoundBox::subOctant
(
    const point& mid,
    const point& sample,
    bool& onEdge
)
{
    direction octant = 0;

    onEdge = false;
    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }
    else if (sample.x() == mid.x())
    {
        onEdge = true;
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }
    else if (sample.y() == mid.y())
    {
        onEdge = true;
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }
    else if (sample.z() == mid.z())
    {
        onEdge = true;
    }

    return octant;
}


// Returns octant in which intersection resides.
// Precalculated midpoint. If the sample is on the dividing line between
// the octants the direction vector determines which octant to use
// (i.e. in which octant the sample would be if it were moved along dir)
inline direction treeBoundBox::subOctant
(
    const point& mid,
    const vector& dir,
    const point& sample,
    bool& onEdge
)
{
    direction octant = 0;

    onEdge = false;

    if (sample.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }
    else if (sample.x() == mid.x())
    {
        onEdge = true;

        if (dir.x() > 0)
        {
            octant |= treeBoundBox::RIGHTHALF;
        }
    }

    if (sample.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }
    else if (sample.y() == mid.y())
    {
        onEdge = true;

        if (dir.y() > 0)
        {
            octant |= treeBoundBox::TOPHALF;
        }
    }

    if (sample.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }
    else if (sample.z() == mid.z())
    {
        onEdge = true;

        if (dir.z() > 0)
        {
            octant |= treeBoundBox::FRONTHALF;
        }
    }

    return octant;
}


// true if bb's intersect or overlap.
// Note: <= to make sure we catch all.
inline bool treeBoundBox::intersects(const treeBoundBox& bb) const
{
    return boundBox::intersects(bb);
}


inline bool treeBoundBox::contains(const point& sample) const
{
    return
    (
        (sample.x() >= min().x()) &&
        (sample.y() >= min().y()) &&
        (sample.z() >= min().z()) &&
        (sample.x() <= max().x()) &&
        (sample.y() <= max().y()) &&
        (sample.z() <= max().z())
    );
}


//- Return slightly wider bounding box
inline treeBoundBox treeBoundBox::extend(Random& rndGen, const scalar s) const
{
    treeBoundBox bb(*this);

    const vector span(bb.max() - bb.min());

    bb.min() -= cmptMultiply(s*rndGen.vector01(), span);
    bb.max() += cmptMultiply(s*rndGen.vector01(), span);

    return bb;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
