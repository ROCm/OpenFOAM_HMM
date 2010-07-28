/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2.
inline bool Foam::pointEdgeStructuredWalk::update
(
    const pointEdgeStructuredWalk& w2,
    const scalar tol
)
{
    if (!valid())
    {
        // current not yet set. Walked from w2 to here (=point0)
        dist_ = w2.dist_ + mag(point0_-w2.previousPoint_);
        previousPoint_ = point0_;
        data_ = w2.data_;

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::pointEdgeStructuredWalk::pointEdgeStructuredWalk()
:
    point0_(vector::max),
    previousPoint_(vector::max),
    dist_(0),
    data_(vector::zero)
{}


// Construct from origin, distance
inline Foam::pointEdgeStructuredWalk::pointEdgeStructuredWalk
(
    const point& point0,
    const point& previousPoint,
    const scalar dist,
    const vector& data
)
:
    point0_(point0),
    previousPoint_(previousPoint),
    dist_(dist),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::pointEdgeStructuredWalk::inZone() const
{
    return point0_ != vector::max;
}


//inline const Foam::point& Foam::pointEdgeStructuredWalk::previousPoint() const
//{
//    return previousPoint_;
//}


inline Foam::scalar Foam::pointEdgeStructuredWalk::dist() const
{
    return dist_;
}


inline const Foam::vector& Foam::pointEdgeStructuredWalk::data() const
{
    return data_;
}


inline bool Foam::pointEdgeStructuredWalk::valid() const
{
    return previousPoint_ != vector::max;
}


// Checks for cyclic points
inline bool Foam::pointEdgeStructuredWalk::sameGeometry
(
    const pointEdgeStructuredWalk& w2,
    const scalar tol
) const
{
    scalar diff = Foam::mag(dist() - w2.dist());

    if (diff < SMALL)
    {
        return true;
    }
    else
    {
        if ((dist() > SMALL) && ((diff/dist()) < tol))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


inline void Foam::pointEdgeStructuredWalk::leaveDomain
(
    const polyPatch& patch,
    const label patchPointI,
    const point& coord
)
{
    previousPoint_ -= coord;
}


inline void Foam::pointEdgeStructuredWalk::transform(const tensor& rotTensor)
{
    previousPoint_ = Foam::transform(rotTensor, previousPoint_);
}


// Update absolute geometric quantities. Note that distance (dist_)
// is not affected by leaving/entering domain.
inline void Foam::pointEdgeStructuredWalk::enterDomain
(
    const polyPatch& patch,
    const label patchPointI,
    const point& coord
)
{
    // back to absolute form
    previousPoint_ += coord;
}


// Update this with information from connected edge
inline bool Foam::pointEdgeStructuredWalk::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const label edgeI,
    const pointEdgeStructuredWalk& edgeInfo,
    const scalar tol
)
{
    if (inZone())
    {
        return update(edgeInfo, tol);
    }
    else
    {
        return false;
    }
}


// Update this with new information on same point
inline bool Foam::pointEdgeStructuredWalk::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const pointEdgeStructuredWalk& newPointInfo,
    const scalar tol
)
{
    if (inZone())
    {
        return update(newPointInfo, tol);
    }
    else
    {
        return false;
    }
}


// Update this with new information on same point. No extra information.
inline bool Foam::pointEdgeStructuredWalk::updatePoint
(
    const pointEdgeStructuredWalk& newPointInfo,
    const scalar tol
)
{
    return update(newPointInfo, tol);
}


// Update this with information from connected point
inline bool Foam::pointEdgeStructuredWalk::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointI,
    const pointEdgeStructuredWalk& pointInfo,
    const scalar tol
)
{
    if (inZone())
    {
        return update(pointInfo, tol);
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointEdgeStructuredWalk::operator==
(
    const Foam::pointEdgeStructuredWalk& rhs
) const
{
    return previousPoint_ == rhs.previousPoint_;
}


inline bool Foam::pointEdgeStructuredWalk::operator!=
(
    const Foam::pointEdgeStructuredWalk& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
