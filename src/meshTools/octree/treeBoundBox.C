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

Description

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "treeBoundBox.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

treeBoundBox treeBoundBox::greatBox
(
    vector(-GREAT, -GREAT, -GREAT),
    vector(GREAT, GREAT, GREAT)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as the bounding box of the given pointField
treeBoundBox::treeBoundBox(const pointField& points)
:
    boundBox()
{
    if (points.size() == 0)
    {
        WarningIn("treeBoundBox::treeBoundBox(const pointField& points)")
            << "cannot find bounding box for zero sized pointField"
            << "returning zero" << endl;

        return;
    }

    min() = points[0];
    max() = points[0];

    forAll(points, i)
    {
        min() = ::Foam::min(min(), points[i]);
        max() = ::Foam::max(max(), points[i]);
    }
}


// Construct from Istream
treeBoundBox::treeBoundBox(Istream& is)
:
    boundBox(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar treeBoundBox::minDim() const
{
    return ::Foam::min
    (
        max().x() - min().x(),
        ::Foam::min
        (
            max().y() - min().y(),
            max().z() - min().z()
        )
    );
}


scalar treeBoundBox::maxDim() const
{
    return ::Foam::max
    (
        max().x() - min().x(),
        ::Foam::max
        (
            max().y() - min().y(),
            max().z() - min().z()
        )
    );
}


scalar treeBoundBox::avgDim() const
{
    return
    (
        (max().x() - min().x()) +
        (max().y() - min().y()) +
        (max().z() - min().z())
    )/3.0;
}


scalar treeBoundBox::typDim() const
{
    return avgDim();
}


point treeBoundBox::mid() const
{
    return 0.5*(min() + max());
}


pointField treeBoundBox::points() const
{
    pointField points(8);
    label pointI = 0;

    points[pointI++] = min();
    points[pointI++] = point(min().x(), max().y(), min().z());
    points[pointI++] = point(max().x(), max().y(), min().z());
    points[pointI++] = point(max().x(), min().y(), min().z());

    points[pointI++] = point(min().x(), min().y(), max().z());
    points[pointI++] = point(min().x(), max().y(), max().z());
    points[pointI++] = max();
    points[pointI++] = point(max().x(), min().y(), max().z());

    return points;
}


edgeList treeBoundBox::edges() const
{
    edgeList edges(12);
    label edgeI = 0;

    // bottom face
    edges[edgeI++] = edge(0, 1);
    edges[edgeI++] = edge(1, 2);
    edges[edgeI++] = edge(2, 3);
    edges[edgeI++] = edge(3, 0);

    // top face
    edges[edgeI++] = edge(4, 5);
    edges[edgeI++] = edge(5, 6);
    edges[edgeI++] = edge(6, 7);
    edges[edgeI++] = edge(7, 4);

    // side edges
    edges[edgeI++] = edge(0, 4);
    edges[edgeI++] = edge(1, 5);
    edges[edgeI++] = edge(2, 6);
    edges[edgeI++] = edge(3, 7);

    return edges;
}


// Octant to bounding box
treeBoundBox treeBoundBox::subBbox(const direction octant) const
{
    if (octant > 7)
    {
        FatalErrorIn
        (
            "treeBoundBox::subCube(const direction)"
        )   << "octant should be [0..7]"
            << abort(FatalError);
    }

    scalar leftx, lefty, leftz;
    scalar rightx, righty, rightz;

    scalar midx=0.5*(min().x() + max().x());
    scalar midy=0.5*(min().y() + max().y());
    scalar midz=0.5*(min().z() + max().z());

    // X half
    if (octant & treeBoundBox::RIGHTHALF)
    {
        leftx = midx;
        rightx = max().x();
    }
    else
    {
        leftx = min().x();
        rightx = midx;
    }

    // Y half
    if (octant & treeBoundBox::TOPHALF)
    {
        lefty = midy;
        righty = max().y();
    }
    else
    {
        lefty = min().y();
        righty = midy;
    }

    // Z half
    if (octant & treeBoundBox::FRONTHALF)
    {
        leftz = midz;
        rightz = max().z();
    }
    else
    {
        leftz = min().z();
        rightz = midz;
    }

    return treeBoundBox
    (
        point(leftx, lefty, leftz),
        point(rightx, righty, rightz)
    );
}


// Octant to bounding box using permutation only.
treeBoundBox treeBoundBox::subBbox(const point& mid, const direction octant)
 const
{
    if (octant > 7)
    {
        FatalErrorIn
        (
            "treeBoundBox::subCube(const point&, const direction)"
        )   << "octant should be [0..7]"
            << abort(FatalError);
    }

    treeBoundBox subBb;
    point& subMin = subBb.min();
    point& subMax = subBb.max();

    if (octant & treeBoundBox::RIGHTHALF)
    {
        subMin.x() = mid.x();
        subMax.x() = max().x();
    }
    else
    {
        subMin.x() = min().x();
        subMax.x() = mid.x();
    }
    if (octant & treeBoundBox::TOPHALF)
    {
        subMin.y() = mid.y();
        subMax.y() = max().y();
    }
    else
    {
        subMin.y() = min().y();
        subMax.y() = mid.y();
    }
    if (octant & treeBoundBox::FRONTHALF)
    {
        subMin.z() = mid.z();
        subMax.z() = max().z();
    }
    else
    {
        subMin.z() = min().z();
        subMax.z() = mid.z();
    }

    return subBb;
}


// line intersection. Returns true if line (start to end) inside
// bb or intersects bb. Sets pt to intersection.
//
// Sutherlands algorithm:
//   loop
//     - start = intersection of line with one of the planes bounding
//       the bounding box
//     - stop if start inside bb (return true)
//     - stop if start and end in same 'half' (e.g. both above bb)
//       (return false)
//
// Uses posBits to efficiently determine 'half' in which start and end
// point are.
//
// Note:
//   - sets coordinate to exact position: e.g. pt.x() = min().x()
//     since plane intersect routine might have truncation error.
//     This makes sure that posBits tests 'inside'
bool treeBoundBox::intersects
(
    const point& start,
    const point& end,
    point& pt
) const
{
    vector vec(end - start);

    pt = start;

    const direction endBits = posBits(end);

    while(true)
    {
        direction ptBits = posBits(pt);

        if (ptBits == 0)
        {
            // pt inside bb
            return true;
        }

        if ((ptBits & endBits) != 0)
        {
            // pt and end in same block outside of bb
            return false;
        }

        if (ptBits & LEFTBIT)
        {
            // Intersect with plane V=min, n=-1,0,0
            if (Foam::mag(vec.x()) > VSMALL)
            {
                scalar s = (min().x() - pt.x())/vec.x();
                pt.x() = min().x();
                pt.y() = pt.y() + vec.y()*s;
                pt.z() = pt.z() + vec.z()*s;
            }
        }
        if (ptBits & RIGHTBIT)
        {
            // Intersect with plane V=max, n=1,0,0
            if (Foam::mag(vec.x()) > VSMALL)
            {
                scalar s = (max().x() - pt.x())/vec.x();
                pt.x() = max().x();
                pt.y() = pt.y() + vec.y()*s;
                pt.z() = pt.z() + vec.z()*s;
            }
        }

        if (ptBits & BELOWBIT)
        {
            // Intersect with plane V=min, n=0,-1,0
            if (Foam::mag(vec.y()) > VSMALL)
            {
                scalar s = (min().y() - pt.y())/vec.y();
                pt.x() = pt.x() + vec.x()*s;
                pt.y() = min().y();
                pt.z() = pt.z() + vec.z()*s;
            }
        }
        if (ptBits & ABOVEBIT)
        {
            // Intersect with plane V=max, n=0,1,0
            if (Foam::mag(vec.y()) > VSMALL)
            {
                scalar s = (max().y() - pt.y())/vec.y();
                pt.x() = pt.x() + vec.x()*s;
                pt.y() = max().y();
                pt.z() = pt.z() + vec.z()*s;
            }
        }

        if (ptBits & BEHINDBIT)
        {
            // Intersect with plane V=min, n=0,0,-1
            if (Foam::mag(vec.z()) > VSMALL)
            {
                scalar s = (min().z() - pt.z())/vec.z();
                pt.x() = pt.x() + vec.x()*s;
                pt.y() = pt.y() + vec.y()*s;
                pt.z() = min().z();
            }
        }
        if (ptBits & INFRONTBIT)
        {
            // Intersect with plane V=max, n=0,0,1
            if (Foam::mag(vec.z()) > VSMALL)
            {
                scalar s = (max().z() - pt.z())/vec.z();
                pt.x() = pt.x() + vec.x()*s;
                pt.y() = pt.y() + vec.y()*s;
                pt.z() = max().z();
            }
        }
    }
}


// this.bb fully contains bb
bool treeBoundBox::contains(const treeBoundBox& bb) const
{
    return contains(bb.min()) && contains(bb.max());
}


bool treeBoundBox::containsNarrow(const point& sample) const
{
    return
    (
        (sample.x() > min().x()) &&
        (sample.y() > min().y()) &&
        (sample.z() > min().z()) &&
        (sample.x() < max().x()) &&
        (sample.y() < max().y()) &&
        (sample.z() < max().z())
    );
}

bool treeBoundBox::contains(const vector& dir, const point& sample) const
{
    //
    // Compare all components against min and max of bb
    //

    for (direction cmpt=0; cmpt<3; cmpt++)
    {
        if (sample[cmpt] < min()[cmpt])
        {
            return false;
        }
        else if (sample[cmpt] == min()[cmpt])
        {
            // On edge. Outside if direction points outwards.
            if (dir[cmpt] < 0)
            {
                return false;
            }
        }

        if (sample[cmpt] > max()[cmpt])
        {
            return false;
        }
        else if (sample[cmpt] == max()[cmpt])
        {
            // On edge. Outside if direction points outwards.
            if (dir[cmpt] > 0)
            {
                return false;
            }
        }
    }

    // All components inside bb
    return true;
}


// Code position of point relative to box
direction treeBoundBox::posBits(const point& pt) const
{
    direction posBits = 0;

    if (pt.x() < min().x())
    {
        posBits |= LEFTBIT;
    }
    if (pt.x() > max().x())
    {
        posBits |= RIGHTBIT;
    }

    if (pt.y() < min().y())
    {
        posBits |= BELOWBIT;
    }
    if (pt.y() > max().y())
    {
        posBits |= ABOVEBIT;
    }

    if (pt.z() < min().z())
    {
        posBits |= BEHINDBIT;
    }
    if (pt.z() > max().z())
    {
        posBits |= INFRONTBIT;
    }
    return posBits;
}


// nearest and furthest corner coordinate.
// !names of treeBoundBox::min() and treeBoundBox::max() are confusing!
void treeBoundBox::calcExtremities
(
    const point& sample,
    point& nearest,
    point& furthest
) const
{
    scalar nearX, nearY, nearZ;
    scalar farX, farY, farZ;

    if (Foam::mag(min().x() - sample.x()) < Foam::mag(max().x() - sample.x()))
    {
        nearX = min().x();
        farX = max().x();
    }
    else
    {
        nearX = max().x();
        farX = min().x();
    }

    if (Foam::mag(min().y() - sample.y()) < Foam::mag(max().y() - sample.y()))
    {
        nearY = min().y();
        farY = max().y();
    }
    else
    {
        nearY = max().y();
        farY = min().y();
    }

    if (Foam::mag(min().z() - sample.z()) < Foam::mag(max().z() - sample.z()))
    {
        nearZ = min().z();
        farZ = max().z();
    }
    else
    {
        nearZ = max().z();
        farZ = min().z();
    }

    nearest = point(nearX, nearY, nearZ);
    furthest = point(farX, farY, farZ);
}


scalar treeBoundBox::maxDist(const point& sample) const
{
    point near, far;
    calcExtremities(sample, near, far);

    return Foam::mag(far - sample);
}


// Distance comparator
// Compare all vertices of bounding box against all of other bounding
// box to see if all vertices of one are nearer
label treeBoundBox::distanceCmp
(
    const point& sample,
    const treeBoundBox& other
) const
{
    //
    // Distance sample <-> nearest and furthest away vertex of this
    //

    point nearThis, farThis;

    // get nearest and furthest away vertex
    calcExtremities(sample, nearThis, farThis);

    const scalar minDistThis = 
        sqr(nearThis.x() - sample.x())
     +  sqr(nearThis.y() - sample.y())
     +  sqr(nearThis.z() - sample.z());
    const scalar maxDistThis = 
        sqr(farThis.x() - sample.x())
     +  sqr(farThis.y() - sample.y())
     +  sqr(farThis.z() - sample.z());

    //
    // Distance sample <-> other
    //

    point nearOther, farOther;

    // get nearest and furthest away vertex
    other.calcExtremities(sample, nearOther, farOther);

    const scalar minDistOther = 
        sqr(nearOther.x() - sample.x())
     +  sqr(nearOther.y() - sample.y())
     +  sqr(nearOther.z() - sample.z());
    const scalar maxDistOther = 
        sqr(farOther.x() - sample.x())
     +  sqr(farOther.y() - sample.y())
     +  sqr(farOther.z() - sample.z());

    //
    // Categorize
    //
    if (maxDistThis < minDistOther)
    {
        // All vertices of this are nearer to sample than any vertex of other
        return -1;
    }
    else if (minDistThis > maxDistOther)
    {
        // All vertices of this are further from sample than any vertex of other
        return 1;
    }
    else
    {
        // Mixed bag
        return 0;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool operator==(const treeBoundBox& a, const treeBoundBox& b)
{
    return (a.min() == b.min()) && (a.max() == b.max());
}


bool operator!=(const treeBoundBox& a, const treeBoundBox& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * IOstream Operator  * * * * * * * * * * * * * //

Istream& operator>>(Istream& is, treeBoundBox& bb)
{
    is >> bb.min() >> bb.max();
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
