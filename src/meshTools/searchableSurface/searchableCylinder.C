/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "searchableCylinder.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableCylinder, 0);
addToRunTimeSelectionTable(searchableSurface, searchableCylinder, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableCylinder::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    vector v(sample - point1_);

    // Decompose sample-point1 into normal and parallel component
    scalar parallel = (v & unitDir_);

    // Remove the parallel component
    v -= parallel*unitDir_;
    scalar magV = mag(v);

    if (parallel <= 0)
    {
        // nearest is at point1 end cap. Clip to radius.
        if (magV < ROOTVSMALL)
        {
            info.setPoint(point1_);
        }
        else
        {
            //info.setPoint(point1_ + min(magV, radius_)*v/magV);
            info.setPoint(point1_ + radius_*(v/magV));
        }
    }
    else if (parallel >= magDir_)
    {
        // nearest is at point2
        if (magV < ROOTVSMALL)
        {
            info.setPoint(point2_);
        }
        else
        {
            info.setPoint(point2_ + min(magV, radius_)*v/magV);
        }
    }
    else
    {
        if (magV < ROOTVSMALL)
        {
            info.setPoint(point1_ + parallel*unitDir_);
        }
        else
        {
            // Project onto radius
            info.setPoint(point1_ + parallel*unitDir_ + radius_*v/magV);
        }
    }

    if (magSqr(sample - info.rawPoint()) < nearestDistSqr)
    {
        info.setHit();
        info.setIndex(0);
    }

    return info;
}


// From http://www.gamedev.net/community/forums/topic.asp?topic_id=467789 -
// intersection of cylinder with ray
void Foam::searchableCylinder::findLineAll
(
    const point& start,
    const point& end,
    pointIndexHit& near,
    pointIndexHit& far
) const
{
    near.setMiss();
    far.setMiss();

    // Line as P = start + t*V
    const vector V(end-start);

//Pout<< "point1:" << point1_ << " point2:" << point2_
//    << " start:" << start << " end:" << end << endl;

    const vector x = (start-point1_) ^ unitDir_;
    const vector y = V ^ unitDir_;
    const scalar d = sqr(radius_);

    // Second order equation of the form a*t^2 + b*t + c
    const scalar a = (y&y);
    const scalar b = 2*(x&y);
    const scalar c = (x&x)-d;

    const scalar disc = b*b-4*a*c;

//Pout<< "a:" << a << " b:" << b << " c:" << c << " disc:" << disc
//    << endl;

    if (disc < 0)
    {
        return;
    }
    else if (disc < ROOTVSMALL)
    {
        // Single solution
        if (mag(a) > ROOTVSMALL)
        {
            scalar t = -b/(2*a);
            if (t >= 0 && t <= 1)
            {
                near.setPoint(start + t*V);
                near.setHit();
                near.setIndex(0);
            }
        }
    }
    else
    {
        if (mag(a) > ROOTVSMALL)
        {
            scalar sqrtDisc = sqrt(disc);

            scalar t1 = (-b + sqrtDisc)/2*a;
            scalar t2 = (-b - sqrtDisc)/2*a;

            if (t1 < t2)
            {
                if (t1 >= 0 && t1 <= 1)
                {
                    near.setPoint(start + t1*V);
                    near.setHit();
                    near.setIndex(0);
                }
                if (t2 >= 0 && t2 <= 1)
                {
                    far.setPoint(start + t2*V);
                    far.setHit();
                    far.setIndex(0);
                }
            }
            else
            {
                if (t2 >= 0 && t2 <= 1)
                {
                    near.setPoint(start + t2*V);
                    near.setHit();
                    near.setIndex(0);
                }
                if (t1 >= 0 && t1 <= 1)
                {
                    far.setPoint(start + t1*V);
                    far.setHit();
                    far.setIndex(0);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableCylinder::searchableCylinder
(
    const IOobject& io,
    const point& point1,
    const point& point2,
    const scalar radius
)
:
    searchableSurface(io),
    point1_(point1),
    point2_(point2),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_),
    radius_(radius)
{
    Pout<< "point1_:" << point1_ << endl;
    Pout<< "point2_:" << point2_ << endl;
    Pout<< "magDir_:" << magDir_ << endl;
    Pout<< "unitDir_:" << unitDir_ << endl;
    Pout<< "radius_:" << radius_ << endl;
}


Foam::searchableCylinder::searchableCylinder
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    point1_(dict.lookup("point1")),
    point2_(dict.lookup("point2")),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_),
    radius_(readScalar(dict.lookup("radius")))
{
    Pout<< "point1_:" << point1_ << endl;
    Pout<< "point2_:" << point2_ << endl;
    Pout<< "magDir_:" << magDir_ << endl;
    Pout<< "unitDir_:" << unitDir_ << endl;
    Pout<< "radius_:" << radius_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableCylinder::~searchableCylinder()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableCylinder::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


void Foam::searchableCylinder::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableCylinder::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected take second one.
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableCylinder::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Discard far intersection
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableCylinder::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit> >& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;
        findLineAll(start[i], end[i], near, far);

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].setSize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].setSize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].setSize(1);
                info[i][0] = far;
            }
            else
            {
                info[i].clear();
            }
        }
    }
}


void Foam::searchableCylinder::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableCylinder::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = vector::zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            vector v(info[i].hitPoint() - point1_);

            // Decompose sample-point1 into normal and parallel component
            scalar parallel = v & unitDir_;

            if (parallel < 0)
            {
                normal[i] = -unitDir_;
            }
            else if (parallel > magDir_)
            {
                normal[i] = -unitDir_;
            }
            else
            {
                // Remove the parallel component
                v -= parallel*unitDir_;
                normal[i] = v/mag(v);
            }
        }
    }
}


void Foam::searchableCylinder::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());
    volType = INSIDE;

    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        vector v(pt - point1_);

        // Decompose sample-point1 into normal and parallel component
        scalar parallel = v & unitDir_;

        if (parallel < 0)
        {
            // left of point1 endcap
            volType[pointI] = OUTSIDE;
        }
        else if (parallel > magDir_)
        {
            // right of point2 endcap
            volType[pointI] = OUTSIDE;
        }
        else
        {
            // Remove the parallel component
            v -= parallel*unitDir_;

            if (mag(v) > radius_)
            {
                volType[pointI] = OUTSIDE;
            }
            else
            {
                volType[pointI] = INSIDE;
            }
        }
    }
}


// ************************************************************************* //
