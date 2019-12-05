/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "treeDataPoint.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "triangleFuncs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(treeDataPoint, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeDataPoint::treeDataPoint(const pointField& points)
:
    points_(points),
    useSubset_(false)
{}


Foam::treeDataPoint::treeDataPoint
(
    const pointField& points,
    const labelUList& pointLabels,
    const bool useSubsetPoints
)
:
    points_(points),
    pointLabels_(pointLabels),
    useSubset_(useSubsetPoints)
{}


Foam::treeDataPoint::treeDataPoint
(
    const pointField& points,
    labelList&& pointLabels,
    const bool useSubsetPoints
)
:
    points_(points),
    pointLabels_(std::move(pointLabels)),
    useSubset_(useSubsetPoints)
{}


Foam::treeDataPoint::findNearestOp::findNearestOp
(
    const indexedOctree<treeDataPoint>& tree
)
:
    tree_(tree)
{}


Foam::treeDataPoint::findIntersectOp::findIntersectOp
(
    const indexedOctree<treeDataPoint>& tree
)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::treeDataPoint::shapePoints() const
{
    if (useSubset_)
    {
        return pointField(points_, pointLabels_);
    }

    return points_;
}


Foam::volumeType Foam::treeDataPoint::getVolumeType
(
    const indexedOctree<treeDataPoint>& oc,
    const point& sample
) const
{
    return volumeType::UNKNOWN;
}


bool Foam::treeDataPoint::overlaps
(
    const label index,
    const treeBoundBox& cubeBb
) const
{
    return cubeBb.contains(shapePoint(index));
}


bool Foam::treeDataPoint::overlaps
(
    const label index,
    const point& centre,
    const scalar radiusSqr
) const
{
    return (magSqr(shapePoint(index) - centre) <= radiusSqr);
}


void Foam::treeDataPoint::findNearestOp::operator()
(
    const labelUList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    const treeDataPoint& shape = tree_.shapes();

    for (const label index : indices)
    {
        const point& pt = shape.shapePoint(index);

        const scalar distSqr = magSqr(pt - sample);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = pt;
        }
    }
}


void Foam::treeDataPoint::findNearestOp::operator()
(
    const labelUList& indices,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& minIndex,
    point& linePoint,
    point& nearestPoint
) const
{
    const treeDataPoint& shape = tree_.shapes();

    // Best so far
    scalar nearestDistSqr = GREAT;
    if (minIndex >= 0)
    {
        nearestDistSqr = magSqr(linePoint - nearestPoint);
    }

    for (const label index : indices)
    {
        const point& shapePt = shape.shapePoint(index);

        if (tightest.contains(shapePt))
        {
            // Nearest point on line
            pointHit pHit = ln.nearestDist(shapePt);
            const scalar distSqr = sqr(pHit.distance());

            if (distSqr < nearestDistSqr)
            {
                nearestDistSqr = distSqr;
                minIndex = index;
                linePoint = pHit.rawPoint();
                nearestPoint = shapePt;

                {
                    point& minPt = tightest.min();
                    minPt = min(ln.start(), ln.end());
                    minPt.x() -= pHit.distance();
                    minPt.y() -= pHit.distance();
                    minPt.z() -= pHit.distance();
                }
                {
                    point& maxPt = tightest.max();
                    maxPt = max(ln.start(), ln.end());
                    maxPt.x() += pHit.distance();
                    maxPt.y() += pHit.distance();
                    maxPt.z() += pHit.distance();
                }
            }
        }
    }
}


bool Foam::treeDataPoint::findIntersectOp::operator()
(
    const label index,
    const point& start,
    const point& end,
    point& result
) const
{
    NotImplemented;
    return false;
}


// ************************************************************************* //
