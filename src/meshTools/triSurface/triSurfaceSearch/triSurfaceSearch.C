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

#include "triSurfaceSearch.H"
#include "octree.H"
#include "boolList.H"
#include "octreeDataTriSurface.H"
#include "octreeDataPoint.H"
#include "triSurface.H"
#include "line.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const point triSurfaceSearch::greatPoint(GREAT, GREAT, GREAT);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface. Holds reference!
triSurfaceSearch::triSurfaceSearch(const triSurface& surface)
:
    surface_(surface),
    treePtr_(NULL)
{
    // Wrap surface information into helper object
    octreeDataTriSurface shapes(surface_);

    treeBoundBox treeBb(surface_.localPoints());

    scalar tol = 1E-6*treeBb.avgDim();

    point& bbMin = treeBb.min();
    bbMin.x() -= tol;
    bbMin.y() -= tol;
    bbMin.z() -= tol;

    point& bbMax = treeBb.max();
    bbMax.x() += tol;
    bbMax.y() += tol;
    bbMax.z() += tol;

    treePtr_ = new octree<octreeDataTriSurface>
    (
        treeBb,
        shapes,
        1,
        100.0,
        10.0
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triSurfaceSearch::~triSurfaceSearch()
{
    if (treePtr_)
    {
        delete treePtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Determine inside/outside for samples
boolList triSurfaceSearch::calcInside(const pointField& samples)
    const
{
    boolList inside(samples.size());

    forAll(samples, sampleI)
    {
        const point& sample = samples[sampleI];

        if (!tree().octreeBb().contains(sample))
        {
            inside[sampleI] = false;
        }
        else if
        (
            tree().getSampleType(sample)
         == octree<octreeDataTriSurface>::INSIDE
        )
        {
            inside[sampleI] = true;
        }
        else
        {
            inside[sampleI] = false;
        }
    }
    return inside;
}


labelList triSurfaceSearch::calcNearestTri
(
    const pointField& samples,
    const vector& span
) const
{
    labelList nearest(samples.size());

    forAll(samples, sampleI)
    {
        const point& sample = samples[sampleI];

        treeBoundBox tightest(sample, sample);
        tightest.min() -= span;
        tightest.max() += span;

        scalar tightestDist = Foam::GREAT;

        nearest[sampleI] = tree().findNearest
        (
            sample,
            tightest,
            tightestDist
        );
    }


    // Bit ropy - comparison of coordinate but is just check.
    if (span == greatPoint)
    {
        forAll(nearest, sampleI)
        {
            if (nearest[sampleI] == -1)
            {
                FatalErrorIn
                (
                    "triSurfaceSearch::calcNearestTri(const pointField&)"
                )   << "Did not find point " << samples[sampleI]
                    << " in octree spanning "
                    << tree().octreeBb() << abort(FatalError);
            }
        }
    }

    return nearest;
}


// Nearest point on surface
tmp<pointField> triSurfaceSearch::calcNearest
(
    const pointField& samples,
    const vector& span
) const
{
    const pointField& points = surface_.points();

    tmp<pointField> tnearest(new pointField(samples.size()));
    pointField& nearest = tnearest();

    labelList nearestTri(calcNearestTri(samples, span));

    forAll(nearestTri, sampleI)
    {
        label triI = nearestTri[sampleI];

        if (triI == -1)
        {
            nearest[sampleI] = greatPoint;
        }
        else
        {
            // Unfortunately octree does not return nearest point
            // so we have to recalculate it.
            const labelledTri& f = surface_[triI];

            triPointRef tri
            (
                points[f[0]],
                points[f[1]],
                points[f[2]]
            );

            nearest[sampleI] = tri.nearestPoint(samples[sampleI]).rawPoint();
        }
    }

    return tnearest;
}


pointIndexHit triSurfaceSearch::nearest(const point& pt, const vector& span)
 const
{
    pointIndexHit pInter;

    treeBoundBox tightest(pt, pt);
    tightest.min() -= span;
    tightest.max() += span;

    scalar tightestDist = Foam::GREAT;

    label triI = tree().findNearest(pt, tightest, tightestDist);

    if (triI == -1)
    {
        pInter.setMiss();
    }
    else
    {
        pInter.setHit();
        pInter.setIndex(triI);

        const labelledTri& f = surface_[triI];

        const pointField& points = surface_.points();

        triPointRef tri
        (
            points[f[0]],
            points[f[1]],
            points[f[2]]
        );

        pInter.setPoint(tri.nearestPoint(pt).rawPoint());
    }

    return pInter;  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
