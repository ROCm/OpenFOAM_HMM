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

\*----------------------------------------------------------------------------*/

#include "triSurfaceMeshes.H"
#include "Random.H"
#include "Time.H"
#include "SortableList.H"
#include "IOmanip.H"
#include "plane.H"
#include "SortableList.H"
#include "triSurfaceTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(triSurfaceMeshes, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate sum of distance to surfaces.
Foam::scalar Foam::triSurfaceMeshes::sumDistSqr
(
    const labelList& surfacesToTest,
    const scalar initDistSqr,
    const point& pt
) const
{
    scalar sum = 0;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        pointIndexHit hit(operator[](surfI).findNearest(pt, initDistSqr));

        // Note: make it fall over if not hit.
        sum += magSqr(hit.hitPoint()-pt);
    }
    return sum;
}


// Reflects the point furthest away around the triangle centre by a factor fac.
// (triangle centre is the average of all points but the ihi. pSum is running
//  sum of all points)
Foam::scalar Foam::triSurfaceMeshes::tryMorphTet
(
    const labelList& surfacesToTest,
    const scalar initDistSqr,
    List<vector>& p,
    List<scalar>& y,
    vector& pSum,
    const label ihi,
    const scalar fac
) const
{
    scalar fac1 = (1.0-fac)/vector::nComponents;
    scalar fac2 = fac1-fac;

    vector ptry = pSum*fac1-p[ihi]*fac2;

    scalar ytry = sumDistSqr(surfacesToTest, initDistSqr, ptry);

    if (ytry < y[ihi])
    {
        y[ihi] = ytry;
        pSum += ptry - p[ihi];
        p[ihi] = ptry;
    }
    return ytry;
}


bool Foam::triSurfaceMeshes::morphTet
(
    const labelList& surfacesToTest,
    const scalar initDistSqr,
    const scalar convergenceDistSqr,
    const label maxIter,
    List<vector>& p,
    List<scalar>& y
) const
{
    vector pSum = sum(p);

    autoPtr<OFstream> str;
    label vertI = 0;
    if (debug)
    {
        Pout<< "triSurfaceMeshes::morphTet : intersection of "
            << IndirectList<fileName>(names(), surfacesToTest)()
            << " starting from points:" << p << endl;
        str.reset(new OFstream("track.obj"));
        meshTools::writeOBJ(str(), p[0]);
        vertI++;
    }

    for (label iter = 0; iter < maxIter; iter++)
    {
        // Get the indices of highest, second-highest and lowest values.
        label ihi, inhi, ilo;
        {
            SortableList<scalar> sortedY(y);
            ilo = sortedY.indices()[0];
            ihi = sortedY.indices()[sortedY.size()-1];
            inhi = sortedY.indices()[sortedY.size()-2];
        }

        if (debug)
        {
            Pout<< "Iteration:" << iter
                << " lowest:" << y[ilo] << " highest:" << y[ihi]
                << " points:" << p << endl;

            meshTools::writeOBJ(str(), p[ilo]);
            vertI++;
            str()<< "l " << vertI-1 << ' ' << vertI << nl;
        }

        if (y[ihi] < convergenceDistSqr)
        {
            // Get point on 0th surface.
            Swap(p[0], p[ilo]);
            Swap(y[0], y[ilo]);
            return true;
        }

        // Reflection: point furthest away gets reflected.
        scalar ytry = tryMorphTet
        (
            surfacesToTest,
            10*y[ihi],             // search box.
            p,
            y,
            pSum,
            ihi,
            -1.0
        );

        if (ytry <= y[ilo])
        {
            // If in right direction (y lower) expand by two.
            ytry = tryMorphTet(surfacesToTest, 10*y[ihi], p, y, pSum, ihi, 2.0);
        }
        else if (ytry >= y[inhi])
        {
            // If inside tet try contraction.

            scalar ysave = y[ihi];

            ytry = tryMorphTet(surfacesToTest, 10*y[ihi], p, y, pSum, ihi, 0.5);

            if (ytry >= ysave)
            {
                // Contract around lowest point.
                forAll(p, i)
                {
                    if (i != ilo)
                    {
                        p[i] = 0.5*(p[i] + p[ilo]);
                        y[i] = sumDistSqr(surfacesToTest, y[ihi], p[i]);
                    }
                }
                pSum = sum(p);
            }
        }
    }

    if (debug)
    {
        meshTools::writeOBJ(str(), p[0]);
        vertI++;
        str()<< "l " << vertI-1 << ' ' << vertI << nl;
    }

    // Failure to converge. Return best guess so far.
    label ilo = findMin(y);
    Swap(p[0], p[ilo]);
    Swap(y[0], y[ilo]);
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::triSurfaceMeshes::triSurfaceMeshes
(
    const IOobject& io,
    const fileNameList& names
)
:
    PtrList<triSurfaceMesh>(names.size()),
    allSurfaces_(identity(names.size()))
{
    forAll(names, i)
    {
        autoPtr<IOobject> surfaceIO = io.clone();
        surfaceIO().rename(names[i]);

        Info<< "Loading surface " << surfaceIO().name() << endl;

        //fileName fullPath = surfaceIO().filePath();
        //
        //if (fullPath.size() == 0)
        //{
        //    FatalErrorIn
        //    (
        //        "triSurfaceMeshes::triSurfaceMeshes"
        //        "(const IOobject&, const fileNameList&)"
        //    )   << "Cannot load surface " << surfaceIO().name()
        //        << " starting from path " << surfaceIO().path()
        //        << exit(FatalError);
        //}

        set(i, new triSurfaceMesh(surfaceIO()));

        if (Pstream::master())
        {
            string oldPrefix(Pout.prefix());
            Pout.prefix() += "    ";
            operator[](i).writeStats(Pout);
            Pout.prefix() = oldPrefix;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileNameList Foam::triSurfaceMeshes::allNames(const IOobject& io)
{
    return readDir(io.path(), fileName::FILE);
}


Foam::fileNameList Foam::triSurfaceMeshes::names() const
{
    fileNameList surfNames(size());

    forAll(surfNames, surfI)
    {
        surfNames[surfI] = operator[](surfI).IOobject::name();
    }
    return surfNames;
}


// Find any intersection
Foam::label Foam::triSurfaceMeshes::findAnyIntersection
(
    const labelList& surfaces,
    const point& start,
    const point& end,
    pointIndexHit& hitInfo
) const
{
    forAll(surfaces, i)
    {
        label surfI = surfaces[i];

        hitInfo = operator[](surfI).findLineAny(start, end);

        if (hitInfo.hit())
        {
            return surfI;
        }
    }
    return -1;
}


Foam::label Foam::triSurfaceMeshes::findAnyIntersection
(
    const point& start,
    const point& end,
    pointIndexHit& hitInfo
) const
{
    return findAnyIntersection
    (
        allSurfaces_,
        start,
        end,
        hitInfo
    );
}


// Find intersections of edge nearest to both endpoints.
void Foam::triSurfaceMeshes::findAllIntersections
(
    const labelList& surfaces,
    const point& start,
    const point& end,

    labelList& surfacesIndex,
    List<pointIndexHit>& surfaceHitInfo
) const
{
    DynamicList<label> hitSurfaces(surfaces.size());
    DynamicList<pointIndexHit> hitInfos(surfaces.size());
    DynamicList<scalar> hitDistSqr(surfaces.size());

    const vector dirVec(end-start);
    const scalar magSqrDirVec(magSqr(dirVec));
    const vector smallVec(Foam::sqrt(SMALL)*dirVec);

    forAll(surfaces, i)
    {
        label surfI = surfaces[i];

        // Current starting point of ray.
        point pt = start;

        while (true)
        {
            // See if any intersection between pt and end
            pointIndexHit inter = operator[](surfI).findLine(pt, end);

            if (!inter.hit())
            {
                break;
            }
            hitSurfaces.append(surfI);
            hitInfos.append(inter);
            hitDistSqr.append(magSqr(inter.hitPoint() - start));

            pt = inter.hitPoint() + smallVec;

            if (((pt-start)&dirVec) > magSqrDirVec)
            {
                // Adding smallVec has taken us beyond end
                break;
            }
        }
    }
    surfacesIndex.setSize(hitSurfaces.size());
    surfaceHitInfo.setSize(hitSurfaces.size());

    if (hitSurfaces.size() > 0)
    {
        // Sort and transfer to arguments

        hitSurfaces.shrink();
        hitInfos.shrink();
        hitDistSqr.shrink();

        // Sort from start to end.
        SortableList<scalar> sortedDist(hitDistSqr);

        forAll(sortedDist.indices(), newPos)
        {
            label oldPos = sortedDist.indices()[newPos];
            surfacesIndex[newPos] = hitSurfaces[oldPos];
            surfaceHitInfo[newPos] = hitInfos[oldPos];
        }
    }
}


void Foam::triSurfaceMeshes::findAllIntersections
(
    const point& start,
    const point& end,

    labelList& surfacesIndex,
    List<pointIndexHit>& surfaceHitInfo
) const
{
    findAllIntersections
    (
        allSurfaces_,
        start,
        end,
        surfacesIndex,
        surfaceHitInfo
    );
}


// Find intersections of edge nearest to both endpoints.
void Foam::triSurfaceMeshes::findNearestIntersection
(
    const labelList& surfaces,
    const point& start,
    const point& end,

    label& surface1,
    pointIndexHit& hit1,
    label& surface2,
    pointIndexHit& hit2
) const
{
    surface1 = -1;
    // Initialize to endpoint
    hit1 = pointIndexHit(false, end, -1);

    forAll(surfaces, i)
    {
        label surfI = surfaces[i];

        if (hit1.rawPoint() == start)
        {
            break;
        }

        // See if any intersection between start and current nearest
        pointIndexHit inter = operator[](surfI).findLine
        (
            start,
            hit1.rawPoint()
        );

        if (inter.hit())
        {
            hit1 = inter;
            surface1 = surfI;
        }
    }


    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = pointIndexHit(hit1);

    if (hit1.hit())
    {
        // Test from the end side.
        forAll(surfaces, i)
        {
            label surfI = surfaces[i];

            if (hit2.rawPoint() == end)
            {
                break;
            }

            // See if any intersection between end and current nearest
            pointIndexHit inter = operator[](surfI).findLine
            (
                end,
                hit2.rawPoint()
            );

            if (inter.hit())
            {
                hit2 = inter;
                surface2 = surfI;
            }
        }
    }
}


void Foam::triSurfaceMeshes::findNearestIntersection
(
    const point& start,
    const point& end,

    label& surface1,
    pointIndexHit& hit1,
    label& surface2,
    pointIndexHit& hit2
) const
{
    findNearestIntersection
    (
        allSurfaces_,
        start,
        end,
        surface1,
        hit1,
        surface2,
        hit2
    );
}


// Find nearest. Return -1 or nearest point
Foam::label Foam::triSurfaceMeshes::findNearest
(
    const labelList& surfaces,
    const point& pt,
    const scalar nearestDistSqr,
    pointIndexHit& nearestHit
) const
{
    // nearest surface
    label minSurface = -1;
    scalar minDistSqr = Foam::sqr(GREAT);

    forAll(surfaces, i)
    {
        label surfI = surfaces[i];

        pointIndexHit hit
        (
            operator[](surfI).findNearest(pt, nearestDistSqr)
        );

        if (hit.hit())
        {
            scalar distSqr = magSqr(hit.hitPoint()-pt);

            if (distSqr < minDistSqr)
            {
                minDistSqr = distSqr;
                minSurface = surfI;
                nearestHit = hit;
            }
        }
    }

    if (minSurface == -1)
    {
        // maxLevel unchanged. No interesting surface hit.
        nearestHit.setMiss();
    }

    return minSurface;
}


// Find nearest. Return -1 or nearest point
Foam::label Foam::triSurfaceMeshes::findNearest
(
    const point& pt,
    const scalar nearestDistSqr,
    pointIndexHit& nearestHit
) const
{
    return findNearest
    (
        allSurfaces_,
        pt,
        nearestDistSqr,
        nearestHit
    );
}


Foam::label Foam::triSurfaceMeshes::findNearestAndClassify
(
    const labelList& surfacesToTest,
    const point& pt,
    const scalar nearestDistSqr,
    surfaceLocation& nearest
) const
{
    pointIndexHit nearestInfo;
    label surfI = findNearest
    (
        surfacesToTest,
        pt,
        nearestDistSqr,
        nearestInfo
    );

    if (surfI != -1)
    {
        nearest = triSurfaceTools::classify
        (
            operator[](surfI),
            nearestInfo.index(),        // triangle
            nearestInfo.rawPoint()      // point on edge/inside triangle
        );
    }

    return surfI;
}


//- Calculate point which is on a set of surfaces.
Foam::surfaceLocation Foam::triSurfaceMeshes::facesIntersectionTrack
(
    const labelList& surfacesToTest,
    const scalar initialDistSqr,
    const scalar convergenceDistSqr,
    const point& start
) const
{
    autoPtr<OFstream> str;
    label vertI = 0;

    if (debug)
    {
        str.reset(new OFstream("track.obj"));
    }

    surfaceLocation current
    (
        pointIndexHit
        (
            false,
            start,
            -1
        ),
        triPointRef::NONE,
        -1
    );

    // Dump start point
    if (str.valid())
    {
        Pout<< "Starting at " << current.info()
            << " written as point " << vertI
            << endl;
        meshTools::writeOBJ(str(), current.rawPoint());
        vertI++;
    }

    // Now slide current along all faces until it settles
    label iter = 0;

    const label maxIter = 10;

    for (; iter < maxIter; iter++)
    {
        // - store old position
        // - slide it along all faces
        // - if it hasn't changed position exit loop

        if (debug)
        {
            Pout<< endl
                << "Iteration " << iter << endl
                << "-----------" << endl;
        }

        point oldCurrent = current.rawPoint();

        forAll(surfacesToTest, i)
        {
            label surfI = surfacesToTest[i];

            // Project pt onto surf
            // ~~~~~~~~~~~~~~~~~~~~
            point copy(current.rawPoint()); // need to work on copy
            if
            (
                findNearestAndClassify
                (
                    labelList(1, surfI),
                    copy,
                    initialDistSqr,             // initialDistSqr
                    current
                )
             == -1
            )
            {
                FatalErrorIn("triSurfaceMeshes::facesIntersectionTrack(..)")
                    << "Did not find a point on surface " << surfI
                    << " which is within sqrt(" << initialDistSqr
                    << ") of " << copy
                    << abort(FatalError);
            }


            if (debug)
            {
                Pout<< "Nearest onto surface " << surfI
                    << ' ' << operator[](surfI).IOobject::name() << " : "
                    << current.info() << " written as point " << vertI
                    << endl;
            }

            // Dump current
            if (str.valid())
            {
                meshTools::writeOBJ(str(), current.rawPoint());
                vertI++;
                str()<< "l " << vertI-1 << ' ' << vertI << nl;
            }

            // Find corresponding point on next surface
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            label nextSurfI = surfacesToTest[surfacesToTest.fcIndex(i)];

            surfaceLocation next;
            findNearestAndClassify
            (
                labelList(1, nextSurfI),
                point(current.rawPoint()),
                initialDistSqr,                 // initialDistSqr
                next
            );

            // Since next is on a different surface we cannot compare
            // indices (like in snapToEnd) so make sure this does not
            // happen.
            next.elementType() = triPointRef::NONE;
            next.setIndex(-123);

            if (debug)
            {
                Pout<< "Nearest onto next surface " << nextSurfI
                    << ' ' << operator[](nextSurfI).IOobject::name() << " : "
                    << next.info() << endl;
            }

            if
            (
                magSqr(current.rawPoint() - next.rawPoint())
              > convergenceDistSqr
            )
            {
                // Track from current to next
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~
                //vector n = current.normal(surfaces[surfI]);

                current = triSurfaceTools::trackToEdge
                (
                    operator[](surfI),
                    current,
                    next,
                    plane(start, current.rawPoint(), next.rawPoint())
                );

                if (debug)
                {
                    Pout<< "Sliding along surface "
                        << surfI << ' ' << operator[](surfI).IOobject::name()
                        << " in direction of "
                        << nextSurfI << ' '
                        << operator[](nextSurfI).IOobject::name()
                        << " stopped at " << current.info()
                        << " written as point " << vertI << endl;
                }
                if (str.valid())
                {
                    meshTools::writeOBJ(str(), current.rawPoint());
                    vertI++;
                    str()<< "l " << vertI-1 << ' ' << vertI << nl;
                }
            }
        }

        scalar distSqr = magSqr(oldCurrent - current.rawPoint());

        if (debug)
        {
            Pout<< "distSqr:" << distSqr
                << " convergenceDistSqr:" << convergenceDistSqr << endl;
        }

        if (distSqr < convergenceDistSqr)
        {
            break;
        }
    }

    if (iter == maxIter)
    {
        FatalErrorIn("triSurfaceMeshes::facesIntersectionTrack(..)")
            << "Did not converge in " << iter
            << " iterations to find a point which is on surfaces "
            << IndirectList<fileName>(names(), surfacesToTest)() << nl
            << "Start point:" << start << nl
            << "Current nearest:" << current.info() << nl
            << "Please check that your surfaces actually intersect and"
            << " that the starting point is close to the intersection point."
            << abort(FatalError);
    }

    return current;
}


//- Calculate point which is on a set of surfaces.
Foam::pointIndexHit Foam::triSurfaceMeshes::facesIntersection
(
    const labelList& surfacesToTest,
    const scalar initDistSqr,
    const scalar convergenceDistSqr,
    const point& start
) const
{
    // Get four starting points. Take these as the projection of the
    // starting point onto the surfaces and the mid point
    List<point> nearest(surfacesToTest.size()+1);

    point sumNearest = vector::zero;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        pointIndexHit hit
        (
            operator[](surfI).findNearest(start, initDistSqr)
        );

        if (hit.hit())
        {
            nearest[i] = hit.hitPoint();
            sumNearest += nearest[i];
        }
        else
        {
            FatalErrorIn
            (
                "triSurfaceMeshes::facesIntersection"
                "(const labelList&, const scalar, const scalar, const point&)"
            )   << "Did not find point within distance "
                << initDistSqr << " of starting point " << start
                << " on surface " << operator[](surfI).IOobject::name()
                << abort(FatalError);
        }
    }

    nearest[nearest.size()-1] = sumNearest / surfacesToTest.size();


    // Get the sum of distances (initial evaluation)
    List<scalar> nearestDist(nearest.size());

    forAll(nearestDist, i)
    {
        nearestDist[i] = sumDistSqr(surfacesToTest, initDistSqr, nearest[i]);
    }


    //- Downhill Simplex method

    bool converged = morphTet
    (
        surfacesToTest,
        initDistSqr,
        convergenceDistSqr,
        2000,
        nearest,
        nearestDist
    );


    pointIndexHit intersection;

    if (converged)
    {
        // Project nearest onto 0th surface.
        intersection = operator[](surfacesToTest[0]).findNearest
        (
            nearest[0],
            nearestDist[0]
        );
    }

    //if (!intersection.hit())
    //{
    //    // Restart
    //    scalar smallDist = Foam::sqr(convergenceDistSqr);
    //    nearest[0] = intersection.hitPoint();
    //    nearest[1] = nearest[0];
    //    nearest[1].x() += smallDist;
    //    nearest[2] = nearest[0];
    //    nearest[2].y() += smallDist;
    //    nearest[3] = nearest[0];
    //    nearest[3].z() += smallDist;
    //
    //    forAll(nearestDist, i)
    //    {
    //        nearestDist[i] = sumDistSqr
    //        (
    //            surfacesToTest,
    //            initDistSqr,
    //            nearest[i]
    //        );
    //    }
    //
    //    intersection = morphTet
    //    (
    //        surfacesToTest,
    //        initDistSqr,
    //        convergenceDistSqr,
    //        1000,
    //        nearest,
    //        nearestDist
    //    );
    //}

    return intersection;
}



// ************************************************************************* //
