/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "ListOps.H"
#include "point.H"
#include "Field.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PointList>
Foam::label Foam::mergePoints
(
    const PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    typename PointList::const_reference origin
)
{
    typedef typename PointList::value_type point_type;

    const label nPoints = points.size();

    // Create an old to new point mapping array
    pointMap.resize_nocopy(nPoints);
    pointMap = -1;

    if (!nPoints)
    {
        return 0;
    }

    point_type compareOrigin = origin;
    if (origin == point_type::max)
    {
        // Use average of input points to define a comparison origin.
        // Same as sum(points)/nPoints, but handles different list types
        compareOrigin = points[0];
        for (label pointi=1; pointi < nPoints; ++pointi)
        {
            compareOrigin += points[pointi];
        }
        compareOrigin /= nPoints;
    }

    // We're comparing distance squared to origin first.
    // Say if starting from two close points:
    //     x, y, z
    //     x+mergeTol, y+mergeTol, z+mergeTol
    // Then the magSqr of both will be
    //     x^2+y^2+z^2
    //     x^2+y^2+z^2 + 2*mergeTol*(x+z+y) + mergeTol^2*...
    // so the difference will be 2*mergeTol*(x+y+z)

    const scalar mergeTolSqr = Foam::sqr(mergeTol);

    // Sort points by magSqr
    List<scalar> magSqrDist(nPoints);
    forAll(points, pointi)
    {
        magSqrDist[pointi] = magSqr(points[pointi] - compareOrigin);
    }
    labelList order(Foam::sortedOrder(magSqrDist));


    Field<scalar> sortedTol(nPoints);
    forAll(order, sortI)
    {
        const point_type& pt = points[order[sortI]];

        // Use scalar precision
        sortedTol[sortI] =
            2*mergeTol*
            (
                mag(scalar(pt.x() - compareOrigin.x()))
              + mag(scalar(pt.y() - compareOrigin.y()))
              + mag(scalar(pt.z() - compareOrigin.z()))
            );
    }

    label newPointi = 0;

    // Handle 0th point separately (is always unique)
    label pointi = order[0];
    pointMap[pointi] = newPointi++;

    for (label sortI = 1; sortI < order.size(); ++sortI)
    {
        // Get original point index
        const label pointi = order[sortI];
        const scalar mag2 = magSqrDist[order[sortI]];

        // Convert to scalar precision
        // NOTE: not yet using point_type template parameter
        const point pt
        (
            scalar(points[pointi].x()),
            scalar(points[pointi].y()),
            scalar(points[pointi].z())
        );


        // Compare to previous points to find equal one.
        label equalPointi = -1;

        for
        (
            label prevSortI = sortI - 1;
            prevSortI >= 0
         && (mag(magSqrDist[order[prevSortI]] - mag2) <= sortedTol[sortI]);
            --prevSortI
        )
        {
            const label prevPointi = order[prevSortI];

            // Convert to scalar precision
            // NOTE: not yet using point_type template parameter
            const point prevPt
            (
                scalar(points[prevPointi].x()),
                scalar(points[prevPointi].y()),
                scalar(points[prevPointi].z())
            );

            if (magSqr(pt - prevPt) <= mergeTolSqr)
            {
                // Found match.
                equalPointi = prevPointi;

                break;
            }
        }


        if (equalPointi != -1)
        {
            // Same coordinate as equalPointi. Map to same new point.
            pointMap[pointi] = pointMap[equalPointi];

            if (verbose)
            {
                Pout<< "Foam::mergePoints : Merging points "
                    << pointi << " and " << equalPointi
                    << " with coordinates:" << points[pointi]
                    << " and " << points[equalPointi]
                    << endl;
            }
        }
        else
        {
            // Differs. Store new point.
            pointMap[pointi] = newPointi++;
        }
    }

    if (verbose)
    {
        Pout<< "Foam::mergePoints : "
            << newPointi << " of " << points.size() << " unique points"
            << endl;
    }

    return newPointi;
}


template<class PointList>
bool Foam::mergePoints
(
    const PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    List<typename PointList::value_type>& newPoints,
    typename PointList::const_reference origin
)
{
    const label nUnique = mergePoints
    (
        points,
        mergeTol,
        verbose,
        pointMap,
        origin
    );

    newPoints.setSize(nUnique);
    forAll(pointMap, pointi)
    {
        newPoints[pointMap[pointi]] = points[pointi];
    }

    return (nUnique != points.size());
}


// ************************************************************************* //
