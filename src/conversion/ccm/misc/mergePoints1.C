/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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


// template<class Type, template<class> class ListType=UList>
template<class Type>
Foam::label Foam::mergePoints
(
    const bool dummy,
    const UIndirectList<Type>& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    const Type& origin
)
{
    // Create a old to new point mapping array
    pointMap.setSize(points.size());
    pointMap = -1;

    if (points.empty())
    {
        return 0;
    }

    // No normal field operations on UIndirectList.
    // Use a tmp pointField instead.
    tmp<Field<Type>> tPoints(new pointField(points));

    Type compareOrigin = origin;
    if (origin == Type::max)
    {
        compareOrigin = sum(tPoints())/points.size();
    }

    // We're comparing distance squared to origin first.
    // Say if starting from two close points:
    //     x, y, z
    //     x+mergeTol, y+mergeTol, z+mergeTol
    // Then the magSqr of both will be
    //     x^2+y^2+z^2
    //     x^2+y^2+z^2 + 2*mergeTol*(x+z+y) + mergeTol^2*...
    // so the difference will be 2*mergeTol*(x+y+z)

    const scalar mergeTolSqr = Foam::sqr(scalar(mergeTol));

    // Sort points by magSqr
    const Field<Type> d(tPoints - compareOrigin);

    List<scalar> magSqrD(d.size());
    forAll(d, pointI)
    {
        magSqrD[pointI] = magSqr(d[pointI]);
    }
    labelList order;
    sortedOrder(magSqrD, order);


    Field<scalar> sortedTol(points.size());
    forAll(order, sortI)
    {
        label pointI = order[sortI];

        // Convert to scalar precision
        const point pt
        (
            scalar(d[pointI].x()),
            scalar(d[pointI].y()),
            scalar(d[pointI].z())
        );
        sortedTol[sortI] = 2*mergeTol*(mag(pt.x())+mag(pt.y())+mag(pt.z()));
    }

    label newPointI = 0;

    // Handle 0th point separately (is always unique)
    label pointI = order[0];
    pointMap[pointI] = newPointI++;


    for (label sortI = 1; sortI < order.size(); sortI++)
    {
        // Get original point index
        label pointI = order[sortI];
        const scalar mag2 = magSqrD[order[sortI]];
        // Convert to scalar precision
        const point pt
        (
            scalar(points[pointI].x()),
            scalar(points[pointI].y()),
            scalar(points[pointI].z())
        );


        // Compare to previous points to find equal one.
        label equalPointI = -1;

        for
        (
            label prevSortI = sortI - 1;
            prevSortI >= 0
         && (mag(magSqrD[order[prevSortI]] - mag2) <= sortedTol[sortI]);
            prevSortI--
        )
        {
            label prevPointI = order[prevSortI];
            const point prevPt
            (
                scalar(points[prevPointI].x()),
                scalar(points[prevPointI].y()),
                scalar(points[prevPointI].z())
            );

            if (magSqr(pt - prevPt) <= mergeTolSqr)
            {
                // Found match.
                equalPointI = prevPointI;

                break;
            }
        }


        if (equalPointI != -1)
        {
            // Same coordinate as equalPointI. Map to same new point.
            pointMap[pointI] = pointMap[equalPointI];

            if (verbose)
            {
                Pout<< "Foam::mergePoints : Merging points "
                    << pointI << " and " << equalPointI
                    << " with coordinates:" << points[pointI]
                    << " and " << points[equalPointI]
                    << endl;
            }
        }
        else
        {
            // Differs. Store new point.
            pointMap[pointI] = newPointI++;
        }
    }

    return newPointI;
}


// ************************************************************************* //
