/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "faceAreaIntersect.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::faceAreaIntersect::triSliceWithPlane
(
    const triPoints& tri,
    const plane& p,
    FixedList<triPoints, 10>& tris,
    label& nTris
)
{
    // distance to cutting plane
    FixedList<scalar, 3> d;

    // determine how many of the points are above the cutting plane
    label nCoPlanar = 0;
    label nPos = 0;
    label posI = -1;
    label negI = -1;
    bool coPlanar;
    forAll(tri, i)
    {
        d[i] = ((tri[i] - p.refPoint()) & p.normal());

        if (mag(d[i]) < 1e-10)
        {
            coPlanar = true;
            nCoPlanar++;
        }
        else
        {
            coPlanar = false;
        }

        if (!coPlanar)
        {
            if (d[i] > 0)
            {
                nPos++;
                posI = i;
            }
            else
            {
                negI = i;
            }
        }
    }

    if ((nPos == 3) || ((nPos == 2) && (nCoPlanar == 1)))
    {
        // all points above cutting plane - add triangle to list
        tris[nTris++] = tri;
    }
    else if ((nPos == 2) || ((nPos == 1) && (nCoPlanar == 1)))
    {
        // 2 points above plane, 1 below
        // resulting quad above plane split into 2 triangles

        // point under the plane
        label i0 = negI;

        // indices of remaining points
        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);

        // determine the two intersection points
        point p01 = planeIntersection(d, tri, i0, i1);
        point p02 = planeIntersection(d, tri, i0, i2);

        // forget triangle below plane
        // - decompose quad above plane into 2 triangles and add to list
        setTriPoints(tri[i1], tri[i2], p02, nTris, tris);
        setTriPoints(tri[i1], p02, p01, nTris, tris);
    }
    else if ((nPos == 1) && (nCoPlanar != 1))
    {
        // 1 point above plane, 2 below
        // resulting quad below plane split into 2 triangles

        // point above the plane
        label i0 = posI;

        // indices of remaining points
        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);

        // determine the two intersection points
        point p01 = planeIntersection(d, tri, i0, i1);
        point p02 = planeIntersection(d, tri, i0, i2);

        // forget quad below plane
        // - add triangle above plane to list
        setTriPoints(tri[i0], p01, p02, nTris, tris);
    }
    else
    {
        // all points below cutting plane - forget
    }
}


Foam::scalar Foam::faceAreaIntersect::triangleIntersect
(
    const triPoints& src,
    const triPoints& tgt,
    const vector& n
)
{
    // Work storage
    FixedList<triPoints, 10> cutTris;
    label nCutTris = 0;

    FixedList<triPoints, 10> tris;
    label nTris = 0;

    // cut source triangle with all inwards pointing faces of target triangle
    // - triangles in cutTris are inside target triangle

    // edge 0
    {
        // cut triangle src with plane and put resulting sub-triangles in
        // cutTris list

        plane pl0(tgt[0], tgt[1], tgt[1] + n);
        triSliceWithPlane(src, pl0, cutTris, nCutTris);
    }

    if (nCutTris == 0)
    {
        return 0.0;
    }

    // edge1
    {
        // cut cutTris with plane and put resulting sub-triangles in
        // tris list (re-use tris storage)

        plane pl1(tgt[1], tgt[2], tgt[2] + n);

        nTris = 0;

        for (label i = 0; i < nCutTris; i++)
        {
            triSliceWithPlane(cutTris[i], pl1, tris, nTris);
        }

        if (nTris == 0)
        {
            return 0.0;
        }
    }

    // edge2
    {
        // cut tris with plane and put resulting sub-triangles in
        // cutTris list (re-use cutTris storage)

        plane pl2(tgt[2], tgt[0], tgt[0] + n);

        nCutTris = 0;

        for (label i = 0; i < nTris; i++)
        {
            triSliceWithPlane(tris[i], pl2, cutTris, nCutTris);
        }

        if (nCutTris == 0)
        {
            return 0.0;
        }
        else
        {
            // calculate area of sub-triangles
            scalar area = 0.0;
            for (label i = 0; i < nCutTris; i++)
            {
                const triPoints& t = cutTris[i];
                area += mag(0.5*((t[1] - t[0])^(t[2] - t[0])));
            }

            return area;
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::faceAreaIntersect::faceAreaIntersect
(
    const pointField& pointsA,
    const pointField& pointsB
)
:
    pointsA_(pointsA),
    pointsB_(pointsB)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::faceAreaIntersect::calc
(
    const face& faceA,
    const face& faceB,
    const vector& n
)
{
    // split faces into triangles
    DynamicList<face> trisA;
    DynamicList<face> trisB;

//    if (useTriangleFan)
//    {
//        triangleFan(faceA, trisA);
//        triangleFan(faceB, trisB);
//    }
//    else
//    {
//        faceA.triangles(pointsA_, trisA);
//        faceB.triangles(pointsB_, trisB);
//    }

    faceA.triangles(pointsA_, trisA);
    faceB.triangles(pointsB_, trisB);

    // intersect triangles
    scalar totalArea = 0.0;
    forAll(trisA, tA)
    {
        triPoints tpA = getTriPoints(pointsA_, trisA[tA], false);

        if (triArea(tpA) > ROOTVSMALL)
        {
            forAll(trisB, tB)
            {
                triPoints tpB = getTriPoints(pointsB_, trisB[tB], true);

                if (triArea(tpB) > ROOTVSMALL)
                {
                    totalArea += triangleIntersect(tpA, tpB, n);
                }
            }
        }
    }

    return totalArea;
}


// ************************************************************************* //
