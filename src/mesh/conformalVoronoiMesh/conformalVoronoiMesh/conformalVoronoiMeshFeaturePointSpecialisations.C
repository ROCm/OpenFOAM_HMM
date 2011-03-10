/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::conformalVoronoiMesh::insertSpecialisedFeaturePoint
(
    const extendedFeatureEdgeMesh& feMesh,
    label ptI
)
{
    const labelList& pEds(feMesh.pointEdges()[ptI]);

    if (pEds.size() != 3)
    {
        // Only three edge specialisations available

        return false;
    }

    label nExternal = 0;
    label nInternal = 0;
    label nFlat = 0;
    label nOpen = 0;
    label nMultiple = 0;

    List<extendedFeatureEdgeMesh::edgeStatus> allEdStat(pEds.size());

    forAll(pEds, i)
    {
        label edgeI = pEds[i];

        extendedFeatureEdgeMesh::edgeStatus& eS = allEdStat[i];

        eS = feMesh.getEdgeStatus(edgeI);

        switch (eS)
        {
            case extendedFeatureEdgeMesh::EXTERNAL:
            {
                nExternal++;
                break;
            }
            case extendedFeatureEdgeMesh::INTERNAL:
            {
                nInternal++;
                break;
            }
            case extendedFeatureEdgeMesh::FLAT:
            {
                nFlat++;
                break;
            }
            case extendedFeatureEdgeMesh::OPEN:
            {
                nOpen++;
                break;
            }
            case extendedFeatureEdgeMesh::MULTIPLE:
            {
                nMultiple++;
                break;
            }
            case extendedFeatureEdgeMesh::NONE:
            {
                break;
            }
        }
    }

    if (nExternal == 2 && nInternal == 1)
    {
        // if (!geometryToConformTo_.positionOnThisProc(pt))
        // {
        //     continue;
        // }

        // Info<< "nExternal == 2 && nInternal == 1" << endl;

        // const Foam::point& featPt = feMesh.points()[ptI];

        // scalar ppDist = pointPairDistance(featPt);

        // const vectorField& normals = feMesh.normals();

        // const labelListList& edgeNormals = feMesh.edgeNormals();

        // label concaveEdgeI = pEds
        // [
        //     findIndex(allEdStat, extendedFeatureEdgeMesh::INTERNAL)
        // ];

        // // // Find which planes are joined to the concave edge

        // // List<label> concaveEdgePlanes(2,label(-1));

        // // label concaveEdgeI = concaveEdges[0];

        // // // Pick up the two faces adjacent to the concave feature edge
        // // const labelList& eFaces = qSurf_.edgeFaces()[concaveEdgeI];

        // // label faceA = eFaces[0];

        // // vector nA = qSurf_.faceNormals()[faceA];

        // // scalar maxNormalDotProduct = -SMALL;

        // // forAll(uniquePlaneNormals, uPN)
        // // {
        // //     scalar normalDotProduct = nA & uniquePlaneNormals[uPN];

        // //     if (normalDotProduct > maxNormalDotProduct)
        // //     {
        // //         maxNormalDotProduct = normalDotProduct;

        // //         concaveEdgePlanes[0] = uPN;
        // //     }
        // // }

        // // label faceB = eFaces[1];
        // // vector nB = qSurf_.faceNormals()[faceB];

        // // maxNormalDotProduct = -SMALL;

        // // forAll(uniquePlaneNormals, uPN)
        // // {
        // //     scalar normalDotProduct = nB & uniquePlaneNormals[uPN];

        // //     if (normalDotProduct > maxNormalDotProduct)
        // //     {
        // //         maxNormalDotProduct = normalDotProduct;

        // //         concaveEdgePlanes[1] = uPN;
        // //     }
        // // }

        // // const vector& concaveEdgePlaneANormal =
        // // uniquePlaneNormals[concaveEdgePlanes[0]];

        // // const vector& concaveEdgePlaneBNormal =
        // // uniquePlaneNormals[concaveEdgePlanes[1]];

        // // label convexEdgesPlaneI;

        // // if (findIndex(concaveEdgePlanes, 0) == -1)
        // // {
        // //     convexEdgesPlaneI = 0;
        // // }
        // // else if (findIndex(concaveEdgePlanes, 1) == -1)
        // // {
        // //     convexEdgesPlaneI = 1;
        // // }
        // // else
        // // {
        // //     convexEdgesPlaneI = 2;
        // // }

        // // const vector& convexEdgesPlaneNormal =
        // // uniquePlaneNormals[convexEdgesPlaneI];

        // // const edge& concaveEdge = edges[concaveEdgeI];

        // // // Check direction of edge, if the feature point is at the end()
        // // // the reverse direction.

        // // scalar edgeDirection = 1.0;

        // // if (ptI == concaveEdge.end())
        // // {
        // //     edgeDirection *= -1.0;
        // // }


        // const vector& concaveEdgePlaneANormal =
        //     normals[edgeNormals[concaveEdgeI][0]];

        // const vector& concaveEdgePlaneBNormal =
        //     normals[edgeNormals[concaveEdgeI][1]];

        // // Intersect planes parallel to the concave edge planes offset
        // // by ppDist and the plane defined by featPt and the edge
        // // vector.
        // plane planeA
        // (
        //     featPt + ppDist*concaveEdgePlaneANormal,
        //     concaveEdgePlaneANormal
        // );

        // plane planeB
        // (
        //     featPt + ppDist*concaveEdgePlaneBNormal,
        //     concaveEdgePlaneBNormal
        // );

        // const vector& concaveEdgeDir = feMesh.edgeDirection
        // (
        //     concaveEdgeI,
        //     ptI
        // );

        // Foam::point concaveEdgeLocalFeatPt = featPt + ppDist*concaveEdgeDir;

        // // Finding the nearest point on the intersecting line to the
        // // edge point.  Floating point errors often encountered using
        // // planePlaneIntersect

        // plane planeF(concaveEdgeLocalFeatPt, concaveEdgeDir);

        // Foam::point concaveEdgeExternalPt = planeF.planePlaneIntersect
        // (
        //     planeA,
        //     planeB
        // );

        // label concaveEdgeExternalPtI = number_of_vertices() + 4;

        // // Redefine planes to be on the feature surfaces to project
        // // through

        // planeA = plane(featPt, concaveEdgePlaneANormal);

        // planeB = plane(featPt, concaveEdgePlaneBNormal);

        // Foam::point internalPtA =
        //     concaveEdgeExternalPt
        //   - 2*planeA.distance(concaveEdgeExternalPt)
        //     *concaveEdgePlaneANormal;

        // label internalPtAI = insertPoint
        // (
        //     internalPtA,
        //     concaveEdgeExternalPtI
        // );

        // Foam::point internalPtB =
        //     concaveEdgeExternalPt
        //   - 2*planeB.distance(concaveEdgeExternalPt)
        //    *concaveEdgePlaneBNormal;

        // label internalPtBI = insertPoint
        // (
        //     internalPtB,
        //     concaveEdgeExternalPtI
        // );

        // // TEMPORARY VARIABLE TO TEST
        // vector convexEdgesPlaneNormal = -concaveEdgeDir;

        // plane planeC(featPt, convexEdgesPlaneNormal);
        // Foam::point externalPtD =
        //     internalPtA
        //   + 2*planeC.distance(internalPtA)*convexEdgesPlaneNormal;

        // insertPoint(externalPtD, internalPtAI);

        // Foam::point externalPtE =
        //     internalPtB
        //   + 2*planeC.distance(internalPtB)*convexEdgesPlaneNormal;

        // insertPoint(externalPtE, internalPtBI);

        // insertPoint(concaveEdgeExternalPt, internalPtAI);

        // Info<< nl << "# featPt " << endl;
        // meshTools::writeOBJ(Info, featPt);
        // Info<< "# internalPtA" << endl;
        // meshTools::writeOBJ(Info, internalPtA);
        // Info<< "# internalPtB" << endl;
        // meshTools::writeOBJ(Info, internalPtB);
        // Info<< "# externalPtD" << endl;
        // meshTools::writeOBJ(Info, externalPtD);
        // Info<< "# externalPtE" << endl;
        // meshTools::writeOBJ(Info, externalPtE);

        // scalar totalAngle = radToDeg
        // (
        //     constant::mathematical::pi +
        //     acos(mag(concaveEdgePlaneANormal & concaveEdgePlaneBNormal))
        // );

        // if (totalAngle > cvMeshControls().maxQuadAngle())
        // {
        //     // Add additional mitering points

        //     scalar angleSign = 1.0;

        //     if
        //     (
        //         geometryToConformTo_.outside
        //         (
        //             featPt - convexEdgesPlaneNormal*ppDist
        //         )
        //     )
        //     {
        //         angleSign = -1.0;
        //     }

        //     scalar phi =
        //         angleSign*acos(concaveEdgeDir & -convexEdgesPlaneNormal);

        //     scalar guard =
        //     (
        //         1 + sin(phi)*ppDist/mag
        //         (
        //             concaveEdgeLocalFeatPt - concaveEdgeExternalPt
        //         )
        //     )/cos(phi) - 1;

        //     Foam::point internalPtF =
        //         concaveEdgeExternalPt
        //       + (2 + guard)*(concaveEdgeLocalFeatPt - concaveEdgeExternalPt);

        //     label internalPtFI =
        //         insertPoint(internalPtF, number_of_vertices() + 1);

        //     Foam::point externalPtG =
        //         internalPtF
        //       + 2*planeC.distance(internalPtF) * convexEdgesPlaneNormal;

        //     insertPoint(externalPtG, internalPtFI);
        // }

        // return true;

        return false;
    }
    else if (nExternal == 1 && nInternal == 2)
    {
        // Info<< "nExternal == 1 && nInternal == 2" << endl;

        return false;
    }

    return false;
}


// ************************************************************************* //
