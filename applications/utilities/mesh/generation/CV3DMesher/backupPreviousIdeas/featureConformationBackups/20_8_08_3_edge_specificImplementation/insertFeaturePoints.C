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

\*----------------------------------------------------------------------------*/

#include "CV3D.H"
#include "plane.H"
#include "triSurfaceTools.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV3D::insertFeaturePoints()
{
    Info<< nl << "Inserting feature points" << endl;

    const edgeList& edges = qSurf_.edges();
    const pointField& localPts = qSurf_.localPoints();

    labelList featPoints(0);
    labelListList featPointFeatEdges(0);

    qSurf_.extractFeatures
    (
        controls_.featAngle,
        featPoints,
        featPointFeatEdges
    );

    scalar planeErrorAngle = 0.1*controls_.featAngle;

    scalar planeErrorAngleCos = cos(mathematicalConstant::pi*planeErrorAngle/180.0);

    forAll(featPoints, i)
    {
        label ptI = featPoints[i];
        const point& featPt = localPts[ptI];

        const labelList& featEdges = featPointFeatEdges[i];

        DynamicList<label> convexEdges(0);
        DynamicList<label> concaveEdges(0);

        List<vector> planeNormals(2*featEdges.size());

        // Classify the edges around the feature point.
        forAll(featEdges, fE)
        {
            label edgeI = featEdges[fE];

            // Pick up the two faces adjacent to the feature edge
            const labelList& eFaces = qSurf_.edgeFaces()[edgeI];

            label faceA = eFaces[0];
            vector nA = qSurf_.faceNormals()[faceA];

            label faceB = eFaces[1];
            vector nB = qSurf_.faceNormals()[faceB];

            point faceAVert =
                localPts[triSurfaceTools::oppositeVertex(qSurf_, faceA, edgeI)];

            // Determine convex or concave angle
            if (((faceAVert - featPt) & nB) < 0)
            {
                // Convex feature edge
                convexEdges.append(edgeI);
            }
            else
            {
                // Concave feature edge
                concaveEdges.append(edgeI);
            }

            planeNormals[2*fE] = nA;
            planeNormals[2*fE + 1] = nB;
        }

        convexEdges.shrink();
        concaveEdges.shrink();

        // Identify the normals of the faces attached to feature edges to
        // identify the unique planes to be reconstructed.  The triangulated
        // surface will usually not mean that two feature edges that should
        // bound a plane are attached to the same face.

        List<vector> uniquePlaneNormals(featEdges.size());

        List<bool> planeMatchedStatus(2*featEdges.size(), bool(false));

        label uniquePlaneNormalI = 0;

        // Examine the plane normals to identify unique planes.
        forAll(planeNormals, nA)
        {
            const vector& normalA = planeNormals[nA];

            scalar maxNormalDotProduct = -SMALL;

            label matchingNormal = -1;

            if (!planeMatchedStatus[nA])
            {
                forAll(planeNormals, nB)
                {
                    if (nA == nB)
                    {
                        continue;
                    }

                    if (!planeMatchedStatus[nB])
                    {
                        const vector& normalB = planeNormals[nB];

                        scalar normalDotProduct = normalA & normalB;

                        if (normalDotProduct > maxNormalDotProduct)
                        {
                            maxNormalDotProduct = normalDotProduct;

                            matchingNormal = nB;
                        }
                    }
                }
            }

            if (matchingNormal >= 0)
            {
                if (maxNormalDotProduct < planeErrorAngleCos)
                {
                    FatalErrorIn("insertFeaturePoints")
                        << "Matching planes are not similar enough "
                        << "at point located at "
                        << featPt << nl
                        << exit(FatalError);
                }

                const vector& normalB = planeNormals[matchingNormal];

                uniquePlaneNormals[uniquePlaneNormalI] =
                (normalA + normalB)/(mag(normalA + normalB) + VSMALL);

                uniquePlaneNormalI++;

                planeMatchedStatus[nA] = true;
                planeMatchedStatus[matchingNormal] = true;
            }
        }

        if (concaveEdges.size() + convexEdges.size() == 3)
        {
            Info<< "3 edge feature at " << featPt << endl;

            if (concaveEdges.size() == 3)
            {
                Info << tab << "Concave feature" << endl;

                const vector& n0 = uniquePlaneNormals[0];
                const vector& n1 = uniquePlaneNormals[1];
                const vector& n2 = uniquePlaneNormals[2];

                // Interest planes outside the corner to find the common point
                plane plane0(featPt +  tols_.ppDist*n0, n0);
                plane plane1(featPt +  tols_.ppDist*n1, n1);
                plane plane2(featPt +  tols_.ppDist*n2, n2);

                point externalPt = plane2.planePlaneIntersect(plane0,plane1);
                label externalPtIndex = number_of_vertices() + 3;

                // Redfine planes to be on the corner surfaces to project through

                plane0 = plane(featPt, n0);
                plane1 = plane(featPt, n1);
                plane2 = plane(featPt, n2);

                point internalPt0 = externalPt - 2*plane0.distance(externalPt)*n0;
                label internalPt0Index = insertPoint(internalPt0, externalPtIndex);

                point internalPt1 = externalPt - 2*plane1.distance(externalPt)*n1;
                insertPoint(internalPt1, externalPtIndex);

                point internalPt2 = externalPt - 2*plane2.distance(externalPt)*n2;
                insertPoint(internalPt2, externalPtIndex);

                insertPoint(externalPt,internalPt0Index);
            }
            else if (convexEdges.size() == 3)
            {
                Info << tab << "Convex feature" << endl;

                const vector& n0 = uniquePlaneNormals[0];
                const vector& n1 = uniquePlaneNormals[1];
                const vector& n2 = uniquePlaneNormals[2];

                // Intersect planes inside the corner to find the common point
                plane plane0(featPt -  tols_.ppDist*n0, n0);
                plane plane1(featPt -  tols_.ppDist*n1, n1);
                plane plane2(featPt -  tols_.ppDist*n2, n2);

                point internalPt = plane1.planePlaneIntersect(plane2,plane0);
                label internalPtIndex = insertPoint(internalPt, number_of_vertices() + 1);

                // Redefine planes to be on the corner surfaces to project through

                plane0 = plane(featPt, n0);
                plane1 = plane(featPt, n1);
                plane2 = plane(featPt, n2);

                point externalPt0 = internalPt + 2*plane0.distance(internalPt)*n0;
                insertPoint(externalPt0, internalPtIndex);

                point externalPt1 = internalPt + 2*plane1.distance(internalPt)*n1;
                insertPoint(externalPt1, internalPtIndex);

                point externalPt2 = internalPt + 2*plane2.distance(internalPt)*n2;
                insertPoint(externalPt2, internalPtIndex);
            }
            else
            {
                Info<< tab << "Mixed feature: convex, concave = "
                    << convexEdges.size() << ", " << concaveEdges.size() << endl;

                if (concaveEdges.size() > 1)
                {
                    FatalErrorIn("insertFeaturePoints")
                        << "Assumption that only one concave edge possible is wrong."
                        << exit(FatalError);
                }

                // Find which planes are joined to the concave edge

                List<label> concaveEdgePlanes(2,label(-1));

                label concaveEdgeI = concaveEdges[0];

                // Pick up the two faces adjacent to the concave feature edge
                const labelList& eFaces = qSurf_.edgeFaces()[concaveEdgeI];

                label faceA = eFaces[0];
                vector nA = qSurf_.faceNormals()[faceA];

                scalar maxNormalDotProduct = -SMALL;

                forAll(uniquePlaneNormals, uPN)
                {
                    scalar normalDotProduct = nA & uniquePlaneNormals[uPN];

                    if (normalDotProduct > maxNormalDotProduct)
                    {
                        maxNormalDotProduct = normalDotProduct;

                        concaveEdgePlanes[0] = uPN;
                    }
                }

                label faceB = eFaces[1];
                vector nB = qSurf_.faceNormals()[faceB];

                maxNormalDotProduct = -SMALL;

                forAll(uniquePlaneNormals, uPN)
                {
                    scalar normalDotProduct = nB & uniquePlaneNormals[uPN];

                    if (normalDotProduct > maxNormalDotProduct)
                    {
                        maxNormalDotProduct = normalDotProduct;

                        concaveEdgePlanes[1] = uPN;
                    }
                }

                const vector& concaveEdgePlaneANormal =
                    uniquePlaneNormals[concaveEdgePlanes[0]];
                const vector& concaveEdgePlaneBNormal =
                    uniquePlaneNormals[concaveEdgePlanes[1]];

                label convexEdgesPlaneI;

                if (findIndex(concaveEdgePlanes, 0) == -1)
                {
                    convexEdgesPlaneI = 0;
                }
                else if (findIndex(concaveEdgePlanes, 1) == -1)
                {
                    convexEdgesPlaneI = 1;
                }
                else
                {
                    convexEdgesPlaneI = 2;
                }

                const vector& convexEdgesPlaneNormal =
                    uniquePlaneNormals[convexEdgesPlaneI];

                const edge& concaveEdge = edges[concaveEdgeI];

                // Check direction of edge, if the feature point is at the end()
                // the reverse direction.

                scalar edgeDirection = 1.0;

                if (ptI == concaveEdge.end())
                {
                    edgeDirection *= -1.0;
                }

                // Intersect planes parallel to the concave edge planes offset by ppDist
                // and the plane defined by featPt and the edge vector.
                plane planeA
                (
                    featPt + tols_.ppDist*concaveEdgePlaneANormal,
                    concaveEdgePlaneANormal
                );

                plane planeB
                (
                    featPt + tols_.ppDist*concaveEdgePlaneBNormal,
                    concaveEdgePlaneBNormal
                );

                point concaveEdgeLocalFeatPt = featPt
                  + tols_.ppDist*edgeDirection
                  * concaveEdge.vec(localPts)/concaveEdge.mag(localPts);

                plane planeF(concaveEdgeLocalFeatPt, concaveEdge.vec(localPts));

                point concaveEdgeExternalPt = planeF.planePlaneIntersect(planeA,planeB);
                label concaveEdgeExternalPtI = number_of_vertices() + 4;

                // Redfine planes to be on the corner surfaces to project through

                planeA = plane(featPt, concaveEdgePlaneANormal);
                planeB = plane(featPt, concaveEdgePlaneBNormal);

                point internalPtA = concaveEdgeExternalPt -
                    2*planeA.distance(concaveEdgeExternalPt) * concaveEdgePlaneANormal;
                label internalPtAI = insertPoint(internalPtA, concaveEdgeExternalPtI);

                point internalPtB = concaveEdgeExternalPt -
                    2*planeB.distance(concaveEdgeExternalPt) * concaveEdgePlaneBNormal;
                label internalPtBI = insertPoint(internalPtB, concaveEdgeExternalPtI);

                plane planeC(featPt, convexEdgesPlaneNormal);

                point externalPtD = internalPtA +
                    2*planeC.distance(internalPtA) * convexEdgesPlaneNormal;
                insertPoint(externalPtD, internalPtAI);

                point externalPtE = internalPtB +
                    2*planeC.distance(internalPtB) * convexEdgesPlaneNormal;
                insertPoint(externalPtE, internalPtBI);

                insertPoint(concaveEdgeExternalPt, internalPtAI);

                scalar totalAngle = 180/mathematicalConstant::pi *
                (
                    mathematicalConstant::pi +
                        acos(mag(concaveEdgePlaneANormal & concaveEdgePlaneBNormal))
                );

                if (totalAngle > controls_.maxQuadAngle)
                {
                    // Add additional mitering points

                    // Redefine planes to be inside the surface, further away from
                    // the surface than the main points.

                    // scalar guardFactor = max(3 - totalAngle/150.0, 1);
                    // // Ad-hoc function, do some analysis to determine properly

                    // planeA = plane
                    // (
                    //     featPt - guardFactor*tols_.ppDist*concaveEdgePlaneANormal,
                    //     concaveEdgePlaneANormal
                    // );

                    // planeB = plane
                    // (
                    //     featPt - guardFactor*tols_.ppDist*concaveEdgePlaneBNormal,
                    //     concaveEdgePlaneBNormal
                    // );

                    // point internalPtF = planeF.planePlaneIntersect(planeA,planeB);
                    // label internalPtFI = insertPoint(internalPtF, number_of_vertices() + 1);

                    // point externalPtG = internalPtF +
                    //     2*planeC.distance(internalPtF) * convexEdgesPlaneNormal;
                    // insertPoint(externalPtG, internalPtFI);

                    vector concaveEdgeNormal =
                        edgeDirection*concaveEdge.vec(localPts)/concaveEdge.mag(localPts);

                    scalar phi = acos(concaveEdgeNormal & -convexEdgesPlaneNormal);

                    scalar guard =
                    (
                        1 + sin(phi)*tols_.ppDist/mag
                        (
                            concaveEdgeLocalFeatPt - concaveEdgeExternalPt
                        )
                    )/cos(phi) - 1;

                    Info<< guard << endl;

                    point internalPtF = concaveEdgeExternalPt + (2 + guard) *
                    (
                        concaveEdgeLocalFeatPt - concaveEdgeExternalPt
                    );
                    label internalPtFI = insertPoint(internalPtF, number_of_vertices() + 1);

                    point externalPtG = internalPtF +
                        2*planeC.distance(internalPtF) * convexEdgesPlaneNormal;
                    insertPoint(externalPtG, internalPtFI);
                }
            }
        }
        else
        {
            Info<< concaveEdges.size() + convexEdges.size() << " edge feature."
                << " NOT IMPLEMENTED." << endl;
        }
    }

    if (controls_.writeFeatureTriangulation)
    {
        writePoints("feat_allPoints.obj", false);
        writePoints("feat_points.obj", true);
        writeTriangles("feat_triangles.obj", true);
    }
}


// ************************************************************************* //
