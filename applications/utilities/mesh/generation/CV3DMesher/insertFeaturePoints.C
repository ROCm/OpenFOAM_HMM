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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV3D::insertFeaturePoints()
{
    Info<< nl << "Inserting feature points" << endl;

    const edgeList& edges = qSurf_.edges();
    const pointField& localPts = qSurf_.localPoints();

    const labelList& featPoints = qSurf_.features().featurePoints();
    labelListList featPointFeatEdges = qSurf_.featurePointFeatureEdges();

    scalar planeErrorAngle = 0.1*(180.0 - controls_.includedAngle);

    scalar planeErrorAngleCos =
    cos(mathematicalConstant::pi*planeErrorAngle/180.0);

    forAll(featPoints, i)
    {
        label ptI = featPoints[i];
        const point& featPt = localPts[ptI];

        Info<< nl <<"Feature at " << featPt << endl;

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

        if (concaveEdges.size() == 0)
        {
            Info<< tab << "Convex feature, "
                << convexEdges.size() << " edges." << endl;

            vector cornerNormal = sum(uniquePlaneNormals);
            cornerNormal /= mag(cornerNormal);

            point internalPt =  featPt - tols_.ppDist*cornerNormal;
            label internalPtIndex =
            insertPoint(internalPt, number_of_vertices() + 1);

            forAll (uniquePlaneNormals, uPN)
            {
                const vector& n = uniquePlaneNormals[uPN];

                plane planeN = plane(featPt, n);

                point externalPt =
                internalPt + 2.0 * planeN.distance(internalPt) * n;

                insertPoint(externalPt, internalPtIndex);
            }

        }
        else if (convexEdges.size() == 0)
        {
            Info<< tab << "Concave feature, "
                << concaveEdges.size() << " edges." << endl;

            vector cornerNormal = sum(uniquePlaneNormals);
            cornerNormal /= mag(cornerNormal);

            point externalPt = featPt +  tols_.ppDist*cornerNormal;

            label externalPtIndex = number_of_vertices() + concaveEdges.size();

            label internalPtIndex = -1;

            forAll (uniquePlaneNormals, uPN)
            {
                const vector& n = uniquePlaneNormals[uPN];

                plane planeN = plane(featPt, n);

                point internalPt = externalPt
                - 2.0 * planeN.distance(externalPt) * n;

                internalPtIndex = insertPoint(internalPt, externalPtIndex);
            }

            insertPoint(externalPt,internalPtIndex);
        }
        else
        {
            Info<< tab << "Mixed feature: convex, concave = "
                << convexEdges.size() << ", " << concaveEdges.size() << endl;

            if (convexEdges.size() + concaveEdges.size() > 3)
            {
                Info<< concaveEdges.size() + convexEdges.size()
                    << " mixed edge feature."
                    << " NOT IMPLEMENTED." << endl;
            }
            else if (convexEdges.size() > concaveEdges.size())
            {
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

                // Intersect planes parallel to the concave edge planes offset
                // by ppDist and the plane defined by featPt and the edge
                // vector.
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

                // Finding the nearest point on the intersecting line to the
                // edge point.  Floating point errors often encountered using
                // planePlaneIntersect

                plane planeF(concaveEdgeLocalFeatPt, concaveEdge.vec(localPts));

                point concaveEdgeExternalPt =
                planeF.planePlaneIntersect(planeA,planeB);

                label concaveEdgeExternalPtI = number_of_vertices() + 4;

                // Redefine planes to be on the feature surfaces to project
                // through

                planeA = plane(featPt, concaveEdgePlaneANormal);

                planeB = plane(featPt, concaveEdgePlaneBNormal);

                point internalPtA = concaveEdgeExternalPt
                - 2*planeA.distance(concaveEdgeExternalPt)
                *concaveEdgePlaneANormal;

                label internalPtAI =
                insertPoint(internalPtA, concaveEdgeExternalPtI);

                point internalPtB = concaveEdgeExternalPt
                - 2*planeB.distance(concaveEdgeExternalPt)
                *concaveEdgePlaneBNormal;

                label internalPtBI =
                insertPoint(internalPtB, concaveEdgeExternalPtI);

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

                    vector concaveEdgeNormal =
                    edgeDirection*concaveEdge.vec(localPts)
                    /concaveEdge.mag(localPts);

                    scalar angleSign = 1.0;

                    if
                    (
                        qSurf_.outside
                        (
                            featPt - convexEdgesPlaneNormal*tols_.ppDist
                        )
                    )
                    {
                        angleSign = -1.0;
                    }

                    scalar phi = angleSign
                    *acos(concaveEdgeNormal & -convexEdgesPlaneNormal);

                    scalar guard =
                    (
                        1 + sin(phi)*tols_.ppDist/mag
                        (
                            concaveEdgeLocalFeatPt - concaveEdgeExternalPt
                        )
                    )/cos(phi) - 1;

                    point internalPtF = concaveEdgeExternalPt + (2 + guard)
                    *(concaveEdgeLocalFeatPt - concaveEdgeExternalPt);

                    label internalPtFI =
                    insertPoint(internalPtF, number_of_vertices() + 1);

                    point externalPtG = internalPtF +
                    2*planeC.distance(internalPtF) * convexEdgesPlaneNormal;
                    insertPoint(externalPtG, internalPtFI);
                }
            }
            else
            {
                // Find which planes are joined to the convex edge

                List<label> convexEdgePlanes(2,label(-1));

                label convexEdgeI = convexEdges[0];

                // Pick up the two faces adjacent to the convex feature edge
                const labelList& eFaces = qSurf_.edgeFaces()[convexEdgeI];

                label faceA = eFaces[0];
                vector nA = qSurf_.faceNormals()[faceA];

                scalar maxNormalDotProduct = -SMALL;

                forAll(uniquePlaneNormals, uPN)
                {
                    scalar normalDotProduct = nA & uniquePlaneNormals[uPN];

                    if (normalDotProduct > maxNormalDotProduct)
                    {
                        maxNormalDotProduct = normalDotProduct;

                        convexEdgePlanes[0] = uPN;
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

                        convexEdgePlanes[1] = uPN;
                    }
                }

                const vector& convexEdgePlaneANormal =
                    uniquePlaneNormals[convexEdgePlanes[0]];
                const vector& convexEdgePlaneBNormal =
                    uniquePlaneNormals[convexEdgePlanes[1]];

                label concaveEdgesPlaneI;

                if (findIndex(convexEdgePlanes, 0) == -1)
                {
                    concaveEdgesPlaneI = 0;
                }
                else if (findIndex(convexEdgePlanes, 1) == -1)
                {
                    concaveEdgesPlaneI = 1;
                }
                else
                {
                    concaveEdgesPlaneI = 2;
                }

                const vector& concaveEdgesPlaneNormal =
                    uniquePlaneNormals[concaveEdgesPlaneI];

                const edge& convexEdge = edges[convexEdgeI];

                // Check direction of edge, if the feature point is at the end()
                // the reverse direction.

                scalar edgeDirection = 1.0;

                if (ptI == convexEdge.end())
                {
                    edgeDirection *= -1.0;
                }

                // Intersect planes parallel to the convex edge planes offset by
                // ppDist and the plane defined by featPt and the edge vector.
                plane planeA
                (
                    featPt - tols_.ppDist*convexEdgePlaneANormal,
                    convexEdgePlaneANormal
                );

                plane planeB
                (
                    featPt - tols_.ppDist*convexEdgePlaneBNormal,
                    convexEdgePlaneBNormal
                );

                point convexEdgeLocalFeatPt = featPt
                  + tols_.ppDist*edgeDirection
                  * convexEdge.vec(localPts)/convexEdge.mag(localPts);


                // Finding the nearest point on the intersecting line to the
                // edge point.  Floating point errors often encountered using
                // planePlaneIntersect

                plane planeF(convexEdgeLocalFeatPt, convexEdge.vec(localPts));

                point convexEdgeInternalPt =
                planeF.planePlaneIntersect(planeA,planeB);

                planeA = plane(featPt, convexEdgePlaneANormal);

                planeB = plane(featPt, convexEdgePlaneBNormal);

                point externalPtA = convexEdgeInternalPt
                + 2*planeA.distance(convexEdgeInternalPt)
                * convexEdgePlaneANormal;
                label externalPtAI = number_of_vertices() + 3;

                point externalPtB = convexEdgeInternalPt
                + 2*planeB.distance(convexEdgeInternalPt)
                * convexEdgePlaneBNormal;

                label externalPtBI = number_of_vertices() + 4;

                label convexEdgeInternalPtI = insertPoint
                (
                    convexEdgeInternalPt,
                    externalPtAI
                );

                plane planeC(featPt, concaveEdgesPlaneNormal);

                point internalPtD = externalPtA -
                    2*planeC.distance(externalPtA) * concaveEdgesPlaneNormal;
                insertPoint(internalPtD, externalPtAI);

                point internalPtE = externalPtB -
                    2*planeC.distance(externalPtB) * concaveEdgesPlaneNormal;
                insertPoint(internalPtE, externalPtBI);

                insertPoint(externalPtA, convexEdgeInternalPtI);

                insertPoint(externalPtB, convexEdgeInternalPtI);
            }
        }
    }

    featureConstrainingVertices_.setSize(number_of_vertices());

    label featPtI = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        // featureConstrainingVertices_[featPtI] = vit;

        featureConstrainingVertices_[featPtI] = Vb(vit->point());

        featureConstrainingVertices_[featPtI].index() = vit->index();

        featureConstrainingVertices_[featPtI].type() = vit->type();

        featPtI++;
    }

    if (controls_.writeFeatureTriangulation)
    {
        writePoints("feat_allPoints.obj", false);
        writePoints("feat_points.obj", true);
        writeTriangles("feat_triangles.obj", true);
    }
}


void Foam::CV3D::reinsertFeaturePoints()
{
    if (featureConstrainingVertices_.size())
    {

        forAll(featureConstrainingVertices_, f)
        {
            const Point& fPt(featureConstrainingVertices_[f].point());

            uint nVert = number_of_vertices();

            Vertex_handle vh = insert(fPt);

            if (nVert == number_of_vertices())
            {
                FatalErrorIn("Foam::CV3D::reinsertFeaturePoints")
                << "Failed to reinsert feature point " << topoint(fPt)
                    << endl;
            }

            vh->index() = featureConstrainingVertices_[f].index();
            vh->type() = featureConstrainingVertices_[f].type();
        }
    }
    else
    {
        WarningIn("Foam::CV3D::reinsertFeaturePoints")
            << "No stored feature points to reinsert." << endl;
    }
}

// ************************************************************************* //
