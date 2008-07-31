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

        List<Pair<vector> > planeNormalPairs(featEdges.size());

        // Classify the edges around the feature point.
        forAll(featEdges, fE)
        {
            label edgeI = featEdges[fE];
            const edge& featEdge = edges[edgeI];

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

            // Identify the normals of the faces attached to feature edges to
            // identify the unique planes to be reconstructed.  The triangulated
            // surface will usually not mean that two feature edges that should
            // bound a plane are attached to the same face.

            planeNormalPairs[fE].first() = nA;
            planeNormalPairs[fE].second() = nB;
        }

        convexEdges.shrink();
        concaveEdges.shrink();

        Info<< nl << convexEdges
            << nl << concaveEdges
            << nl << planeNormalPairs
            << endl;

        List<vector> uniquePlaneNormals(featEdges.size());
        label uniquePlaneNormalI = 0;

        List<Pair<bool> > planeMatchedStatus(featEdges.size(), Pair<bool>(false,false));

        Info<< planeMatchedStatus << endl;

        // Examine the plane normals to identify unique planes.
        forAll(planeNormalPairs, nA)
        {
            const Pair<vector>& normalPairA = planeNormalPairs[nA];

            if (!planeMatchedStatus[nA].first())
            {
                const vector& nAf = normalPairA.first();

                scalar minNormalDotProduct = 1 + SMALL;

                forAll(planeNormalPairs, nB)
                {
                    if (nA == nB)
                    {
                        continue;
                    }

                    const Pair<vector>& normalPairB = planeNormalPairs[nB];

                    if (!planeMatchedStatus[nB].first())
                    {
                        const vector& nBf = normalPairB.first();

                        scalar normalDotProduct = nAf & nBf;

                        if (normalDotProduct < minNormalDotProduct)
                        {
                            minNormalDotProduct = normalDotProduct;
                        }
                    }

                    if (!planeMatchedStatus[nB].second())
                    {
                        const vector& nBs = normalPairA.second();
                    }
                }
            }

            if (!planeMatchedStatus[nA].second())
            {
                const vector& nAs = normalPairA.second();
            }

            // if (minNormalDotProduct > planeErrorAngleCos)
            // {
            //     FatalErrorIn("insertFeaturePoints")
            //         << "Could not find unique planes matching feature edges "
            //         << "at point located at "
            //         << featPt << nl
            //         << exit(FatalError);
            // }

            // const vector& matchingPNB = planeNormalsB[matchingPB];

            // uniquePlaneNormals[pA] =
            // (pNA + matchingPNB)/(mag(pNA + matchingPNB) + VSMALL);
        }

        Info<< uniquePlaneNormals << endl;
    }


    if (controls_.writeFeatureTriangulation)
    {
        writePoints("feat_allPoints.obj", false);
        writePoints("feat_points.obj", true);
        writeTriangles("feat_triangles.obj", true);
    }
}


// ************************************************************************* //
