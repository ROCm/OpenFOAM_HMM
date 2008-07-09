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
    labelList featEdges(qSurf_.extractFeatures3D(controls_.featAngle));

    const pointField& localPts = qSurf_.localPoints();

//     forAll(featEdges, i)
//     {
//         label edgeI = featEdges[i];
//         const edge& featEdge = qSurf_.edges()[edgeI];
// 
//         // Get the feature point as the mid-point of the edge and convert to 2D
//         point featPt = featEdge.centre(qSurf_.localPoints());
// 
//         // Pick up the two faces adjacent to the feature edge
//         const labelList& eFaces = qSurf_.edgeFaces()[edgeI];
// 
//         label faceA = eFaces[0];
//         vector nA = qSurf_.faceNormals()[faceA];
// 
//         label faceB = eFaces[1];
//         vector nB = qSurf_.faceNormals()[faceB];
// 
//         // Intersect planes parallel to faceA and faceB offset by ppDist.
//         plane planeA(featPt - tols_.ppDist*nA, nA);
//         plane planeB(featPt - tols_.ppDist*nB, nB);
//         plane::ray interLine(planeA.planeIntersect(planeB));
// 
//         // The reference point is where ... TODO
//         point refPt =
//         (
//             interLine.refPoint()
//         );
// 
//         point faceAVert =
//         (
//             localPts[triSurfaceTools::oppositeVertex(qSurf_, faceA, edgeI)]
//         );
// 
//         // Determine convex or concave angle
//         if (((faceAVert - featPt) & nB) < 0)
//         {
//             // Convex. So refPt will be inside domain and hence a master point
// 
//             // Insert the master point refering the first slave
//             label masterPtIndex = insertPoint(refPt, number_of_vertices() + 1);
// 
//             // Insert the slave points by reflecting refPt in both faces.
//             // with each slave refering to the master
// 
//             point reflectedA = refPt + 2*((featPt - refPt) & nA)*nA;
//             insertPoint(reflectedA, masterPtIndex);
// 
//             point reflectedB = refPt + 2*((featPt - refPt) & nB)*nB;
//             insertPoint(reflectedB, masterPtIndex);
//         }
//         else
//         {
//             // Concave. master and reflected points inside the domain.
//             // Generate reflected master to be outside.
//             point reflMasterPt = refPt + 2*(featPt - refPt);
// 
//             // Reflect refPt in both faces.
//             point reflectedA = 
//                 reflMasterPt + 2*((featPt - reflMasterPt) & nA)*nA;
// 
//             point reflectedB =
//                 reflMasterPt + 2*((featPt - reflMasterPt) & nB)*nB;
// 
//             // Total angle around the concave feature
// 
//             scalar totalAngle = 
//                 180*(mathematicalConstant::pi + acos(mag(nA & nB)))
//                /mathematicalConstant::pi;
// 
// 
//             // Number of quadrants the angle should be split into
//             int nQuads = int(totalAngle/controls_.maxQuadAngle) + 1;
// 
//             // The number of additional master points needed to obtain the
//             // required number of quadrants.
//             int nAddPoints = min(max(nQuads - 2, 0), 2);
// 
//             // index of reflMaster
//             label reflectedMaster = number_of_vertices() + 2 + nAddPoints;
// 
//             // Master A is inside.
//             label reflectedAI = insertPoint(reflectedA, reflectedMaster);
// 
//             // Master B is inside.
//             insertPoint(reflectedB, reflectedMaster);
// 
//             if (nAddPoints == 1)
//             {
//                 // One additinal point is the reflection of the slave point,
//                 // i.e. the original reference point
//                 insertPoint(refPt, reflectedMaster);
//             }
//             else if (nAddPoints == 2)
//             {
//                 point reflectedAa = refPt - ((featPt - reflMasterPt) & nB)*nB;
//                 insertPoint(reflectedAa, reflectedMaster);
// 
//                 point reflectedBb = refPt - ((featPt - reflMasterPt) & nA)*nA;
//                 insertPoint(reflectedBb, reflectedMaster);
//             }
// 
//             // Slave is outside.
//             insertPoint(reflMasterPt, reflectedAI);
//         }
//     }

    if (controls_.writeFeatureTriangulation)
    {
        writePoints("feat_allPoints.obj", false);
        writePoints("feat_points.obj", true);
//         writeFaces("feat_allFaces.obj", false);
//         writeFaces("feat_faces.obj", true);
        writeTriangles("feat_triangles.obj", true);
    }
}


// ************************************************************************* //
