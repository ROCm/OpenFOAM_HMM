/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "surfaceIntersection.H"
#include "triSurfaceSearch.H"
#include "OFstream.H"
#include "labelPairHashes.H"
#include "triSurface.H"
#include "pointIndexHit.H"
#include "mergePoints.H"
#include "plane.H"
#include "edgeIntersections.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(surfaceIntersection, 0);
}

const Foam::scalar Foam::surfaceIntersection::defaultTolerance = 1e-3;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Checks if there exists a special topological situation that causes
// edge and the face it hit not to be recognized.
//
// For now if the face shares a point with the edge
bool Foam::surfaceIntersection::excludeEdgeHit
(
    const triSurface& surf,
    const label edgeI,
    const label facei,
    const scalar
)
{
    const triSurface::FaceType& f = surf.localFaces()[facei];
    const edge& e = surf.edges()[edgeI];

    forAll(f, fp)
    {
        if (f[fp] == e.start() || f[fp] == e.end())
        {
            return true;
        }
    }

// {
//        // Get edge vector
//        vector eVec = e.vec(surf.localPoints());
//        eVec /= mag(eVec) + VSMALL;
//
//        const labelList& eLabels = surf.faceEdges()[facei];
//
//        // Get edge vector of 0th edge of face
//        vector e0Vec = surf.edges()[eLabels[0]].vec(surf.localPoints());
//        e0Vec /= mag(e0Vec) + VSMALL;
//
//        vector n = e0Vec ^ eVec;
//
//        if (mag(n) < SMALL)
//        {
//            // e0 is aligned with e. Choose next edge of face.
//            vector e1Vec = surf.edges()[eLabels[1]].vec(surf.localPoints());
//            e1Vec /= mag(e1Vec) + VSMALL;
//
//            n = e1Vec ^ eVec;
//
//            if (mag(n) < SMALL)
//            {
//                // Problematic triangle. Two edges aligned with edgeI. Give
//                // up.
//                return true;
//            }
//        }
//
//        // Check if same as faceNormal
//        if (mag(n & surf.faceNormals()[facei]) > 1-tol)
//        {
//
//            Pout<< "edge:" << e << "  face:" << facei
//                << "  e0Vec:" << e0Vec << "  n:" << n
//                << "  normalComponent:" << (n & surf.faceNormals()[facei])
//                << "  tol:" << tol << endl;
//
//            return true;
//        }
//        else
//        {
//            return false;
//        }
// }

    return false;
}


//// Find intersection of plane with edges of hitFacei. Returns
//// - edgeI
//// - intersection point
//Foam::pointIndexHit Foam::surfaceIntersection::faceEdgeIntersection
//(
//    const triSurface& surf,
//    const label hitFacei,
//
//    const vector& n,
//    const point& eStart,
//    const point& eEnd
//)
//{
//    pointIndexHit pInter;
//
//    const pointField& points = surf.points();
//
//    const triSurface::FaceType& f = surf.localFaces()[hitFacei];
//
//    // Plane for intersect test.
//    plane pl(eStart, n);
//
//    forAll(f, fp)
//    {
//        label fp1 = f.fcIndex(fp);
//
//        const point& start = points[f[fp]];
//        const point& end = points[f[fp1]];
//
//        vector eVec(end - start);
//
//        scalar s = pl.normalIntersect(start, eVec);
//
//        if (s < 0 || s > 1)
//        {
//            pInter.setPoint(start + s*eVec);
//
//            // Check if is correct one: orientation walking
//            //  eStart - eEnd - hitPoint should be opposite n
//            vector n2(triPointRef(start, end, pInter.hitPoint()).normal());
//
//            Pout<< "plane normal:" << n
//                << "  start:" << start << "  end:" << end
//                << "  hit at:" << pInter.hitPoint()
//                << "  resulting normal:" << n2 << endl;
//
//            if ((n2 & n) < 0)
//            {
//                pInter.setHit();
//
//                // Find corresponding edge between f[fp] f[fp1]
//                label edgeI =
//                    meshTools::findEdge
//                    (
//                        surf.edges(),
//                        surf.faceEdges()[hitFacei],
//                        f[fp],
//                        f[fp1]
//                    );
//
//                pInter.setIndex(edgeI);
//
//                return pInter;
//            }
//        }
//    }
//
//    FatalErrorInFunction
//        << "Did not find intersection of plane " << pl
//        << " with edges of face " << hitFacei << " verts:" << f
//        << abort(FatalError);
//
//    return pInter;
//}


void Foam::surfaceIntersection::storeIntersection
(
    const enum originatingType cutFrom,
    const labelList& facesA,
    const label faceB,
    DynamicList<edge>& allCutEdges,
    DynamicList<point>& allCutPoints
)
{
    // Our lookup for two faces - populate with faceB (invariant)
    // Normally always have face from the first surface as first element
    labelPair twoFaces(faceB, faceB);

    const label pointId = allCutPoints.size()-1;

    forAll(facesA, facesAI)
    {
        const label faceA = facesA[facesAI];

        switch (cutFrom)
        {
            case surfaceIntersection::FIRST:
            {
                // faceA from 1st, faceB from 2nd
                twoFaces.first() = faceA;
                break;
            }
            case surfaceIntersection::SECOND:
            {
                // faceA from 2nd, faceB from 1st
                twoFaces.second() = faceA;
                break;
            }
            case surfaceIntersection::SELF:
            {
                // Lookup should be commutativity - use sorted order
                if (faceA < faceB)
                {
                    twoFaces.first()  = faceA;
                    twoFaces.second() = faceB;
                }
                else
                {
                    twoFaces.first()  = faceB;
                    twoFaces.second() = faceA;
                }
                break;
            }
        }

        labelPairLookup::const_iterator iter = facePairToVertex_.find(twoFaces);

        if (iter == facePairToVertex_.end())
        {
            // New intersection. Store face-face intersection.
            facePairToVertex_.insert(twoFaces, pointId);
        }
        else
        {
            // Second occurrence of surf1-surf2 intersection.
            // Or rather the face on surf1 intersects a face on
            // surface2 twice -> we found edge.

            // Check whether perhaps degenerate
            const point& prevHit = allCutPoints[*iter];
            const point& thisHit = allCutPoints.last();

            if (mag(prevHit - thisHit) < SMALL)
            {
                WarningInFunction
                    << "Encountered degenerate edge between face "
                    << twoFaces[0] << " on first surface and face "
                    << twoFaces[1] << " on second surface" << endl
                    << "Point on first surface:" << prevHit << endl
                    << "Point on second surface:" << thisHit << endl
                    << endl;
            }
            else
            {
                allCutEdges.append(edge(*iter, pointId));

                // Record intersection of faces on surf
                facePairToEdge_.insert(twoFaces, allCutEdges.size()-1);
            }
        }
    }
}


// Classify cut of edge of surface1 with surface2:
// 1- point of edge hits point on surface2
// 2- edge pierces point on surface2
// 3- point of edge hits edge on surface2
// 4- edge pierces edge on surface2
// 5- point of edge hits face on surface2
// 6- edge pierces face on surface2
//
// Note that handling of 2 and 4 should be the same but with surface1 and
// surface2 reversed.
void Foam::surfaceIntersection::classifyHit
(
    const triSurface& surf1,
    const scalarField& surf1PointTol,
    const triSurface& surf2,
    const enum originatingType cutFrom,
    const label edgeI,
    const pointIndexHit& pHit,

    DynamicList<edge>& allCutEdges,
    DynamicList<point>& allCutPoints,
    List<DynamicList<label>>& surfEdgeCuts
)
{
    const edge& e = surf1.edges()[edgeI];

    const labelList& facesA = surf1.edgeFaces()[edgeI];

    // Label of face on surface2 edgeI intersected
    label surf2Facei = pHit.index();

    // Classify point on surface2

    const triSurface::FaceType& f2 = surf2.localFaces()[surf2Facei];
    const pointField& surf2Pts = surf2.localPoints();

    label nearType, nearLabel;

    f2.nearestPointClassify(pHit.hitPoint(), surf2Pts, nearType, nearLabel);

    // Classify points on edge of surface1
    const label edgeEnd =
        classify
        (
            surf1PointTol[e.start()],
            surf1PointTol[e.end()],
            pHit.hitPoint(),
            e,
            surf1.localPoints()
        );

    if (nearType == triPointRef::POINT)
    {
        if (edgeEnd >= 0)
        {
            // 1. Point hits point. Do nothing.
            if (debug & 2)
            {
                Pout<< pHit.hitPoint() << " is surf1:"
                    << " end point of edge " << e
                    << " surf2: vertex " << f2[nearLabel]
                    << " coord:" << surf2Pts[f2[nearLabel]] << endl;
            }
        }
        else
        {
            // 2. Edge hits point. Cut edge with new point.
            if (debug & 2)
            {
                Pout<< pHit.hitPoint() << " is surf1:"
                    << " somewhere on edge " << e
                    << " surf2: vertex " << f2[nearLabel]
                    << " coord:" << surf2Pts[f2[nearLabel]] << endl;
            }

            allCutPoints.append(pHit.hitPoint());
            surfEdgeCuts[edgeI].append(allCutPoints.size()-1);

            const labelList& facesB = surf2.pointFaces()[f2[nearLabel]];

            forAll(facesB, faceBI)
            {
                storeIntersection
                (
                    cutFrom,
                    facesA,
                    facesB[faceBI],
                    allCutEdges,
                    allCutPoints
                );
            }
        }
    }
    else if (nearType == triPointRef::EDGE)
    {
        if (edgeEnd >= 0)
        {
            // 3. Point hits edge. Do nothing on this side. Reverse
            // is handled by 2 (edge hits point)
            label edge2I = getEdge(surf2, surf2Facei, nearLabel);
            const edge& e2 = surf2.edges()[edge2I];

            if (debug & 2)
            {
                Pout<< pHit.hitPoint() << " is surf1:"
                    << " end point of edge " << e
                    << " surf2: edge " << e2
                    << " coords:" << surf2Pts[e2.start()]
                    << surf2Pts[e2.end()] << endl;
            }
        }
        else
        {
            // 4. Edge hits edge.

            // Cut edge with new point (creates duplicates when
            // doing the surf2 with surf1 intersection but these
            // are merged later on)

            label edge2I = getEdge(surf2, surf2Facei, nearLabel);
            const edge& e2 = surf2.edges()[edge2I];

            if (debug & 2)
            {
                Pout<< pHit.hitPoint() << " is surf1:"
                    << " somewhere on edge " << e
                    << " surf2: edge " << e2
                    << " coords:" << surf2Pts[e2.start()]
                    << surf2Pts[e2.end()] << endl;
            }

            allCutPoints.append(pHit.hitPoint());
            surfEdgeCuts[edgeI].append(allCutPoints.size()-1);

            // edge hits all faces on surf2 connected to the edge

            if (cutFrom == surfaceIntersection::FIRST)
            {
                // edge-edge intersection is symmetric, store only
                // once.
                // edge hits all faces on surf2 connected to the
                // edge

                const labelList& facesB = surf2.edgeFaces()[edge2I];

                forAll(facesB, faceBI)
                {
                    storeIntersection
                    (
                        cutFrom,
                        facesA,
                        facesB[faceBI],
                        allCutEdges,
                        allCutPoints
                    );
                }
            }
        }
    }
    else
    {
        if (edgeEnd >= 0)
        {
            // 5. Point hits face. Do what? Introduce
            // point & triangulation in face?
            if (debug & 2)
            {
                Pout<< pHit.hitPoint() << " is surf1:"
                    << " end point of edge " << e
                    << " surf2: face " << surf2Facei
                    << endl;
            }

            //
            // Look exactly at what side (of surf2) edge is. Leave out ones on
            // inside of surf2 (i.e. on opposite side of normal)
            //

            // Vertex on/near surf2; vertex away from surf2
            const label nearVert  = (edgeEnd == 0 ?  e.start() : e.end());
            const label otherVert = (edgeEnd == 0 ?  e.end() : e.start());

            const point& nearPt  = surf1.localPoints()[nearVert];
            const point& otherPt = surf1.localPoints()[otherVert];

            if (debug)
            {
                Pout
                    << pHit.hitPoint() << " is surf1:"
                    << " end point of edge " << e << " coord:"
                    << nearPt
                    << " surf2: face " << surf2Facei << endl;
            }

            const vector eVec = otherPt - nearPt;

            if ((surf2.faceNormals()[surf2Facei] & eVec) > 0)
            {
                // otherVert on outside of surf2

                // Shift hitPoint a bit along edge.
                //point hitPt = nearPt + 0.1*eVec;
                point hitPt = nearPt;

                if (debug & 2)
                {
                    Pout<< "Shifted " << pHit.hitPoint()
                        << " to " << hitPt
                        << " along edge:" << e
                        << " coords:" << surf1.localPoints()[e.start()]
                        << surf1.localPoints()[e.end()] << endl;
                }

                // Reclassify as normal edge-face pierce (see below)

                allCutPoints.append(hitPt);
                surfEdgeCuts[edgeI].append(allCutPoints.size()-1);

                // edge hits single face only
                storeIntersection
                (
                    cutFrom,
                    facesA,
                    surf2Facei,
                    allCutEdges,
                    allCutPoints
                );
            }
            else
            {
                if (debug & 2)
                {
                    Pout<< "Discarding " << pHit.hitPoint()
                        << " since edge " << e << " on inside of surf2."
                        << " surf2 normal:" << surf2.faceNormals()[surf2Facei]
                        << endl;
                }
            }
        }
        else
        {
            // 6. Edge pierces face. 'Normal' situation.
            if (debug & 2)
            {
                Pout<< pHit.hitPoint() << " is surf1:"
                    << " somewhere on edge " << e
                    << " surf2: face " << surf2Facei
                    << endl;
            }

            // edgeI intersects surf2. Store point.
            allCutPoints.append(pHit.hitPoint());
            surfEdgeCuts[edgeI].append(allCutPoints.size()-1);

            // edge hits single face only
            storeIntersection
            (
                cutFrom,
                facesA,
                surf2Facei,
                allCutEdges,
                allCutPoints
            );
        }
    }
    if (debug & 2)
    {
        Pout<< endl;
    }
}


// Cut all edges of surf1 with surf2. Sets
// - cutPoints          : coordinates of cutPoints
// - cutEdges           : newly created edges between cutPoints
// - facePairToVertex   : hash from face1I and face2I to cutPoint
// - facePairToEdge     : hash from face1I and face2I to cutEdge
// - surfEdgeCuts       : gives for each edge the cutPoints
//                        (in order from start to end)
//
void Foam::surfaceIntersection::doCutEdges
(
    const triSurface& surf1,
    const triSurfaceSearch& querySurf2,
    const enum originatingType cutFrom,

    DynamicList<edge>& allCutEdges,
    DynamicList<point>& allCutPoints,
    List<DynamicList<label>>& surfEdgeCuts
)
{
    const scalar oldTol = intersection::setPlanarTol(planarTol_);

    const pointField& surf1Pts = surf1.localPoints();

    // Calculate local (to point) tolerance based on min edge length.
    scalarField surf1PointTol(surf1Pts.size());

    forAll(surf1PointTol, pointi)
    {
        surf1PointTol[pointi] =
            intersection::planarTol()
          * minEdgeLen(surf1, pointi);
    }

    const indexedOctree<treeDataPrimitivePatch<triSurface>>& searchTree
        = querySurf2.tree();

    if (cutFrom == surfaceIntersection::SELF)
    {
        // An edge may intersect multiple faces
        // - mask out faces that have already been hit before trying again
        // - never intersect with faces attached to the edge itself
        DynamicList<label> maskFaces(8);

        treeDataTriSurface::findAllIntersectOp
            allIntersectOp(searchTree, maskFaces);

        forAll(surf1.edges(), edgeI)
        {
            const edge& e = surf1.edges()[edgeI];
            const vector edgeVec = e.vec(surf1Pts);

            // Extend start/end by tolerance - ensures cleaner cutting
            const point ptStart =
                surf1Pts[e.start()] - surf1PointTol[e.start()]*edgeVec;
            const point ptEnd =
                surf1Pts[e.end()]   + surf1PointTol[e.end()]*edgeVec;

            // Never intersect with faces attached to the edge itself
            maskFaces = surf1.edgeFaces()[edgeI];
            while (true)
            {
                pointIndexHit pHit = searchTree.findLine
                (
                    ptStart,
                    ptEnd,
                    allIntersectOp
                );

                if (!pHit.hit())
                {
                    break;
                }

                maskFaces.append(pHit.index());

                classifyHit
                (
                    surf1,
                    surf1PointTol,
                    surf1,
                    cutFrom,
                    edgeI,
                    pHit,

                    allCutEdges,
                    allCutPoints,
                    surfEdgeCuts
                );
            }
        }
    }
    else
    {
        const triSurface& surf2 = querySurf2.surface();

        forAll(surf1.edges(), edgeI)
        {
            const edge& e = surf1.edges()[edgeI];
            const vector edgeVec = e.vec(surf1Pts);

            const point tolVec = intersection::planarTol()*(edgeVec);
            const scalar tolDim = mag(tolVec);

            point ptStart = surf1Pts[e.start()];
            const point ptEnd = surf1Pts[e.end()];

            bool doTrack = false;
            do
            {
                pointIndexHit pHit = searchTree.findLine(ptStart, ptEnd);

                if (!pHit.hit())
                {
                    break;
                }

                classifyHit
                (
                    surf1,
                    surf1PointTol,
                    surf2,
                    cutFrom,
                    edgeI,
                    pHit,

                    allCutEdges,
                    allCutPoints,
                    surfEdgeCuts
                );

                if (tolDim > 0.0)
                {
                    if (mag(pHit.hitPoint() - ptEnd) < tolDim)
                    {
                        doTrack = false;
                    }
                    else
                    {
                        // continue tracking a bit further on
                        ptStart = pHit.hitPoint() + tolVec;
                        doTrack = true;
                    }
                }
            }
            while (doTrack);  // execute at least once
        }
    }

    intersection::setPlanarTol(oldTol);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceIntersection::surfaceIntersection()
:
    planarTol_(defaultTolerance),
    cutPoints_(0),
    cutEdges_(0),
    facePairToVertex_(0),
    facePairToEdge_(0),
    surf1EdgeCuts_(0),
    surf2EdgeCuts_(0)
{}


Foam::surfaceIntersection::surfaceIntersection
(
    const triSurfaceSearch& query1,
    const triSurfaceSearch& query2,
    const scalar planarTol
)
:
    planarTol_(planarTol),
    cutPoints_(0),
    cutEdges_(0),
    facePairToVertex_(2*max(query1.surface().size(), query2.surface().size())),
    facePairToEdge_(2*max(query1.surface().size(), query2.surface().size())),
    surf1EdgeCuts_(0),
    surf2EdgeCuts_(0)
{
    const triSurface& surf1 = query1.surface();
    const triSurface& surf2 = query2.surface();

    //
    // Cut all edges of surf1 with surf2.
    //
    if (debug)
    {
        Pout<< "Cutting surf1 edges" << endl;
    }


    DynamicList<edge> allCutEdges(surf1.nEdges()/20);
    DynamicList<point> allCutPoints(surf1.nPoints()/20);


    // From edge to cut index on surface1
    List<DynamicList<label>> edgeCuts1(query1.surface().nEdges());

    // 1st surf (for labelPair order)
    doCutEdges
    (
        surf1,
        query2,
        surfaceIntersection::FIRST,
        allCutEdges,
        allCutPoints,
        edgeCuts1
    );
    // Transfer to straight labelListList
    transfer(edgeCuts1, surf1EdgeCuts_);


    //
    // Cut all edges of surf2 with surf1.
    //
    if (debug)
    {
        Pout<< "Cutting surf2 edges" << endl;
    }

    // From edge to cut index
    List<DynamicList<label>> edgeCuts2(query2.surface().nEdges());

    // 2nd surf (for labelPair order)
    doCutEdges
    (
        surf2,
        query1,
        surfaceIntersection::SECOND,
        allCutEdges,
        allCutPoints,
        edgeCuts2
    );

    // Transfer to straight label(List)List
    transfer(edgeCuts2, surf2EdgeCuts_);
    cutEdges_.transfer(allCutEdges);
    cutPoints_.transfer(allCutPoints);

    if (debug)
    {
        Pout<< "surfaceIntersection : Intersection generated:"
            << endl
            << "    points:" << cutPoints_.size() << endl
            << "    edges :" << cutEdges_.size() << endl;

        Pout<< "surfaceIntersection : Writing intersection to intEdges.obj"
            << endl;

        OFstream intStream("intEdges.obj");
        writeOBJ(cutPoints_, cutEdges_, intStream);

        // Dump all cut edges to files
        Pout<< "Dumping cut edges of surface1 to surf1EdgeCuts.obj" << endl;
        OFstream edge1Stream("surf1EdgeCuts.obj");
        writeIntersectedEdges(surf1, surf1EdgeCuts_, edge1Stream);

        Pout<< "Dumping cut edges of surface2 to surf2EdgeCuts.obj" << endl;
        OFstream edge2Stream("surf2EdgeCuts.obj");
        writeIntersectedEdges(surf2, surf2EdgeCuts_, edge2Stream);
    }
}


Foam::surfaceIntersection::surfaceIntersection
(
    const triSurfaceSearch& query1,
    const scalar planarTol
)
:
    planarTol_(planarTol),
    cutPoints_(0),
    cutEdges_(0),
    facePairToVertex_(2*query1.surface().size()),
    facePairToEdge_(2*query1.surface().size()),
    surf1EdgeCuts_(0),
    surf2EdgeCuts_(0)
{
    const triSurface& surf1 = query1.surface();

    //
    // Cut all edges of surf1 with surf1 itself.
    //
    if (debug)
    {
        Pout<< "Cutting surf1 edges" << endl;
    }

    DynamicList<edge> allCutEdges;
    DynamicList<point> allCutPoints;

    // From edge to cut index on surface1
    List<DynamicList<label>> edgeCuts1(query1.surface().nEdges());

    // self-intersection
    doCutEdges
    (
        surf1,
        query1,
        surfaceIntersection::SELF,
        allCutEdges,
        allCutPoints,
        edgeCuts1
    );

    // Transfer to straight label(List)List
    transfer(edgeCuts1, surf1EdgeCuts_);
    cutEdges_.transfer(allCutEdges);
    cutPoints_.transfer(allCutPoints);

    // Shortcut.
    if (cutPoints_.empty() && cutEdges_.empty())
    {
        if (debug)
        {
            Pout<< "Empty intersection" << endl;
        }
        return;
    }

    if (debug)
    {
        Pout<< "surfaceIntersection : Intersection generated and compressed:"
            << endl
            << "    points:" << cutPoints_.size() << endl
            << "    edges :" << cutEdges_.size() << endl;


        Pout<< "surfaceIntersection : Writing intersection to intEdges.obj"
            << endl;

        OFstream intStream("intEdges.obj");
        writeOBJ(cutPoints_, cutEdges_, intStream);

        // Dump all cut edges to files
        Pout<< "Dumping cut edges of surface1 to surf1EdgeCuts.obj" << endl;
        OFstream edge1Stream("surf1EdgeCuts.obj");
        writeIntersectedEdges(surf1, surf1EdgeCuts_, edge1Stream);
    }
}


Foam::surfaceIntersection::surfaceIntersection
(
    const triSurface& surf1,
    const edgeIntersections& intersections1,
    const triSurface& surf2,
    const edgeIntersections& intersections2
)
:
    cutPoints_(0),
    cutEdges_(0),
    facePairToVertex_(2*max(surf1.size(), surf2.size())),
    facePairToEdge_(2*max(surf1.size(), surf2.size())),
    surf1EdgeCuts_(0),
    surf2EdgeCuts_(0)
{

    // All intersection Pout (so for both surfaces)
    DynamicList<edge> allCutEdges((surf1.nEdges() + surf2.nEdges())/20);
    DynamicList<point> allCutPoints((surf1.nPoints() + surf2.nPoints())/20);


    // Cut all edges of surf1 with surf2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "Storing surf1 intersections" << endl;
    }

    {
        // From edge to cut index on surface1
        List<DynamicList<label>> edgeCuts1(surf1.nEdges());

        forAll(intersections1, edgeI)
        {
            const List<pointIndexHit>& intersections = intersections1[edgeI];

            forAll(intersections, i)
            {
                const pointIndexHit& pHit = intersections[i];

                // edgeI intersects surf2. Store point.
                allCutPoints.append(pHit.hitPoint());
                edgeCuts1[edgeI].append(allCutPoints.size()-1);

                storeIntersection
                (
                    surfaceIntersection::FIRST,
                    surf1.edgeFaces()[edgeI],
                    pHit.index(),               // surf2Facei
                    allCutEdges,
                    allCutPoints
                );
            }
        }

        // Transfer to straight labelListList
        transfer(edgeCuts1, surf1EdgeCuts_);
    }



    // Cut all edges of surf2 with surf1
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "Storing surf2 intersections" << endl;
    }

    {
        // From edge to cut index on surface2
        List<DynamicList<label>> edgeCuts2(surf2.nEdges());

        forAll(intersections2, edgeI)
        {
            const List<pointIndexHit>& intersections = intersections2[edgeI];

            forAll(intersections, i)
            {
                const pointIndexHit& pHit = intersections[i];

                // edgeI intersects surf1. Store point.
                allCutPoints.append(pHit.hitPoint());
                edgeCuts2[edgeI].append(allCutPoints.size()-1);

                storeIntersection
                (
                    surfaceIntersection::SECOND,
                    surf2.edgeFaces()[edgeI],
                    pHit.index(),               // surf2Facei
                    allCutEdges,
                    allCutPoints
                );
            }
        }

        // Transfer to surf2EdgeCuts_ (straight labelListList)
        transfer(edgeCuts2, surf2EdgeCuts_);
    }


    // Transfer to straight label(List)List
    cutEdges_.transfer(allCutEdges);
    cutPoints_.transfer(allCutPoints);


    if (debug)
    {
        Pout<< "surfaceIntersection : Intersection generated:"
            << endl
            << "    points:" << cutPoints_.size() << endl
            << "    edges :" << cutEdges_.size() << endl;

        Pout<< "surfaceIntersection : Writing intersection to intEdges.obj"
            << endl;

        OFstream intStream("intEdges.obj");
        writeOBJ(cutPoints_, cutEdges_, intStream);

        // Dump all cut edges to files
        Pout<< "Dumping cut edges of surface1 to surf1EdgeCuts.obj" << endl;
        OFstream edge1Stream("surf1EdgeCuts.obj");
        writeIntersectedEdges(surf1, surf1EdgeCuts_, edge1Stream);

        Pout<< "Dumping cut edges of surface2 to surf2EdgeCuts.obj" << endl;
        OFstream edge2Stream("surf2EdgeCuts.obj");
        writeIntersectedEdges(surf2, surf2EdgeCuts_, edge2Stream);
    }


    // Debugging stuff
    {
        // Check all facePairToVertex is used.
        labelHashSet usedPoints;

        forAllConstIter(labelPairLookup, facePairToEdge_, iter)
        {
            label edgeI = iter();

            const edge& e = cutEdges_[edgeI];

            usedPoints.insert(e[0]);
            usedPoints.insert(e[1]);
        }

        forAllConstIter(labelPairLookup, facePairToVertex_, iter)
        {
            label pointi = iter();

            if (!usedPoints.found(pointi))
            {
                WarningInFunction
                    << "Problem: cut point:" << pointi
                    << " coord:" << cutPoints_[pointi]
                    << " not used by any edge" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::surfaceIntersection::cutPoints() const
{
    return cutPoints_;
}


const Foam::edgeList& Foam::surfaceIntersection::cutEdges() const
{
    return cutEdges_;
}


const Foam::labelPairLookup& Foam::surfaceIntersection::facePairToVertex() const
{
    return facePairToVertex_;
}


const Foam::labelPairLookup& Foam::surfaceIntersection::facePairToEdge() const
{
    return facePairToEdge_;
}


const Foam::labelListList& Foam::surfaceIntersection::edgeCuts
(
    const bool isFirstSurf
) const
{
    if (isFirstSurf)
    {
        return surf1EdgeCuts_;
    }
    else
    {
        return surf2EdgeCuts_;
    }
}


const Foam::labelListList& Foam::surfaceIntersection::surf1EdgeCuts() const
{
    return surf1EdgeCuts_;
}


const Foam::labelListList& Foam::surfaceIntersection::surf2EdgeCuts() const
{
    return surf2EdgeCuts_;
}


// ************************************************************************* //
