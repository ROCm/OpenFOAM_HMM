/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
#include "vectorTools.H"

using namespace Foam::vectorTools;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::insertBoundingPoints()
{
    pointField farPts = geometryToConformTo_.globalBounds().points();

    // Shift corners of bounds relative to origin
    farPts -= geometryToConformTo_.globalBounds().midpoint();

    // Scale the box up
    farPts *= 10.0;

    // Shift corners of bounds back to be relative to midpoint
    farPts += geometryToConformTo_.globalBounds().midpoint();

    limitBounds_ = treeBoundBox(farPts);

    forAll(farPts, fPI)
    {
        insertPoint(farPts[fPI], Vb::vtFar);
    }
}


void Foam::conformalVoronoiMesh::createEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    label edgeI = edHit.index();

    extendedFeatureEdgeMesh::edgeStatus edStatus = feMesh.getEdgeStatus(edgeI);

    switch (edStatus)
    {
        case extendedFeatureEdgeMesh::EXTERNAL:
        {
            createExternalEdgePointGroup(feMesh, edHit, pts, indices, types);
            break;
        }
        case extendedFeatureEdgeMesh::INTERNAL:
        {
            createInternalEdgePointGroup(feMesh, edHit, pts, indices, types);
            break;
        }
        case extendedFeatureEdgeMesh::FLAT:
        {
            createFlatEdgePointGroup(feMesh, edHit, pts, indices, types);
            break;
        }
        case extendedFeatureEdgeMesh::OPEN:
        {
            createOpenEdgePointGroup(feMesh, edHit, pts, indices, types);
            break;
        }
        case extendedFeatureEdgeMesh::MULTIPLE:
        {
            createMultipleEdgePointGroup(feMesh, edHit, pts, indices, types);
            break;
        }
        case extendedFeatureEdgeMesh::NONE:
        {
            break;
        }
    }
}


void Foam::conformalVoronoiMesh::createExternalEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    const Foam::point& edgePt = edHit.hitPoint();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    // As this is an external edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    if (areParallel(nA, nB))
    {
        // The normals are nearly parallel, so this is too sharp a feature to
        // conform to.

        return;
    }

    // Normalised distance of reference point from edge point
    vector refVec((nA + nB)/(1 + (nA & nB)));

    if (magSqr(refVec) > sqr(5.0))
    {
        // Limit the size of the conformation
        ppDist *= 5.0/mag(refVec);

        // Pout<< nl << "createExternalEdgePointGroup limit "
        //     << "edgePt " << edgePt << nl
        //     << "refVec " << refVec << nl
        //     << "mag(refVec) " << mag(refVec) << nl
        //     << "ppDist " << ppDist << nl
        //     << "nA " << nA << nl
        //     << "nB " << nB << nl
        //     << "(nA & nB) " << (nA & nB) << nl
        //     << endl;
    }

    // Convex. So refPt will be inside domain and hence a master point
    Foam::point refPt = edgePt - ppDist*refVec;

    // Result when the points are eventually inserted.
    // Add number_of_vertices() at insertion of first vertex to all numbers:
    // pt           index type
    // refPt        0     1
    // reflectedA   1     0
    // reflectedB   2     0

    // Insert the master point pairing the the first slave

    pts.append(refPt);
    indices.append(0);
    types.append(1);

    // Insert the slave points by reflecting refPt in both faces.
    // with each slave refering to the master

    Foam::point reflectedA = refPt + 2*ppDist*nA;
    pts.append(reflectedA);
    indices.append(0);
    types.append(-1);

    Foam::point reflectedB = refPt + 2*ppDist*nB;
    pts.append(reflectedB);
    indices.append(0);
    types.append(-2);
}


void Foam::conformalVoronoiMesh::createInternalEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    const Foam::point& edgePt = edHit.hitPoint();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    // As this is an external edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    if (areParallel(nA, nB))
    {
        // The normals are nearly parallel, so this is too sharp a feature to
        // conform to.

        return;
    }

    // Normalised distance of reference point from edge point
    vector refVec((nA + nB)/(1 + (nA & nB)));

    if (magSqr(refVec) > sqr(5.0))
    {
        // Limit the size of the conformation
        ppDist *= 5.0/mag(refVec);

        // Pout<< nl << "createInternalEdgePointGroup limit "
        //     << "edgePt " << edgePt << nl
        //     << "refVec " << refVec << nl
        //     << "mag(refVec) " << mag(refVec) << nl
        //     << "ppDist " << ppDist << nl
        //     << "nA " << nA << nl
        //     << "nB " << nB << nl
        //     << "(nA & nB) " << (nA & nB) << nl
        //     << endl;
    }

    // Concave. master and reflected points inside the domain.
    Foam::point refPt = edgePt - ppDist*refVec;

    // Generate reflected master to be outside.
    Foam::point reflMasterPt = refPt + 2*(edgePt - refPt);

    // Reflect reflMasterPt in both faces.
    Foam::point reflectedA = reflMasterPt - 2*ppDist*nA;

    Foam::point reflectedB = reflMasterPt - 2*ppDist*nB;

    scalar totalAngle =
        radToDeg(constant::mathematical::pi + radAngleBetween(nA, nB));

    // Number of quadrants the angle should be split into
    int nQuads = int(totalAngle/cvMeshControls().maxQuadAngle()) + 1;

    // The number of additional master points needed to obtain the
    // required number of quadrants.
    int nAddPoints = min(max(nQuads - 2, 0), 2);

    // index of reflMaster
    label reflectedMaster = 2 + nAddPoints;

    // Add number_of_vertices() at insertion of first vertex to all numbers:
    // Result for nAddPoints 1 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     2
    // reflectedB   1     2
    // reflMasterPt 2     0

    // Result for nAddPoints 1 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     3
    // reflectedB   1     3
    // refPt        2     3
    // reflMasterPt 3     0

    // Result for nAddPoints 2 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     4
    // reflectedB   1     4
    // reflectedAa  2     4
    // reflectedBb  3     4
    // reflMasterPt 4     0

    // Master A is inside.
    pts.append(reflectedA);
    indices.append(0);
    types.append(reflectedMaster--);

    // Master B is inside.
    pts.append(reflectedB);
    indices.append(0);
    types.append(reflectedMaster--);

    if (nAddPoints == 1)
    {
        // One additinal point is the reflection of the slave point,
        // i.e. the original reference point
        pts.append(refPt);
        indices.append(0);
        types.append(reflectedMaster--);
    }
    else if (nAddPoints == 2)
    {
        Foam::point reflectedAa = refPt + ppDist*nB;
        pts.append(reflectedAa);
        indices.append(0);
        types.append(reflectedMaster--);

        Foam::point reflectedBb = refPt + ppDist*nA;
        pts.append(reflectedBb);
        indices.append(0);
        types.append(reflectedMaster--);
    }

    // Slave is outside.
    pts.append(reflMasterPt);
    indices.append(0);
    types.append(-(nAddPoints + 2));
}


void Foam::conformalVoronoiMesh::createFlatEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    const Foam::point& edgePt = edHit.hitPoint();

    const scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    // As this is a flat edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    // Average normal to remove any bias to one side, although as this
    // is a flat edge, the normals should be essentially the same
    const vector n = 0.5*(nA + nB);

    // Direction along the surface to the control point, sense of edge
    // direction not important, as +s and -s can be used because this
    // is a flat edge
    vector s = ppDist*(feMesh.edgeDirections()[edHit.index()] ^ n);

    createPointPair(ppDist, edgePt + s, n, pts, indices, types);
    createPointPair(ppDist, edgePt - s, n, pts, indices, types);
}


void Foam::conformalVoronoiMesh::createOpenEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    Info<< "NOT INSERTING OPEN EDGE POINT GROUP, NOT IMPLEMENTED" << endl;
}


void Foam::conformalVoronoiMesh::createMultipleEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    Info<< "NOT INSERTING MULTIPLE EDGE POINT GROUP, NOT IMPLEMENTED" << endl;
}


void Foam::conformalVoronoiMesh::reinsertBoundingPoints()
{
    tmp<pointField> farPts = limitBounds_.points();

    forAll(farPts(), fPI)
    {
        insertPoint(farPts()[fPI], Vb::vtFar);
    }
}


void Foam::conformalVoronoiMesh::reinsertFeaturePoints(bool distribute)
{
    Info<< nl << "Reinserting stored feature points" << endl;

    label preReinsertionSize(number_of_vertices());

    if (distribute)
    {
        DynamicList<Foam::point> pointsToInsert;
        DynamicList<label> indices;
        DynamicList<label> types;

        for
        (
            List<Vb>::iterator vit=featureVertices_.begin();
            vit != featureVertices_.end();
            ++vit
        )
        {
            pointsToInsert.append(topoint(vit->point()));
            indices.append(vit->index());
            types.append(vit->type());
        }

        // Insert distributed points
        insertPoints(pointsToInsert, indices, types, true);

        // Save points in new distribution
        featureVertices_.clear();
        featureVertices_.setSize(pointsToInsert.size());

        forAll(pointsToInsert, pI)
        {
            featureVertices_[pI] =
                Vb(toPoint(pointsToInsert[pI]), indices[pI], types[pI]);
        }
    }
    else
    {
        for
        (
            List<Vb>::iterator vit=featureVertices_.begin();
            vit != featureVertices_.end();
            ++vit
        )
        {
            // Assuming that all of the reinsertions are pair points, and that
            // the index and type are relative, i.e. index 0 and type relative
            // to it.
            insertPoint
            (
                vit->point(),
                vit->index() + number_of_vertices(),
                vit->type() + number_of_vertices()
            );
        }
    }

    Info<< "    Reinserted "
        << returnReduce
        (
            label(number_of_vertices()) - preReinsertionSize,
            sumOp<label>()
        )
        << " vertices" << endl;
}


void Foam::conformalVoronoiMesh::createConvexFeaturePoints
(
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    const PtrList<extendedFeatureEdgeMesh>& feMeshes
    (
        geometryToConformTo_.features()
    );

    forAll(feMeshes, i)
    {
        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);

        for
        (
            label ptI = feMesh.convexStart();
            ptI < feMesh.concaveStart();
            ptI++
        )
        {
            const Foam::point& featPt = feMesh.points()[ptI];

            if (!positionOnThisProc(featPt))
            {
                continue;
            }

            const vectorField& featPtNormals = feMesh.featurePointNormals(ptI);
            const scalar ppDist = - pointPairDistance(featPt);

            vector cornerNormal = sum(featPtNormals);
            cornerNormal /= mag(cornerNormal);

            Foam::point internalPt = featPt + ppDist*cornerNormal;

            // Result when the points are eventually inserted (example n = 4)
            // Add number_of_vertices() at insertion of first vertex to all
            // numbers:
            // pt           index type
            // internalPt   0     1
            // externalPt0  1     0
            // externalPt1  2     0
            // externalPt2  3     0
            // externalPt3  4     0

            pts.append(internalPt);
            indices.append(0);
            types.append(1);

            label internalPtIndex = -1;

            forAll(featPtNormals, nI)
            {
                const vector& n = featPtNormals[nI];

                plane planeN = plane(featPt, n);

                Foam::point externalPt =
                    internalPt + 2.0 * planeN.distance(internalPt) * n;

                pts.append(externalPt);
                indices.append(0);
                types.append(internalPtIndex--);
            }
        }
    }
}


void Foam::conformalVoronoiMesh::createConcaveFeaturePoints
(
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    const PtrList<extendedFeatureEdgeMesh>& feMeshes
    (
        geometryToConformTo_.features()
    );

    forAll(feMeshes, i)
    {
        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);

        for
        (
            label ptI = feMesh.concaveStart();
            ptI < feMesh.mixedStart();
            ptI++
        )
        {
            const Foam::point& featPt = feMesh.points()[ptI];

            if (!positionOnThisProc(featPt))
            {
                continue;
            }

            const vectorField& featPtNormals = feMesh.featurePointNormals(ptI);
            const scalar ppDist = pointPairDistance(featPt);

            vector cornerNormal = sum(featPtNormals);
            cornerNormal /= mag(cornerNormal);

            Foam::point externalPt = featPt + ppDist*cornerNormal;

            label externalPtIndex = featPtNormals.size();

            // Result when the points are eventually inserted (example n = 5)
            // Add number_of_vertices() at insertion of first vertex to all
            // numbers:
            // pt           index type
            // internalPt0  0     5
            // internalPt1  1     5
            // internalPt2  2     5
            // internalPt3  3     5
            // internalPt4  4     5
            // externalPt   5     4

            forAll(featPtNormals, nI)
            {
                const vector& n = featPtNormals[nI];

                const plane planeN = plane(featPt, n);

                const Foam::point internalPt =
                    externalPt - 2.0 * planeN.distance(externalPt) * n;

                pts.append(internalPt);
                indices.append(0);
                types.append(externalPtIndex--);
            }

            pts.append(externalPt);
            indices.append(0);
            types.append(-1);
        }
    }
}


void Foam::conformalVoronoiMesh::createMixedFeaturePoints
(
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    if (cvMeshControls().mixedFeaturePointPPDistanceCoeff() < 0)
    {
        Info<< nl << "Skipping specialised handling for mixed feature points"
            << endl;
        return;
    }

    const PtrList<extendedFeatureEdgeMesh>& feMeshes
    (
        geometryToConformTo_.features()
    );

    forAll(feMeshes, i)
    {
        const extendedFeatureEdgeMesh& feMesh = feMeshes[i];
        const labelListList& pointsEdges = feMesh.pointEdges();
        const pointField& points = feMesh.points();

        for
        (
            label ptI = feMesh.mixedStart();
            ptI < feMesh.nonFeatureStart();
            ptI++
        )
        {
            const labelList& pEds = pointsEdges[ptI];

            pointFeatureEdgesTypes pFEdgeTypes(ptI);

            const List<extendedFeatureEdgeMesh::edgeStatus> allEdStat
                = calcPointFeatureEdgesTypes(feMesh, pEds, pFEdgeTypes);

            bool specialisedSuccess = false;

            if (cvMeshControls().specialiseFeaturePoints())
            {
                specialisedSuccess = createSpecialisedFeaturePoint
                (
                    feMesh, pEds, pFEdgeTypes, allEdStat, ptI,
                    pts, indices, types
                );
            }

            if (!specialisedSuccess)
            {
                // Specialisations available for some mixed feature points.  For
                // non-specialised feature points, inserting mixed internal and
                // external edge groups at feature point.

                // Skipping unsupported mixed feature point types
                bool skipEdge = false;

                forAll(pEds, e)
                {
                    const label edgeI = pEds[e];

                    const extendedFeatureEdgeMesh::edgeStatus edStatus
                        = feMesh.getEdgeStatus(edgeI);

                    if
                    (
                        edStatus == extendedFeatureEdgeMesh::OPEN
                     || edStatus == extendedFeatureEdgeMesh::MULTIPLE
                    )
                    {
                        Info<< "Edge type " << edStatus
                            << " found for mixed feature point " << ptI
                            << ". Not supported."
                            << endl;

                        skipEdge = true;
                    }
                }

                if (skipEdge)
                {
                    Info<< "Skipping point " << ptI << nl << endl;

                    continue;
                }

                const Foam::point& pt = points[ptI];

                if (!positionOnThisProc(pt))
                {
                    continue;
                }

                const scalar edgeGroupDistance = mixedFeaturePointDistance(pt);

                forAll(pEds, e)
                {
                    const label edgeI = pEds[e];

                    const Foam::point edgePt =
                        pt + edgeGroupDistance*feMesh.edgeDirection(edgeI, ptI);

                    const pointIndexHit edgeHit(true, edgePt, edgeI);

                    createEdgePointGroup(feMesh, edgeHit, pts, indices, types);
                }
            }
        }
    }
}


void Foam::conformalVoronoiMesh::insertFeaturePoints()
{
    Info<< nl << "Conforming to feature points" << endl;

    DynamicList<Foam::point> pts;
    DynamicList<label> indices;
    DynamicList<label> types;

    const label preFeaturePointSize = number_of_vertices();

    createConvexFeaturePoints(pts, indices, types);

    createConcaveFeaturePoints(pts, indices, types);

    createMixedFeaturePoints(pts, indices, types);

    // Insert the created points, distributing to the appropriate processor
    insertPoints(pts, indices, types, true);

    if (cvMeshControls().objOutput())
    {
        writePoints("featureVertices.obj", pts);
    }

    label nFeatureVertices = number_of_vertices() - preFeaturePointSize;

    if (Pstream::parRun())
    {
        reduce(nFeatureVertices, sumOp<label>());
    }

    if (nFeatureVertices > 0)
    {
        Info<< "    Inserted " << nFeatureVertices
            << " feature vertices" << endl;
    }

    featureVertices_.clear();
    featureVertices_.setSize(pts.size());

    forAll(pts, pI)
    {
        featureVertices_[pI] = Vb(toPoint(pts[pI]), indices[pI], types[pI]);
    }

    constructFeaturePointLocations();
}


void Foam::conformalVoronoiMesh::constructFeaturePointLocations()
{
    DynamicList<Foam::point> ftPtLocs;

    const PtrList<extendedFeatureEdgeMesh>& feMeshes
    (
        geometryToConformTo_.features()
    );

    forAll(feMeshes, i)
    {
        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);

        if (cvMeshControls().mixedFeaturePointPPDistanceCoeff() < 0)
        {
            // Ignoring mixed feature points
            for
            (
                label ptI = feMesh.convexStart();
                ptI < feMesh.mixedStart();
                ptI++
            )
            {
                ftPtLocs.append(feMesh.points()[ptI]);
            }
        }
        else
        {
            for
            (
                label ptI = feMesh.convexStart();
                ptI < feMesh.nonFeatureStart();
                ptI++
            )
            {
                ftPtLocs.append(feMesh.points()[ptI]);
            }
        }
    }

    featurePointLocations_.transfer(ftPtLocs);
}
