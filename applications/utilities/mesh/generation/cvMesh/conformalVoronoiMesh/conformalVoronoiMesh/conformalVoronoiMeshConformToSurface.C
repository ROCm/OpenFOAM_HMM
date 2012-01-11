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
#include "backgroundMeshDecomposition.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::conformToSurface()
{
    reconformationMode reconfMode = reconformationControl();

    if (reconfMode == rmNone)
    {
        // Reinsert stored surface conformation
        reinsertSurfaceConformation();

        buildParallelInterface("move_" + runTime_.timeName());
    }
    else
    {
        // Rebuild, insert and store new surface conformation
        buildSurfaceConformation(reconfMode);

        if (distributeBackground())
        {
            // distributeBackground has destroyed all referred vertices, so the
            // parallel interface needs to be rebuilt.
            buildParallelInterface("rebuild");

            // Use storeSizesAndAlignments with no feed points because all
            // background points may have been distributed.
            storeSizesAndAlignments();
        }

        // Do not store the surface conformation until after it has been
        // (potentially) redistributed.
        storeSurfaceConformation();
    }

    // reportSurfaceConformationQuality();
}


Foam::conformalVoronoiMesh::reconformationMode
Foam::conformalVoronoiMesh::reconformationControl() const
{
    if (!runTime_.run())
    {
        Info<< nl << "Rebuilding surface conformation for final output"
            << endl;

        return rmFine;
    }
    else if
    (
        runTime_.timeIndex()
      % cvMeshControls().surfaceConformationRebuildFrequency() == 0
    )
    {
        Info<< nl << "Rebuilding surface conformation for more iterations"
            << endl;

        return rmCoarse;
    }

    return rmNone;
}


void Foam::conformalVoronoiMesh::buildSurfaceConformation
(
    reconformationMode reconfMode
)
{
    timeCheck("Start buildSurfaceConformation");

    if (reconfMode == rmCoarse)
    {
        Info<< nl << "Build coarse surface conformation" << endl;
    }
    else if (reconfMode == rmFine)
    {
        Info<< nl << "Build fine surface conformation" << endl;
    }
    else if (reconfMode == rmNone)
    {
        WarningIn("buildSurfaceConformation(reconformationMode reconfMode)")
            << "reconformationMode rmNone specified, not building conformation"
            << endl;

        return;
    }
    else
    {
        WarningIn("buildSurfaceConformation(reconformationMode reconfMode)")
            << "Unknown reconformationMode " << reconfMode
            << " not building conformation" << endl;

        return;
    }

    // Initialise containers to store the edge conformation locations
    DynamicList<Foam::point> newEdgeLocations;

    pointField existingEdgeLocations(0);

    autoPtr<indexedOctree<treeDataPoint> > edgeLocationTree;

    // Initialise the edgeLocationTree
    buildEdgeLocationTree(edgeLocationTree, existingEdgeLocations);

    label initialTotalHits = 0;


    // Surface protrusion conformation is done in two steps.
    // 1. the dual edges (of all internal vertices) can stretch to
    //    'infinity' so any intersection would be badly behaved. So
    //    just find the nearest point on the geometry and insert point
    //    pairs.
    // Now most of the surface conformation will be done with some
    // residual protrusions / incursions.
    // 2. find any segments of dual edges outside the geometry. Shoot
    //    ray from Delaunay vertex to middle of this segment and introduce
    //    point pairs. This will handle e.g.

    // protruding section of face:
    //
    //     internal
    // \             /
    // -+-----------+-- boundary
    //   \         /
    //     --------
    //
    // Shoot ray and find intersection with outside segment (x) and
    // introduce pointpair (..)
    //
    //        |
    // \      .      /
    // -+-----|-----+-- boundary
    //   \    .    /
    //     ---x----



    // Initial surface protrusion conformation - nearest surface point
    {
        const scalar edgeSearchDistCoeffSqr =
            cvMeshControls().edgeSearchDistCoeffSqrInitial(reconfMode);

        const scalar surfacePtReplaceDistCoeffSqr =
            cvMeshControls().surfacePtReplaceDistCoeffSqrInitial(reconfMode);

        DynamicList<pointIndexHit> surfaceHits;
        DynamicList<label> hitSurfaces;

        DynamicList<pointIndexHit> featureEdgeHits;
        DynamicList<label> featureEdgeFeaturesHit;

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->internalPoint())
            {
                const Foam::point vert = topoint(vit->point());

                DynamicList<pointIndexHit> surfHitList;
                DynamicList<label> hitSurfaceList;

                if
                (
                    dualCellSurfaceAllIntersections
                    (
                        vit,
                        surfHitList,
                        hitSurfaceList
                    )
                )
                {
                    // This used to be just before this if statement.
                    // Moved because a point is only near the boundary if
                    // the dual cell intersects the surface.
                    vit->setNearBoundary();

                    // meshTools::writeOBJ(Pout, vert);
                    // meshTools::writeOBJ(Pout, surfHit.hitPoint());
                    // Pout<< "l cr0 cr1" << endl;

                    addSurfaceAndEdgeHits
                    (
                        vit,
                        vert,
                        surfHitList,
                        hitSurfaceList,
                        surfacePtReplaceDistCoeffSqr,
                        edgeSearchDistCoeffSqr,
                        surfaceHits,
                        hitSurfaces,
                        featureEdgeHits,
                        featureEdgeFeaturesHit,
                        newEdgeLocations,
                        existingEdgeLocations,
                        edgeLocationTree
                    );
                }
            }
        }

        label nVerts = number_of_vertices();
        label nSurfHits = surfaceHits.size();
        label nFeatEdHits = featureEdgeHits.size();

        if (Pstream::parRun())
        {
            reduce(nVerts, sumOp<label>());
            reduce(nSurfHits, sumOp<label>());
            reduce(nFeatEdHits, sumOp<label>());
        }

        Info<< nl << "Initial conformation" << nl
            << "    Number of vertices " << nVerts << nl
            << "    Number of surface hits " << nSurfHits << nl
            << "    Number of edge hits " << nFeatEdHits
            << endl;

        insertSurfacePointPairs
        (
            surfaceHits,
            hitSurfaces,
            "surfaceConformationLocations_initial.obj"
        );

        insertEdgePointGroups
        (
            featureEdgeHits,
            featureEdgeFeaturesHit,
            "edgeConformationLocations_initial.obj"
        );

        // Filter small edges at the boundary by inserting surface point pairs
//        for
//        (
//            Delaunay::Finite_cells_iterator cit = finite_cells_begin();
//            cit != finite_cells_end();
//            cit++
//        )
//        {
//            const Foam::point dualVertex = topoint(dual(cit));
//
//
//        }
//        for
//        (
//            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
//            vit != finite_vertices_end();
//            vit++
//        )
//        {
//            if (vit->nearBoundary())
//            {
//                std::list<Facet> faces;
//                incident_facets(vit, std::back_inserter(faces));
//
//                List<scalar> edgeLengthList(faces.size());
//                scalar totalLength = 0;
//                label count = 0;
//
//                for
//                (
//                    std::list<Facet>::iterator fit=faces.begin();
//                    fit != faces.end();
//                    ++fit
//                )
//                {
//                    if
//                    (
//                        is_infinite(fit->first)
//                     || is_infinite(fit->first->neighbor(fit->second))
//                    )
//                    {
//                        continue;
//                    }
//
//                    const Point dE0 = dual(fit->first);
//                    const Point dE1 = dual(fit->first->neighbor(fit->second));
//
//                    const Segment s(dE0, dE1);
//
//                    const scalar length = Foam::sqrt(s.squared_length());
//
//                    edgeLengthList[count++] = length;
//
//                    totalLength += length;
//
//                    //Info<< length << " / " << totalLength << endl;
//                }
//
//                const scalar averageLength = totalLength/edgeLengthList.size();
//
//                forAll(edgeLengthList, eI)
//                {
//                    const scalar& el = edgeLengthList[eI];
//
//                    if (el < 0.1*averageLength)
//                    {
//                        //Info<< "SMALL" << endl;
//                    }
//                }
//            }
//        }

        timeCheck("After initial conformation");

        initialTotalHits = nSurfHits + nFeatEdHits;
    }

    // Remember which vertices were referred to each processor so only updates
    // are sent.
    List<labelHashSet> referralVertices(Pstream::nProcs());

    // Store the vertices that have been received and added from each processor
    // already so that there is no attempt to add them more than once.
    List<labelHashSet> receivedVertices(Pstream::nProcs());

    // Build the parallel interface the initial surface conformation
    buildParallelInterface(referralVertices, receivedVertices, true, "initial");

    label iterationNo = 0;

    label maxIterations =
        cvMeshControls().maxConformationIterations(reconfMode);

    scalar iterationToInitialHitRatioLimit =
        cvMeshControls().iterationToInitialHitRatioLimit(reconfMode);

    label hitLimit = label(iterationToInitialHitRatioLimit*initialTotalHits);

    Info<< nl << "Stopping iterations when: " << nl
        <<"    total number of hits drops below "
        << iterationToInitialHitRatioLimit << " of initial hits ("
        << hitLimit << ")" << nl
        << " or " << nl
        << "    maximum number of iterations ("
        << maxIterations << ") is reached"
        << endl;

    // Set totalHits to a large enough positive value to enter the while loop on
    // the first iteration
    label totalHits = initialTotalHits;

    while
    (
        totalHits > 0
     && totalHits >= hitLimit
     && iterationNo < maxIterations
    )
    {
        scalar edgeSearchDistCoeffSqr =
            cvMeshControls().edgeSearchDistCoeffSqrIteration(reconfMode);

        scalar surfacePtReplaceDistCoeffSqr =
            cvMeshControls().surfacePtReplaceDistCoeffSqrIteration(reconfMode);

        DynamicList<pointIndexHit> surfaceHits;
        DynamicList<label> hitSurfaces;

        DynamicList<pointIndexHit> featureEdgeHits;
        DynamicList<label> featureEdgeFeaturesHit;

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            // The initial surface conformation has already identified the
            // nearBoundary set of vertices.  Previously inserted boundary
            // points and referred internal vertices from other processors can
            // also generate protrusions and must be assessed too.
            if
            (
                vit->nearBoundary()
             || vit->ppMaster()
             || vit->referredInternalOrBoundaryPoint()
            )
            {
                const Foam::point vert = topoint(vit->point());

                pointIndexHit surfHit;
                label hitSurface;

                // Find segments of dual face outside the geometry and find the
                // the middle of this
                dualCellLargestSurfaceProtrusion(vit, surfHit, hitSurface);

                if (surfHit.hit())
                {
                    DynamicList<pointIndexHit> tmpPIH;
                    tmpPIH.append(surfHit);

                    DynamicList<label> tmpHS;
                    tmpHS.append(hitSurface);

                    addSurfaceAndEdgeHits
                    (
                        vit,
                        vert,
                        tmpPIH,
                        tmpHS,
                        surfacePtReplaceDistCoeffSqr,
                        edgeSearchDistCoeffSqr,
                        surfaceHits,
                        hitSurfaces,
                        featureEdgeHits,
                        featureEdgeFeaturesHit,
                        newEdgeLocations,
                        existingEdgeLocations,
                        edgeLocationTree
                    );
                }
            }
            else if (vit->ppSlave() || vit->referredExternal())
            {
                const Foam::point vert = topoint(vit->point());

                pointIndexHit surfHit;
                label hitSurface;

                // Detect slave (external vertices) whose dual face incurs
                // into nearby (other than originating) geometry
                dualCellLargestSurfaceIncursion(vit, surfHit, hitSurface);

                if (surfHit.hit())
                {
                    DynamicList<pointIndexHit> tmpPIH;
                    tmpPIH.append(surfHit);

                    DynamicList<label> tmpHS;
                    tmpHS.append(hitSurface);

                    addSurfaceAndEdgeHits
                    (
                        vit,
                        vert,
                        tmpPIH,
                        tmpHS,
                        surfacePtReplaceDistCoeffSqr,
                        edgeSearchDistCoeffSqr,
                        surfaceHits,
                        hitSurfaces,
                        featureEdgeHits,
                        featureEdgeFeaturesHit,
                        newEdgeLocations,
                        existingEdgeLocations,
                        edgeLocationTree
                    );
                }
            }
        }

        label nVerts = number_of_vertices();
        label nSurfHits = surfaceHits.size();
        label nFeatEdHits = featureEdgeHits.size();

        if (Pstream::parRun())
        {
            reduce(nVerts, sumOp<label>());
            reduce(nSurfHits, sumOp<label>());
            reduce(nFeatEdHits, sumOp<label>());
        }

        Info<< nl << "Conformation iteration " << iterationNo << nl
            << "    Number of vertices " << nVerts << nl
            << "    Number of surface hits " << nSurfHits << nl
            << "    Number of edge hits " << nFeatEdHits
            << endl;

        totalHits = nSurfHits + nFeatEdHits;

        if (totalHits > 0)
        {
            insertSurfacePointPairs
            (
                surfaceHits,
                hitSurfaces,
                fileName
                (
                    "surfaceConformationLocations_" + name(iterationNo) + ".obj"
                )
            );

            insertEdgePointGroups
            (
                featureEdgeHits,
                featureEdgeFeaturesHit,
                "edgeConformationLocations_" + name(iterationNo) + ".obj"
            );
        }

        timeCheck("Conformation iteration " + name(iterationNo));

        // Update the parallel interface
        buildParallelInterface
        (
            referralVertices,
            receivedVertices,
            false,
            name(iterationNo)
        );

        iterationNo++;

        if (iterationNo == maxIterations)
        {
            WarningIn("conformalVoronoiMesh::conformToSurface()")
                << "Maximum surface conformation iterations ("
                << maxIterations <<  ") reached." << endl;
        }

        if (totalHits < hitLimit)
        {
            Info<< nl << "Total hits (" << totalHits
                << ") less than limit (" << hitLimit
                << "), stopping iterations" << endl;
        }
    }
}


bool Foam::conformalVoronoiMesh::dualCellSurfaceAnyIntersection
(
    const Delaunay::Finite_vertices_iterator& vit
) const
{
    std::list<Facet> facets;
    incident_facets(vit, std::back_inserter(facets));

    for
    (
        std::list<Facet>::iterator fit=facets.begin();
        fit != facets.end();
        ++fit
    )
    {
        if
        (
            is_infinite(fit->first)
         || is_infinite(fit->first->neighbor(fit->second))
        )
        {
            continue;
        }

        Foam::point dE0 = topoint(dual(fit->first));
        Foam::point dE1 = topoint(dual(fit->first->neighbor(fit->second)));

        if (Pstream::parRun())
        {
            Foam::point& a = dE0;
            Foam::point& b = dE1;

            bool inProc = clipLineToProc(a, b);

            // Check for the edge passing through a surface
            if
            (
                inProc
             && geometryToConformTo_.findSurfaceAnyIntersection(a, b)
            )
            {
                // Pout<< "# findSurfaceAnyIntersection" << endl;
                // meshTools::writeOBJ(Pout, a);
                // meshTools::writeOBJ(Pout, b);
                // Pout<< "l cr0 cr1" << endl;

                return true;
            }
        }
        else
        {
            if (geometryToConformTo_.findSurfaceAnyIntersection(dE0, dE1))
            {
                return true;
            }
        }
    }

    return false;
}


bool Foam::conformalVoronoiMesh::dualCellSurfaceAllIntersections
(
    const Delaunay::Finite_vertices_iterator& vit,
    DynamicList<pointIndexHit>& infoList,
    DynamicList<label>& hitSurfaceList
) const
{
    bool flagIntersection = false;

    std::list<Facet> facets;
    incident_facets(vit, std::back_inserter(facets));

    for
    (
        std::list<Facet>::iterator fit=facets.begin();
        fit != facets.end();
        ++fit
    )
    {
        if
        (
            is_infinite(fit->first)
         || is_infinite(fit->first->neighbor(fit->second))
        )
        {
            continue;
        }

        const Foam::point dE0 = topoint(dual(fit->first));

        const Foam::point dE1
            = topoint(dual(fit->first->neighbor(fit->second)));

        pointIndexHit info;
        label hitSurface = -1;

        geometryToConformTo_.findSurfaceAnyIntersection
        (
            dE0,
            dE1,
            info,
            hitSurface
        );

        if (info.hit())
        {
            flagIntersection = true;

            vectorField norm(1);

            allGeometry_[hitSurface].getNormal
            (
                List<pointIndexHit>(1, info),
                norm
            );

            const vector& n = norm[0];

            const Foam::point vertex = topoint(vit->point());

            const plane p(info.hitPoint(), n);

            const plane::ray r(vertex, n);

            const scalar d = p.normalIntersect(r);

            const Foam::point newPoint = vertex + d*n;

            geometryToConformTo_.findSurfaceAnyIntersection
            (
                vertex,
                vertex + 1.1*(newPoint - vertex),
                info,
                hitSurface
            );

            bool rejectPoint = false;

            if (info.hit())
            {
                if (!infoList.empty())
                {
                    forAll(infoList, hitI)
                    {
                        // Reject point if the point is already added
                        if (infoList[hitI].index() == info.index())
                        {
                            rejectPoint = true;
                            break;
                        }

                        const Foam::point& p = infoList[hitI].hitPoint();

                        const scalar separationDistance
                            = mag(p - info.hitPoint());

                        const scalar minSepDist
                            = cvMeshControls().removalDistCoeff()
                             *targetCellSize(p);

                        // Reject the point if it is too close to another
                        // surface point.
                        // Could merge?
                        if (separationDistance < minSepDist)
                        {
                            rejectPoint = true;
                            break;
                        }
                    }
                }
            }

            // The normal ray from the vertex will not always result in a hit
            // because another surface may be in the way.
            if (!rejectPoint && info.hit())
            {
                infoList.append(info);
                hitSurfaceList.append(hitSurface);
            }
        }
    }

    return flagIntersection;
}


bool Foam::conformalVoronoiMesh::clipLineToProc
(
    Foam::point& a,
    Foam::point& b
) const
{
    bool intersects = false;

    if (decomposition_().positionOnThisProcessor(a))
    {
        // a is inside this processor

        if (decomposition_().positionOnThisProcessor(b))
        {
            // both a and b inside, no clip required
            intersects = true;
        }
        else
        {
            // b is outside, clip b to bounding box.

            pointIndexHit info = decomposition_().findLine(b, a);

            if (info.hit())
            {
                intersects = true;
                b = info.hitPoint();
            }
        }
    }
    else
    {
        // a is outside this processor

        if (decomposition_().positionOnThisProcessor(b))
        {
            // b is inside this processor, clip a to processor

            pointIndexHit info = decomposition_().findLine(a, b);

            if (info.hit())
            {
                intersects = true;
                a = info.hitPoint();
            }
        }
        else
        {
            // both a and b outside, but they can still intersect the processor

            pointIndexHit info = decomposition_().findLine(a, b);

            if (info.hit())
            {
                intersects = true;
                a = info.hitPoint();

                info = decomposition_().findLine(b, a);

                if (info.hit())
                {
                    b = info.hitPoint();
                }
            }
        }
    }

    return intersects;
}


void Foam::conformalVoronoiMesh::buildParallelInterface
(
    const word& outputName
)
{
    List<labelHashSet> referralVertices(Pstream::nProcs());
    List<labelHashSet> receivedVertices(Pstream::nProcs());

    buildParallelInterface
    (
        referralVertices,
        receivedVertices,
        true,
        outputName
    );
}


void Foam::conformalVoronoiMesh::buildParallelInterface
(
    List<labelHashSet>& referralVertices,
    List<labelHashSet>& receivedVertices,
    bool initialEdgeReferral,
    const word& outputName
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    Info<< nl << "Parallel interface construction" << endl;

    timeCheck("Before buildParallelInterface");

    // Hard coded switch, can be turned on for testing and debugging purposes -
    // all vertices will be referred to all processors, use with caution for
    // big cases.
    bool allPointReferral = false;

    if (allPointReferral)
    {
        buildParallelInterfaceAll
        (
            referralVertices,
            receivedVertices,
            outputName
        );

        timeCheck("After buildParallelInterfaceAll");
    }

    if (initialEdgeReferral)
    {
        // Used as an initial pass to localise the vertex referring - find
        // vertices whose dual edges pierce nearby processor volumes and refer
        // them to establish a sensible boundary interface region before
        // running a circumsphere assessment.

        buildParallelInterfaceIntersection
        (
            referralVertices,
            receivedVertices,
            outputName
        );

        timeCheck("After buildParallelInterfaceIntersection");
    }

    buildParallelInterfaceInfluence
    (
        referralVertices,
        receivedVertices,
        outputName
    );

    timeCheck("After buildParallelInterface");
}


void Foam::conformalVoronoiMesh::buildParallelInterfaceAll
(
    List<labelHashSet>& referralVertices,
    List<labelHashSet>& receivedVertices,
    const word& outputName
)
{
    // Refer all points to all processors

    DynamicList<Foam::point> parallelAllPoints;
    DynamicList<label> targetProcessor;
    DynamicList<label> parallelAllIndices;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (!vit->farPoint())
        {
            for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {
                if (procI == Pstream::myProcNo())
                {
                    continue;
                }

                label vIndex = vit->index();

                // Using the hashSet to ensure that each vertex is
                // only referred once to each processor
                if (!referralVertices[procI].found(vIndex))
                {
                    referralVertices[procI].insert(vIndex);

                    parallelAllPoints.append
                    (
                        topoint(vit->point())
                    );

                    targetProcessor.append(procI);

                    if (vit->internalOrBoundaryPoint())
                    {
                        parallelAllIndices.append(vIndex);
                    }
                    else
                    {
                        parallelAllIndices.append(-vIndex);
                    }
                }
            }
        }
    }

    referVertices
    (
        targetProcessor,
        parallelAllPoints,
        parallelAllIndices,
        receivedVertices,
        "all",
        outputName
    );
}

void Foam::conformalVoronoiMesh::buildParallelInterfaceIntersection
(
    List<labelHashSet>& referralVertices,
    List<labelHashSet>& receivedVertices,
    const word& outputName
)
{
    DynamicList<Foam::point> parallelIntersectionPoints;
    DynamicList<label> targetProcessor;
    DynamicList<label> parallelIntersectionIndices;

    // End points of dual edges.  Some of these values will not be used,
    // i.e. for edges with non-real vertices.
    DynamicList<Foam::point> dE0;
    DynamicList<Foam::point> dE1;

    PackedBoolList testFacetIntersection(number_of_facets(), false);

    // Index outer (all) Delaunauy facets for whether they are potential
    // intersections, index (inner) the list of tests an results.
    label fIInner = 0;
    label fIOuter = 0;

    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const Cell_handle c2(c1->neighbor(fit->second));

        // If either Delaunay cell at the end of the Dual edge is infinite,
        // skip.
        if (!is_infinite(c1) && !is_infinite(c2))
        {
            // The Delaunauy cells at either end of the dual edge need to be
            // real, i.e. all vertices form part of the internal or boundary
            // definition
            if
            (
                c1->internalOrBoundaryDualVertex()
             && c2->internalOrBoundaryDualVertex()
            )
            {
                Foam::point a = topoint(dual(c1));
                Foam::point b = topoint(dual(c2));

                // Only if the dual edge cuts the boundary of this processor is
                // it going to be counted.
                if (decomposition_().findLineAny(a, b).hit())
                {
                    dE0.append(a);
                    dE1.append(b);

                    testFacetIntersection[fIOuter] = true;
                }
            }
        }

        fIOuter++;
    }

    // Preform intersections in both directions, as there is no sense
    // associated with the Dual edge
    List<List<pointIndexHit> > intersectionForward(intersectsProc(dE0, dE1));
    List<List<pointIndexHit> > intersectionReverse(intersectsProc(dE1, dE0));

    // Reset counter
    fIOuter = 0;

    // Relying on the order of iteration of facets being the same as before
    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const Cell_handle c2(c1->neighbor(fit->second));

        // Pre-tested facet intersection potential
        if (testFacetIntersection[fIOuter])
        {
            const Foam::point a = dE0[fIInner];
            const Foam::point b = dE1[fIInner];

            scalar hitDistSqr = GREAT;
            label closestHitProc = -1;
            scalar closestHitDistSqr = GREAT;

            // Find the closest intersection with the other background meshes
            // of the other processors in each direction, finding the closest.
            forAll(intersectionForward[fIInner], iFI)
            {
                const pointIndexHit& info = intersectionForward[fIInner][iFI];

                if (info.hit())
                {
                    // Line was a -> hit
                    hitDistSqr = magSqr(a - info.hitPoint());
                }

                if (hitDistSqr < closestHitDistSqr)
                {
                    closestHitProc = info.index();
                    closestHitDistSqr = hitDistSqr;
                }
            }

            forAll(intersectionReverse[fIInner], iRI)
            {
                const pointIndexHit& info = intersectionReverse[fIInner][iRI];

                if (info.hit())
                {
                    // Line was b -> hit
                    hitDistSqr = magSqr(b - info.hitPoint());
                }

                if (hitDistSqr < closestHitDistSqr)
                {
                    closestHitProc = info.index();
                    closestHitDistSqr = hitDistSqr;
                }
            }

            if (closestHitProc >= 0)
            {
                // This dual edge pierces a processor, refer all vertices from
                // both Delaunauy cells to it.

                for (int i = 0; i < 4; i++)
                {
                    Vertex_handle v = c1->vertex(i);

                    label vIndex = v->index();

                    if (v->farPoint() || v->referred())
                    {
                        continue;
                    }

                    // Using the hashSet to ensure that each vertex is only
                    // referred once to each processor
                    if (!referralVertices[closestHitProc].found(vIndex))
                    {
                        referralVertices[closestHitProc].insert(vIndex);

                        parallelIntersectionPoints.append(topoint(v->point()));

                        targetProcessor.append(closestHitProc);

                        if (v->internalOrBoundaryPoint())
                        {
                            parallelIntersectionIndices.append(vIndex);
                        }
                        else
                        {
                            parallelIntersectionIndices.append(-vIndex);
                        }

                        // Pout<< "Refer "
                        //     << parallelIntersectionPoints.last()
                        //     << " " << parallelIntersectionIndices.last()
                        //     << " " << closestHitProc
                        //     << endl;
                    }
                }

                for (int i = 0; i < 4; i++)
                {
                    Vertex_handle v = c2->vertex(i);

                    label vIndex = v->index();

                    if (v->farPoint() || v->referred())
                    {
                        continue;
                    }

                    // Using the hashSet to ensure that each vertex is only
                    // referred once to each processor
                    if (!referralVertices[closestHitProc].found(vIndex))
                    {
                        referralVertices[closestHitProc].insert(vIndex);

                        parallelIntersectionPoints.append(topoint(v->point()));

                        targetProcessor.append(closestHitProc);

                        if (v->internalOrBoundaryPoint())
                        {
                            parallelIntersectionIndices.append(vIndex);
                        }
                        else
                        {
                            parallelIntersectionIndices.append(-vIndex);
                        }

                        // Pout<< "Refer "
                        //     << parallelIntersectionPoints.last()
                        //     << " " << parallelIntersectionIndices.last()
                        //     << " " << closestHitProc
                        //     << endl;
                    }
                }
            }

            fIInner++;
        }

        fIOuter++;
    }

    referVertices
    (
        targetProcessor,
        parallelIntersectionPoints,
        parallelIntersectionIndices,
        receivedVertices,
        "intersection",
        outputName
    );
}


void Foam::conformalVoronoiMesh::buildParallelInterfaceInfluence
(
    List<labelHashSet>& referralVertices,
    List<labelHashSet>& receivedVertices,
    const word& outputName
)
{
    DynamicList<Foam::point> parallelInfluencePoints;
    DynamicList<label> targetProcessor;
    DynamicList<label> parallelInfluenceIndices;

    // Some of these values will not be used, i.e. for non-real cells
    DynamicList<Foam::point> circumcentre;
    DynamicList<scalar> circumradiusSqr;

    PackedBoolList testCellInfluence(number_of_cells(), false);

    // Index outer (all) Delaunauy cells for whether they are potential
    // overlaps, index (inner) the list of tests an results.
    label cIInner = 0;
    label cIOuter = 0;

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        // Assess the influence of the circumsphere of each Delaunay cell with
        // the defining volumes for all processors.  Any processor touched by
        // the circumsphere requires all points of the cell to be referred to
        // it.

        // The Delaunay cells to assess have to be real, i.e. all vertices form
        // part of the internal or any part of the boundary definition
        if (cit->real())
        {
            Foam::point cc(topoint(dual(cit)));

            scalar crSqr
            (
                magSqr(cc - topoint(cit->vertex(0)->point()))
            );

            // Only if the circumsphere overlaps the boundary of this processor
            // is there a chance of it overlapping others
            if (decomposition_().overlapsThisProcessor(cc, sqr(1.01)*crSqr))
            {
                circumcentre.append(cc);
                circumradiusSqr.append(crSqr);

                testCellInfluence[cIOuter] = true;
            }
        }

        cIOuter++;
    }

    // Increasing the circumspheres to increase the overlaps and compensate for
    // floating point errors missing some referrals
    labelListList circumsphereOverlaps
    (
        overlapsProc(circumcentre, sqr(1.01)*circumradiusSqr)
    );

    // Reset counters
    cIOuter = 0;

    // Relying on the order of iteration of cells being the same as before
    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        // Pre-tested circumsphere potential influence
        if (testCellInfluence[cIOuter])
        {
            const labelList& citOverlaps = circumsphereOverlaps[cIInner];

            forAll(citOverlaps, cOI)
            {
                label procI = citOverlaps[cOI];

                for (int i = 0; i < 4; i++)
                {
                    Vertex_handle v = cit->vertex(i);

                    label vIndex = v->index();

                    if (v->farPoint() || v->referred())
                    {
                        continue;
                    }

                    // Using the hashSet to ensure that each vertex is only
                    // referred once to each processor
                    if (!referralVertices[procI].found(vIndex))
                    {
                        referralVertices[procI].insert(vIndex);

                        parallelInfluencePoints.append(topoint(v->point()));

                        targetProcessor.append(procI);

                        if (v->internalOrBoundaryPoint())
                        {
                            parallelInfluenceIndices.append(vIndex);
                        }
                        else
                        {
                            parallelInfluenceIndices.append(-vIndex);
                        }
                    }
                }
            }

            cIInner++;
        }

        cIOuter++;
    }

    referVertices
    (
        targetProcessor,
        parallelInfluencePoints,
        parallelInfluenceIndices,
        receivedVertices,
        "influence",
        outputName
    );
}


void Foam::conformalVoronoiMesh::referVertices
(
    const DynamicList<label>& targetProcessor,
    DynamicList<Foam::point>& parallelPoints,
    DynamicList<label>& parallelIndices,
    List<labelHashSet>& receivedVertices,
    const word& stageName,
    const word& outputName
)
{
    timeCheck("Start of referVertices " + stageName);

    if (cvMeshControls().objOutput())
    {
        writePoints
        (
            "parallel_" + stageName +"_pointsToSend_" + outputName + ".obj",
            parallelPoints
        );
    }

    mapDistribute pointMap = backgroundMeshDecomposition::buildMap
    (
        targetProcessor
    );

    label totalVertices = parallelPoints.size();

    reduce(totalVertices, sumOp<label>());

    pointMap.distribute(parallelPoints);

    pointMap.distribute(parallelIndices);

    if (cvMeshControls().objOutput())
    {
        writePoints
        (
            "parallel_" + stageName +"_pointsReceived_" + outputName + ".obj",
            parallelPoints
        );
    }

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        const labelList& constructMap = pointMap.constructMap()[procI];

        if (constructMap.size())
        {
            forAll(constructMap, i)
            {
                label origIndex = parallelIndices[constructMap[i]];

                if (!receivedVertices[procI].found(origIndex))
                {
                    // For the initial referred vertices, the original
                    // processor is the one that is sending it.
                    label encodedProcI = Vb::encodeProcIndex(procI);

                    // Pout<< "Insert "
                    //     << parallelPoints[constructMap[i]]
                    //     << " " << origIndex
                    //     << " " << procI
                    //     << endl;

                    insertPoint
                    (
                        parallelPoints[constructMap[i]],
                        origIndex,
                        encodedProcI
                    );

                    receivedVertices[procI].insert(origIndex);
                }
            }
        }
    }

    Info<< "    Total " << stageName << " vertices " << totalVertices << endl;

    timeCheck("End of referVertices " + stageName);
}


void Foam::conformalVoronoiMesh::dualCellLargestSurfaceProtrusion
(
    const Delaunay::Finite_vertices_iterator& vit,
    pointIndexHit& surfHitLargest,
    label& hitSurfaceLargest
) const
{
    // Set no-hit data
    surfHitLargest = pointIndexHit();
    hitSurfaceLargest = -1;

    std::list<Facet> facets;
    incident_facets(vit, std::back_inserter(facets));

    const Foam::point vert = topoint(vit->point());

    scalar maxProtrusionDistance = maxSurfaceProtrusion(vert);

    for
    (
        std::list<Facet>::iterator fit = facets.begin();
        fit != facets.end();
        ++fit
    )
    {
        if
        (
            !is_infinite(fit->first)
         && !is_infinite(fit->first->neighbor(fit->second))
        )
        {
            const Foam::point edgeMid =
                0.5
               *(
                    topoint(dual(fit->first))
                  + topoint(dual(fit->first->neighbor(fit->second)))
                );

            pointIndexHit surfHit;
            label hitSurface;

            geometryToConformTo_.findSurfaceAnyIntersection
            (
                vert,
                edgeMid,
                surfHit,
                hitSurface
            );

            if (surfHit.hit())
            {
                vectorField norm(1);

                allGeometry_[hitSurface].getNormal
                (
                    List<pointIndexHit>(1, surfHit),
                    norm
                );

                const vector& n = norm[0];

                const scalar normalProtrusionDistance =
                    (edgeMid - surfHit.hitPoint()) & n;

                if (normalProtrusionDistance > maxProtrusionDistance)
                {
                    surfHitLargest = surfHit;
                    hitSurfaceLargest = hitSurface;

                    maxProtrusionDistance = normalProtrusionDistance;
                }
            }
        }
    }

    // Relying on short-circuit evaluation to not call for hitPoint when this
    // is a miss
    if
    (
        surfHitLargest.hit()
     && !positionOnThisProc(surfHitLargest.hitPoint())
    )
    {
        // A protrusion was identified, but not penetrating on this processor,
        // so set no-hit data and allow the other that should have this point
        // referred to generate it.
        surfHitLargest = pointIndexHit();
        hitSurfaceLargest = -1;
    }
}


void Foam::conformalVoronoiMesh::dualCellLargestSurfaceIncursion
(
    const Delaunay::Finite_vertices_iterator& vit,
    pointIndexHit& surfHitLargest,
    label& hitSurfaceLargest
) const
{
    // Set no-hit data
    surfHitLargest = pointIndexHit();
    hitSurfaceLargest = -1;

    std::list<Facet> facets;
    incident_facets(vit, std::back_inserter(facets));

    Foam::point vert(topoint(vit->point()));

    scalar minIncursionDistance = -maxSurfaceProtrusion(vert);

    for
    (
        std::list<Facet>::iterator fit = facets.begin();
        fit != facets.end();
        ++fit
    )
    {
        if
        (
            !is_infinite(fit->first)
         && !is_infinite(fit->first->neighbor(fit->second))
        )
        {
            Foam::point edgeMid =
                0.5
               *(
                    topoint(dual(fit->first))
                  + topoint(dual(fit->first->neighbor(fit->second)))
                );

            pointIndexHit surfHit;
            label hitSurface;

            geometryToConformTo_.findSurfaceAnyIntersection
            (
                vert,
                edgeMid,
                surfHit,
                hitSurface
            );

            if (surfHit.hit())
            {
                vectorField norm(1);

                allGeometry_[hitSurface].getNormal
                (
                    List<pointIndexHit>(1, surfHit),
                    norm
                );

                const vector& n = norm[0];

                scalar normalIncursionDistance =
                    (edgeMid - surfHit.hitPoint()) & n;

                if (normalIncursionDistance < minIncursionDistance)
                {
                    surfHitLargest = surfHit;
                    hitSurfaceLargest = hitSurface;

                    minIncursionDistance = normalIncursionDistance;

                    // Info<< nl << "# Incursion: " << endl;
                    // meshTools::writeOBJ(Info, vert);
                    // meshTools::writeOBJ(Info, edgeMid);
                    // Info<< "l Na Nb" << endl;
                }
            }
        }
    }

    // Relying on short-circuit evaluation to not call for hitPoint when this
    // is a miss
    if
    (
        surfHitLargest.hit()
     && !positionOnThisProc(surfHitLargest.hitPoint())
    )
    {
        // A protrusion was identified, but not penetrating on this processor,
        // so set no-hit data and allow the other that should have this point
        // referred to generate it.
        surfHitLargest = pointIndexHit();
        hitSurfaceLargest = -1;
    }
}


void Foam::conformalVoronoiMesh::reportProcessorOccupancy()
{
    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->real())
        {
            if (!positionOnThisProc(topoint(vit->point())))
            {
                Pout<< topoint(vit->point()) << " is not on this processor "
                    << endl;
            }
        }
    }
}


void Foam::conformalVoronoiMesh::reportSurfaceConformationQuality()
{
    Info<< nl << "Check surface conformation quality" << endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            Foam::point vert(topoint(vit->point()));
            pointIndexHit surfHit;
            label hitSurface;

            dualCellLargestSurfaceProtrusion(vit, surfHit, hitSurface);

            if (surfHit.hit())
            {
                Pout<< nl << "Residual penetration: " << nl
                    << vit->index() << nl
                    << vit->type() << nl
                    << vit->ppMaster() << nl
                    << "nearFeaturePt "
                    << nearFeaturePt(surfHit.hitPoint()) << nl
                    << vert << nl
                    << surfHit.hitPoint()
                    << endl;
            }
        }
    }

    {
        // Assess close surface points

        setVertexSizeAndAlignment();

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->ppMaster())
            {
                std::list<Vertex_handle> adjacentVertices;

                adjacent_vertices(vit, std::back_inserter(adjacentVertices));

                Foam::point pt = topoint(vit->point());

                // Pout<< nl << "vit: " << vit->index() << " "
                //     << topoint(vit->point())
                //     << endl;

                // Pout<< adjacentVertices.size() << endl;

                for
                (
                    std::list<Vertex_handle>::iterator
                    avit = adjacentVertices.begin();
                    avit != adjacentVertices.end();
                    ++avit
                )
                {
                    Vertex_handle avh = *avit;

                    // The lower indexed vertex will perform the assessment
                    if
                    (
                        avh->ppMaster()
                     && vit->index() < avh->index()
                     && vit->type() != avh->type()
                    )
                    {
                        scalar targetSize = 0.2*averageAnyCellSize(vit, avh);

                        // Pout<< "diff " << mag(pt - topoint(avh->point()))
                        //     << " " << targetSize << endl;

                        if
                        (
                            magSqr(pt - topoint(avh->point()))
                          < sqr(targetSize)
                        )
                        {
                            Pout<< nl << "vit: " << vit->index() << " "
                                << topoint(vit->point())
                                << endl;

                            Pout<< "    adjacent too close: "
                                << avh->index() << " "
                                << topoint(avh->point())
                                << endl;
                        }
                    }
                }
            }
        }
    }
}


void Foam::conformalVoronoiMesh::limitDisplacement
(
    const Delaunay::Finite_vertices_iterator& vit,
    vector& displacement,
    label callCount
) const
{
    callCount++;

    // Do not allow infinite recursion
    if (callCount > 7)
    {
        return;
    }

    Foam::point pt = topoint(vit->point());
    Foam::point dispPt = pt + displacement;

    bool limit = false;

    pointIndexHit surfHit;
    label hitSurface;

    if (!geometryToConformTo_.globalBounds().contains(dispPt))
    {
        // If dispPt is outside bounding box then displacement cuts boundary
        limit = true;
    }
    else if (geometryToConformTo_.findSurfaceAnyIntersection(pt, dispPt))
    {
        // Full surface penetration test
        limit = true;
    }
    else
    {
        // Testing if the displaced position is too close to the surface.
        // Within twice the local surface point pair insertion distance is
        // considered "too close"

        scalar searchDistanceSqr = sqr
        (
            2*vit->targetCellSize()
           *cvMeshControls().pointPairDistanceCoeff()
        );

        geometryToConformTo_.findSurfaceNearest
        (
            dispPt,
            searchDistanceSqr,
            surfHit,
            hitSurface
        );

        if (surfHit.hit())
        {
            limit = true;

            if (magSqr(pt - surfHit.hitPoint()) <= searchDistanceSqr)
            {
                // Cannot limit displacement, point closer than tolerance
                return;
            }
        }
    }

    if (limit)
    {
        // Halve the displacement and call this function again.  Will continue
        // recursively until the displacement is small enough.

        displacement *= 0.5;

        limitDisplacement(vit, displacement, callCount);
    }
}


bool Foam::conformalVoronoiMesh::nearFeatureEdgeLocation
(
    const Foam::point& pt,
    DynamicList<Foam::point>& newEdgeLocations,
    pointField& existingEdgeLocations,
    autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree
) const
{
    scalar exclusionRangeSqr = featureEdgeExclusionDistanceSqr(pt);

    // 0.01 and 1000 determined from speed tests, varying the indexedOctree
    // rebuild frequency and balance of additions to queries.

    if
    (
        newEdgeLocations.size()
     >= max(0.01*existingEdgeLocations.size(), 1000)
    )
    {
        existingEdgeLocations.append(newEdgeLocations);

        buildEdgeLocationTree(edgeLocationTree, existingEdgeLocations);

        newEdgeLocations.clear();
    }
    else
    {
        // Search for the nearest point in newEdgeLocations.
        // Searching here first, because the intention is that the value of
        // newEdgeLocationsSizeLimit should make this faster by design.

        if (min(magSqr(newEdgeLocations - pt)) <= exclusionRangeSqr)
        {
            return true;
        }
    }

    // Searching for the nearest point in existingEdgeLocations using the
    // indexedOctree

    pointIndexHit info = edgeLocationTree().findNearest(pt, exclusionRangeSqr);

    // Average the points...

    return info.hit();
}


void Foam::conformalVoronoiMesh::buildEdgeLocationTree
(
    autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree,
    const pointField& existingEdgeLocations
) const
{
    treeBoundBox overallBb
    (
        geometryToConformTo_.globalBounds().extend(rndGen_, 1e-4)
    );

    overallBb.min() -= Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    edgeLocationTree.reset
    (
        new indexedOctree<treeDataPoint>
        (
            treeDataPoint(existingEdgeLocations),
            overallBb,  // overall search domain
            10,         // max levels
            10.0,       // maximum ratio of cubes v.s. cells
            100.0       // max. duplicity; n/a since no bounding boxes.
        )
    );
}


void Foam::conformalVoronoiMesh::buildSizeAndAlignmentTree() const
{
    if (sizeAndAlignmentLocations_.empty())
    {
        FatalErrorIn("buildSizeAndAlignmentTree()")
            << "sizeAndAlignmentLocations empty, must be populated before "
            << "sizeAndAlignmentTree can be built."
            << exit(FatalError);
    }

    treeBoundBox overallBb
    (
        geometryToConformTo_.globalBounds().extend(rndGen_, 1e-4)
    );

    overallBb.min() -= Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    sizeAndAlignmentTreePtr_.reset
    (
        new indexedOctree<treeDataPoint>
        (
            treeDataPoint(sizeAndAlignmentLocations_),
            overallBb,  // overall search domain
            10,         // max levels
            10.0,       // maximum ratio of cubes v.s. cells
            100.0       // max. duplicity; n/a since no bounding boxes.
        )
    );
}


void Foam::conformalVoronoiMesh::addSurfaceAndEdgeHits
(
    const Delaunay::Finite_vertices_iterator& vit,
    const Foam::point& vert,
    const DynamicList<pointIndexHit>& surfHit,
    const DynamicList<label>& hitSurface,
    scalar surfacePtReplaceDistCoeffSqr,
    scalar edgeSearchDistCoeffSqr,
    DynamicList<pointIndexHit>& surfaceHits,
    DynamicList<label>& hitSurfaces,
    DynamicList<pointIndexHit>& featureEdgeHits,
    DynamicList<label>& featureEdgeFeaturesHit,
    DynamicList<Foam::point>& newEdgeLocations,
    pointField& existingEdgeLocations,
    autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree
) const
{
    bool keepSurfacePoint = true;

    const scalar targetCellSizeSqr = sqr(targetCellSize(vert));

    DynamicList<label> pointsUsed;

    // Place edge points for external edges.
    forAll(surfHit, hI)
    {
        if (surfHit.size() < 2)
        {
            break;
        }

        const pointIndexHit& hitInfo = surfHit[hI];
        const label& hitSurf = hitSurface[hI];

        if (!hitInfo.hit())
        {
            continue;
        }

        vectorField norm(1);

        allGeometry_[hitSurf].getNormal
        (
            List<pointIndexHit>(1, hitInfo),
            norm
        );

        const vector& n = norm[0];

        const plane surfPlane(hitInfo.hitPoint(), n);

        pointsUsed.append(hI);

        // Find the new edge points by finding the distance of the other surface
        // points from the the plane that the first point is on.
        forAll(surfHit, hI2)
        {
            bool alreadyUsed = false;
            forAll(pointsUsed, puI)
            {
                if (hI2 == pointsUsed[puI])
                {
                    alreadyUsed = true;
                }
            }

            if (alreadyUsed == true)
            {
                continue;
            }

            const pointIndexHit& hitInfo2 = surfHit[hI2];

            const label& hitSurf2 = hitSurface[hI2];

            if (!hitInfo2.hit())
            {
                continue;
            }

            const Foam::point& surfPoint = hitInfo2.hitPoint();

            const scalar d = surfPlane.distance(surfPoint);

            vectorField norm2(1);

            allGeometry_[hitSurf2].getNormal
            (
                List<pointIndexHit>(1, hitInfo2),
                norm2
            );

            // If the normals are nearly parallel then ignore.
            if (mag(mag(norm2[0] & n) - 1.0) < SMALL)
            {
               continue;
            }

            const Foam::point newEdgePoint = surfPoint + n*d;

            pointIndexHit edHit;
            label featureHit = -1;

            // This is a bit lazy...
            geometryToConformTo_.findEdgeNearest
            (
                newEdgePoint,
                edgeSearchDistCoeffSqr*targetCellSizeSqr,
                edHit,
                featureHit
            );

            if (edHit.hit())
            {
                if (!nearFeaturePt(edHit.hitPoint()))
                {
                    if
                    (
                        !nearFeatureEdgeLocation
                        (
                            edHit.hitPoint(),
                            newEdgeLocations,
                            existingEdgeLocations,
                            edgeLocationTree
                        )
                    )
                    {
                        featureEdgeHits.append(edHit);
                        featureEdgeFeaturesHit.append(featureHit);
                        newEdgeLocations.append(edHit.hitPoint());
                    }
                }
            }
        }
    }

    forAll(surfHit, sI)
    {
        const pointIndexHit& surfHitI = surfHit[sI];
        const label hitSurfaceI = hitSurface[sI];

        if (!surfHitI.hit())
        {
            continue;
        }

        if (nearFeaturePt(surfHitI.hitPoint()))
        {
            keepSurfacePoint = false;

            // If the triggering Vertex is part of a feature point, allow it to
            // conform to the surface
            if (vit->index() < startOfInternalPoints_)
            {
                surfaceHits.append(surfHitI);

                hitSurfaces.append(hitSurfaceI);
            }
        }

        List<pointIndexHit> edHits;

        labelList featuresHit;

        geometryToConformTo_.findEdgeNearestByType
        (
            surfHitI.hitPoint(),
            edgeSearchDistCoeffSqr*targetCellSizeSqr,
            edHits,
            featuresHit
        );

        // Gather edge locations but do not add them to newEdgeLocations inside
        // the loop as they will prevent nearby edge locations of different
        // types being conformed to.

        DynamicList<Foam::point> currentEdgeLocations;

        forAll(edHits, i)
        {
            // Ignore external edges, they have been dealt with.
            if (i == 0 && surfHit.size() > 1)
            {
                continue;
            }

            const pointIndexHit& edHit = edHits[i];

            const label featureHit = featuresHit[i];

            if (edHit.hit())
            {
                if (!positionOnThisProc(edHit.hitPoint()))
                {
                    continue;
                }

                if (!nearFeaturePt(edHit.hitPoint()))
                {
                    if
                    (
                        magSqr(edHit.hitPoint() - surfHitI.hitPoint())
                      < surfacePtReplaceDistCoeffSqr*targetCellSizeSqr
                    )
                    {
                        // If the point is within a given distance of a feature
                        // edge, give control to edge control points instead,
                        // this will prevent "pits" forming.

                        keepSurfacePoint = false;
                    }

                    if
                    (
                        !nearFeatureEdgeLocation
                        (
                            edHit.hitPoint(),
                            newEdgeLocations,
                            existingEdgeLocations,
                            edgeLocationTree
                        )
                    )
                    {
                        // Do not place edge control points too close to a
                        // feature point or existing edge control points

                        featureEdgeHits.append(edHit);
                        featureEdgeFeaturesHit.append(featureHit);

                        currentEdgeLocations.append(edHit.hitPoint());
                    }
                }
            }
        }

        newEdgeLocations.append(currentEdgeLocations);

        if (keepSurfacePoint)
        {
            surfaceHits.append(surfHitI);

            hitSurfaces.append(hitSurfaceI);
        }
    }
}


void Foam::conformalVoronoiMesh::storeSurfaceConformation()
{
    Info<< nl << "Storing surface conformation" << endl;

    surfaceConformationVertices_.clear();

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        // Store points that are not referred, part of a pair, but not feature
        // points
        if
        (
            !vit->referred()
         && vit->pairPoint()
         && vit->index() >= startOfInternalPoints_
        )
        {
            surfaceConformationVertices_.append
            (
                Vb
                (
                    vit->point(),
                    0,                          // index, reset to zero
                    vit->type() - vit->index()  // type, relative to index
                )
            );
        }
    }

    Info<< "    Stored "
        << returnReduce
        (
            label(surfaceConformationVertices_.size()),
            sumOp<label>()
        )
        << " vertices" << endl;
}


void Foam::conformalVoronoiMesh::reinsertSurfaceConformation()
{
    Info<< nl << "Reinserting stored surface conformation" << endl;

    label preReinsertionSize(number_of_vertices());

    // It is assumed that the stored surface conformation is on the correct
    // processor and does not need distributed

    for
    (
        List<Vb>::iterator vit=surfaceConformationVertices_.begin();
        vit != surfaceConformationVertices_.end();
        ++vit
    )
    {
        // Assuming that all of the reinsertions are pair points, and that the
        // index and type are relative, i.e. index 0 and type relative to it.
        insertPoint
        (
            vit->point(),
            vit->index() + number_of_vertices(),
            vit->type() + number_of_vertices()
        );
    }

    Info<< "    Reinserted "
        << returnReduce
        (
            label(number_of_vertices()) - preReinsertionSize,
            sumOp<label>()
        )
        << " vertices" << endl;
}


// ************************************************************************* //
