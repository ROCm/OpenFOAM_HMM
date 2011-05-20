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

void Foam::conformalVoronoiMesh::conformToSurface()
{
    reconformationMode reconfMode = reconformationControl();

    if (reconfMode == rmNone)
    {
        // Reinsert stored surface conformation
        reinsertSurfaceConformation();
    }
    else
    {
        // Rebuild, insert and store new surface conformation
        buildSurfaceConformation(reconfMode);
    }
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

    startOfSurfacePoints_ = number_of_vertices();

    // Initialise containers to store the edge conformation locations
    DynamicList<Foam::point> newEdgeLocations;

    pointField existingEdgeLocations(0);

    autoPtr<indexedOctree<treeDataPoint> > edgeLocationTree;

    // Initialise the edgeLocationTree
    buildEdgeLocationTree(edgeLocationTree, existingEdgeLocations);

    label initialTotalHits = 0;

    // Initial surface protrusion conformation - nearest surface point
    {
        scalar edgeSearchDistCoeffSqr =
            cvMeshControls().edgeSearchDistCoeffSqrInitial(reconfMode);

        scalar surfacePtReplaceDistCoeffSqr =
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
                Foam::point vert(topoint(vit->point()));
                scalar searchDistanceSqr = surfaceSearchDistanceSqr(vert);
                pointIndexHit surfHit;
                label hitSurface;

                geometryToConformTo_.findSurfaceNearest
                (
                    vert,
                    searchDistanceSqr,
                    surfHit,
                    hitSurface
                );

                if (surfHit.hit())
                {
                    vit->setNearBoundary();

                    if (dualCellSurfaceAnyIntersection(vit))
                    {
                        // meshTools::writeOBJ(Pout, vert);
                        // meshTools::writeOBJ(Pout, surfHit.hitPoint());
                        // Pout<< "l cr0 cr1" << endl;

                        addSurfaceAndEdgeHits
                        (
                            vit,
                            vert,
                            surfHit,
                            hitSurface,
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

    scalar iterationToIntialHitRatioLimit =
        cvMeshControls().iterationToIntialHitRatioLimit(reconfMode);

    label hitLimit = label(iterationToIntialHitRatioLimit*initialTotalHits);

    Info<< nl << "Stopping iterations when: " << nl
        <<"    total number of hits drops below "
        << iterationToIntialHitRatioLimit << " of initial hits ("
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
             || vit->referredInternal()
            )
            {
                Foam::point vert(topoint(vit->point()));
                pointIndexHit surfHit;
                label hitSurface;

                dualCellLargestSurfaceProtrusion(vit, surfHit, hitSurface);

                if (surfHit.hit())
                {
                    addSurfaceAndEdgeHits
                    (
                        vit,
                        vert,
                        surfHit,
                        hitSurface,
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
                Foam::point vert(topoint(vit->point()));
                pointIndexHit surfHit;
                label hitSurface;

                dualCellLargestSurfaceIncursion(vit, surfHit, hitSurface);

                if (surfHit.hit())
                {
                    addSurfaceAndEdgeHits
                    (
                        vit,
                        vert,
                        surfHit,
                        hitSurface,
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

    // reportSurfaceConformationQuality();

    storeSurfaceConformation();
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

        // if (Pstream::parRun())
        // {
        //     const treeBoundBoxList& procBbs =
        //         geometryToConformTo_.processorDomains()[Pstream::myProcNo()];

        //     forAll(procBbs, pBI)
        //     {
        //         const treeBoundBox& procBb = procBbs[pBI];

        //         Foam::point a = dE0;
        //         Foam::point b = dE1;

        //         bool inBox = clipLineToBox(a, b, procBb);

        //         // Check for the edge passing through a surface
        //         if
        //         (
        //             inBox
        //          && geometryToConformTo_.findSurfaceAnyIntersection(a, b)
        //         )
        //         {
        //             // Pout<< "# findSurfaceAnyIntersection" << endl;
        //             // meshTools::writeOBJ(Pout, a);
        //             // meshTools::writeOBJ(Pout, b);
        //             // Pout<< "l cr0 cr1" << endl;

        //             return true;
        //         }
        //     }
        // }
        // else
        // {
            if (geometryToConformTo_.findSurfaceAnyIntersection(dE0, dE1))
            {
                return true;
            }
        // }
    }

    return false;
}


bool Foam::conformalVoronoiMesh::clipLineToBox
(
    Foam::point& a,
    Foam::point& b,
    const treeBoundBox& box
) const
{
    Foam::point boxPt = vector::one*GREAT;

    bool intersects = false;

    if (box.posBits(a) == 0)
    {
        // a is inside
        if (box.posBits(b) == 0)
        {
            // both a and b inside.
            intersects = true;
        }
        else
        {
            // b is outside, clip b to bounding box.
            intersects = box.intersects(b, a, boxPt);
            b = boxPt;
        }
    }
    else
    {
        // a is outside
        if (box.posBits(b) == 0)
        {
            // b is inside
            intersects = box.intersects(a, b, boxPt);
            a = boxPt;
        }
        else
        {
            // both a and b outside, but they can still intersect the box

            intersects = box.intersects(a, b, boxPt);

            if (intersects)
            {
                a = boxPt;
                box.intersects(b, a, boxPt);
                b = boxPt;
            }
        }
    }

    return intersects;
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

    // {
    //     // Update the processorMeshBounds

    //     DynamicList<Foam::point> parallelAllPoints;
    //     DynamicList<label> targetProcessor;
    //     DynamicList<label> parallelAllIndices;

    //     Foam::point minPt = geometryToConformTo_.bounds().min();
    //     Foam::point maxPt = geometryToConformTo_.bounds().max();

    //     for
    //     (
    //         Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
    //         vit != finite_vertices_end();
    //         vit++
    //     )
    //     {
    //         if (vit->real())
    //         {
    //             Foam::point v = topoint(vit->point());

    //             minPt = Foam::min(minPt, v);
    //             maxPt = Foam::max(maxPt, v);
    //         }
    //     }

    //     treeBoundBox& procMeshBb =
    //         geometryToConformTo_.processorMeshBounds()[Pstream::myProcNo()];

    //     if (cvMeshControls().objOutput())
    //     {
    //         Pout<< "Before processorMeshBounds update" << procMeshBb << endl;
    //     }

    //     procMeshBb = treeBoundBox(minPt, maxPt);

    //     if (cvMeshControls().objOutput())
    //     {
    //         Pout<< "After processorMeshBounds update" << procMeshBb << endl;

    //         OFstream str
    //         (
    //             runTime_.path()
    //             /"processorMeshBoundsUpdated_"
    //           + name(Pstream::myProcNo())
    //           + "_bounds.obj"
    //         );

    //         Pout<< "Writing " << str.name() << endl;

    //         pointField bbPoints(procMeshBb.points());

    //         forAll(bbPoints, i)
    //         {
    //             meshTools::writeOBJ(str, bbPoints[i]);
    //         }

    //         forAll(treeBoundBox::faces, i)
    //         {
    //             const face& f = treeBoundBox::faces[i];

    //             str << "f"
    //                 << ' ' << f[0] + 1
    //                 << ' ' << f[1] + 1
    //                 << ' ' << f[2] + 1
    //                 << ' ' << f[3] + 1
    //                 << nl;
    //         }
    //     }

    //     Pstream::gatherList(geometryToConformTo_.processorMeshBounds());

    //     Pstream::scatterList(geometryToConformTo_.processorMeshBounds());
    // }

    boolList sendToProc(Pstream::nProcs(), false);

    // Hard coded switch, can be turned on for debugging purposes and all
    // vertices will be referred to all processors.
    bool allPointReferral = true;

    if (allPointReferral)
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
                sendToProc = true;

                sendToProc[Pstream::myProcNo()] = false;

                forAll(sendToProc, procI)
                {
                    if (sendToProc[procI])
                    {
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
        }

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelAllPointsToSend.obj",
                parallelAllPoints
            );
        }

        mapDistribute pointMap = buildReferringMap(targetProcessor);

        label totalAllVertices = parallelAllPoints.size();

        reduce(totalAllVertices, sumOp<label>());

        pointMap.distribute(parallelAllPoints);

        pointMap.distribute(parallelAllIndices);

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelAllPointsReceived_" + outputName + ".obj",
                parallelAllPoints
            );
        }

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            const labelList& constructMap = pointMap.constructMap()[procI];

            if (constructMap.size())
            {
                forAll(constructMap, i)
                {
                    label origIndex =
                        parallelAllIndices[constructMap[i]];

                    if (!receivedVertices[procI].found(origIndex))
                    {
                        // For the initial referred vertices, the original
                        // processor is the one that is sending it.
                        label encodedProcI = -(procI + 1);

                        insertPoint
                        (
                            parallelAllPoints[constructMap[i]],
                            origIndex,
                            encodedProcI
                        );

                        receivedVertices[procI].insert(origIndex);
                    }
                }
            }
        }

        Info<< "totalAllVertices "
            << totalAllVertices << endl;
    }

    if (initialEdgeReferral)
    {
        DynamicList<Foam::point> parallelIntersectionPoints;
        DynamicList<label> targetProcessor;
        DynamicList<label> parallelIntersectionIndices;

        // Initial pass - find vertices whose dual edges pierce nearby
        // processor volumes and refer them to establish a sensible boundary
        // interface region before running a circumsphere assessment.
        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->internalOrBoundaryPoint())
            {
                parallelInterfaceIntersection(vit, sendToProc);

                forAll(sendToProc, procI)
                {
                    if (sendToProc[procI])
                    {
                        label vIndex = vit->index();

                        // Using the hashSet to ensure that each vertex is only
                        // referred once to each processor
                        if (!referralVertices[procI].found(vIndex))
                        {
                            referralVertices[procI].insert(vIndex);

                            parallelIntersectionPoints.append
                            (
                                topoint(vit->point())
                            );

                            targetProcessor.append(procI);

                            if (vit->internalOrBoundaryPoint())
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
                            //     << " " << procI
                            //     << endl;
                        }
                    }
                }
            }
        }

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelIntersectionPointsToSend_" + outputName + ".obj",
                parallelIntersectionPoints
            );
        }

        mapDistribute pointMap = buildReferringMap(targetProcessor);

        label totalIntersectionVertices = parallelIntersectionPoints.size();

        reduce(totalIntersectionVertices, sumOp<label>());

        pointMap.distribute(parallelIntersectionPoints);

        pointMap.distribute(parallelIntersectionIndices);

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelIntersectionPointsReceived_" + outputName + ".obj",
                parallelIntersectionPoints
            );
        }

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            const labelList& constructMap = pointMap.constructMap()[procI];

            if (constructMap.size())
            {
                forAll(constructMap, i)
                {
                    label origIndex =
                        parallelIntersectionIndices[constructMap[i]];

                    if (!receivedVertices[procI].found(origIndex))
                    {
                        // For the initial referred vertices, the original
                        // processor is the one that is sending it.
                        label encodedProcI = -(procI + 1);

                        // Pout<< "Insert "
                        //     << parallelIntersectionPoints[constructMap[i]]
                        //     << " " << origIndex
                        //     << " " << procI
                        //     << endl;

                        insertPoint
                        (
                            parallelIntersectionPoints[constructMap[i]],
                            origIndex,
                            encodedProcI
                        );

                        receivedVertices[procI].insert(origIndex);
                    }
                }
            }
        }

        Info<< "totalIntersectionVertices "
            << totalIntersectionVertices << endl;
    }

    {
        DynamicList<Foam::point> parallelInfluencePoints;
        DynamicList<label> targetProcessor;
        DynamicList<label> parallelInfluenceIndices;

        for
        (
            Delaunay::Finite_cells_iterator cit = finite_cells_begin();
            cit != finite_cells_end();
            ++cit
        )
        {
            parallelInterfaceInfluence(cit, sendToProc);

            forAll(sendToProc, procI)
            {
                if (sendToProc[procI])
                {
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
            }
        }

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelInfluencePointsToSend_" + outputName + ".obj",
                parallelInfluencePoints
            );
        }

        mapDistribute pointMap = buildReferringMap(targetProcessor);

        label totalInfluenceVertices = parallelInfluencePoints.size();

        reduce(totalInfluenceVertices, sumOp<label>());

        pointMap.distribute(parallelInfluencePoints);

        pointMap.distribute(parallelInfluenceIndices);

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelInfluencePointsReceived_" + outputName + ".obj",
                parallelInfluencePoints
            );
        }

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            const labelList& constructMap = pointMap.constructMap()[procI];

            if (constructMap.size())
            {
                forAll(constructMap, i)
                {
                    label origIndex = parallelInfluenceIndices[constructMap[i]];

                    if (!receivedVertices[procI].found(origIndex))
                    {
                        // For the initial referred vertices, the original
                        // processor is the one that is sending it.
                        label encodedProcI = -(procI + 1);

                        // Pout<< "Insert "
                        //     << parallelInfluencePoints[constructMap[i]]
                        //     << " " << origIndex
                        //     << " " << procI
                        //     << endl;

                        insertPoint
                        (
                            parallelInfluencePoints[constructMap[i]],
                            origIndex,
                            encodedProcI
                        );

                        receivedVertices[procI].insert(origIndex);
                    }
                }
            }
        }

        Info<< "totalInfluenceVertices " << totalInfluenceVertices << endl;
    }
}


Foam::mapDistribute Foam::conformalVoronoiMesh::buildReferringMap
(
    const DynamicList<label>& targetProcessor
) const
{
    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // 1. Count
    labelList nSend(Pstream::nProcs(), 0);

    forAll(targetProcessor, i)
    {
        label procI = targetProcessor[i];

        nSend[procI]++;
    }

    // Send over how many I need to receive
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList sendSizes(Pstream::nProcs());

    sendSizes[Pstream::myProcNo()] = nSend;

    combineReduce(sendSizes, UPstream::listEq());

    // 2. Size sendMap
    labelListList sendMap(Pstream::nProcs());

    forAll(nSend, procI)
    {
        sendMap[procI].setSize(nSend[procI]);

        nSend[procI] = 0;
    }

    // 3. Fill sendMap
    forAll(targetProcessor, i)
    {
        label procI = targetProcessor[i];

        sendMap[procI][nSend[procI]++] = i;
    }

    // Determine receive map
    // ~~~~~~~~~~~~~~~~~~~~~

    labelListList constructMap(Pstream::nProcs());

    // Local transfers first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label constructSize = constructMap[Pstream::myProcNo()].size();

    forAll(constructMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            label nRecv = sendSizes[procI][Pstream::myProcNo()];

            constructMap[procI].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[procI][i] = constructSize++;
            }
        }
    }

    return mapDistribute
    (
        constructSize,
        sendMap.xfer(),
        constructMap.xfer()
    );
}

void Foam::conformalVoronoiMesh::parallelInterfaceIntersection
(
    const Delaunay::Finite_vertices_iterator& vit,
    boolList& toProc
) const
{
    toProc = false;

    Foam::point vert(topoint(vit->point()));

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

        Foam::point boxPt = vector::one*GREAT;
        scalar hitDistSqr = GREAT;
        bool intersects = false;

        label closestHitProc = -1;
        scalar closestHitDistSqr = GREAT;

        // forAll(geometryToConformTo_.processorDomains(), procI)
        // {
        //     if (procI == Pstream::myProcNo())
        //     {
        //         continue;
        //     }

        //     const treeBoundBoxList& procBbs =
        //         geometryToConformTo_.processorDomains()[procI];

        //     forAll(procBbs, pBI)
        //     {
        //         const treeBoundBox& procBb = procBbs[pBI];

        //         intersects = procBb.intersects(dE0, dE1, boxPt);

        //         if (intersects)
        //         {
        //             hitDistSqr = magSqr(vert - boxPt);

        //             if (hitDistSqr < closestHitDistSqr)
        //             {
        //                 closestHitProc = procI;
        //                 closestHitDistSqr = hitDistSqr;
        //             }
        //         }

        //         // Perform the query in the opposite direction
        //         intersects = procBb.intersects(dE1, dE0, boxPt);

        //         if (intersects)
        //         {
        //             hitDistSqr = magSqr(vert - boxPt);

        //             if (hitDistSqr < closestHitDistSqr)
        //             {
        //                 closestHitProc = procI;
        //                 closestHitDistSqr = hitDistSqr;
        //             }
        //         }
        //     }
        // }

        if (closestHitProc >= 0)
        {
            toProc[closestHitProc] = true;
        }
    }
}


void Foam::conformalVoronoiMesh::parallelInterfaceInfluence
(
    const Delaunay::Finite_cells_iterator& cit,
    boolList& toProc
) const
{
    // Assess the influence of the circumsphere of each Delaunay cell with the
    // defining volumes for all processors.  Any processor touched by the
    // circumsphere requires all points of the cell to be referred to it.

    toProc = false;

    // The Delaunay cells to assess have to be real, i.e. all vertices form
    // part of the internal or boundary definition
    if
    (
        cit->vertex(0)->internalOrBoundaryPoint()
     || cit->vertex(1)->internalOrBoundaryPoint()
     || cit->vertex(2)->internalOrBoundaryPoint()
     || cit->vertex(3)->internalOrBoundaryPoint()
    )
    {
        Foam::point circumcentre = topoint(dual(cit));

        scalar circumradiusSqr = magSqr
        (
            circumcentre - topoint(cit->vertex(0)->point())
        );

        // Pout<< nl << "# circumradius " << sqrt(circumradiusSqr) << endl;
        // drawDelaunayCell(Pout, cit);

        // forAll(geometryToConformTo_.processorMeshBounds(), procI)
        // {
        //     if (procI == Pstream::myProcNo())
        //     {
        //         continue;
        //     }

        //     const treeBoundBox& procBb =
        //         geometryToConformTo_.processorMeshBounds()[procI];

        //     if (procBb.overlaps(circumcentre, circumradiusSqr))
        //     {
        //         toProc[procI] = true;
        //     }
        // }
    }
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

    Foam::point vert(topoint(vit->point()));

    scalar maxProtrusionDistance = maxSurfaceProtrusion(vert);

    for
    (
        std::list<Facet>::iterator fit=facets.begin();
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

                scalar normalProtrusionDistance =
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
        std::list<Facet>::iterator fit=facets.begin();
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
            if
            (
                vit->index() >= startOfSurfacePoints_
             && vit->internalOrBoundaryPoint()
            )
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
                        avh->index() >= startOfSurfacePoints_
                     && avh->internalOrBoundaryPoint()
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
    const pointIndexHit& surfHit,
    label hitSurface,
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

    if (nearFeaturePt(surfHit.hitPoint()))
    {
        keepSurfacePoint = false;

        if (vit->index() < startOfInternalPoints_)
        {
            surfaceHits.append(surfHit);

            hitSurfaces.append(hitSurface);
        }
    }

    List<pointIndexHit> edHits;

    labelList featuresHit;

    scalar targetCellSizeSqr = sqr(targetCellSize(vert));

    geometryToConformTo_.findEdgeNearestByType
    (
        surfHit.hitPoint(),
        edgeSearchDistCoeffSqr*targetCellSizeSqr,
        edHits,
        featuresHit
    );

    // Gather edge locations but do not add them to newEdgeLocations inside the
    // loop as they will prevent nearby edge locations of different types being
    // conformed to.

    DynamicList<Foam::point> currentEdgeLocations;

    forAll(edHits, i)
    {
        const pointIndexHit& edHit(edHits[i]);

        label featureHit = featuresHit[i];

        if (edHit.hit())
        {
            if(!positionOnThisProc(edHit.hitPoint()))
            {
                continue;
            }

            if (!nearFeaturePt(edHit.hitPoint()))
            {
                if
                (
                    magSqr(edHit.hitPoint() - surfHit.hitPoint())
                  < surfacePtReplaceDistCoeffSqr*targetCellSizeSqr
                )
                {
                    // If the point is within a given distance of a feature
                    // edge, give control to edge control points instead, this
                    // will prevent "pits" forming.

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
                    // Do not place edge control points too close to a feature
                    // point or existing edge control points

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
        surfaceHits.append(surfHit);

        hitSurfaces.append(hitSurface);
    }
}


void Foam::conformalVoronoiMesh::storeSurfaceConformation()
{
    Info<< nl << "Storing surface conformation" << endl;

    surfaceConformationVertices_.setSize
    (
        number_of_vertices() - startOfSurfacePoints_
    );

    label surfPtI = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (!vit->referred() && vit->index() >= startOfSurfacePoints_)
        {
            if (!vit->pairPoint())
            {
                FatalErrorIn("storeSurfaceConformation()")
                    << "Trying to store a vertex that is not a surface point"
                    << exit(FatalError);
            }

            surfaceConformationVertices_[surfPtI] = Vb(vit->point());

            surfaceConformationVertices_[surfPtI].index() =
                vit->index() - startOfSurfacePoints_;

            surfaceConformationVertices_[surfPtI].type() =
                vit->type() - startOfSurfacePoints_;

            surfPtI++;
        }
    }

    Info<< "    Stored " << surfaceConformationVertices_.size()
        << " vertices" << endl;
}


void Foam::conformalVoronoiMesh::reinsertSurfaceConformation()
{
    Info<< nl << "Reinserting stored surface conformation" << endl;

    startOfSurfacePoints_ = number_of_vertices();

    forAll(surfaceConformationVertices_, v)
    {
        insertVb(surfaceConformationVertices_[v], startOfSurfacePoints_);
    }

    Info<< "    Reinserted " << number_of_vertices() - startOfSurfacePoints_
        << " vertices" << endl;
}


// ************************************************************************* //
