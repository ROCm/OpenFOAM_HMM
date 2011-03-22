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

    startOfSurfacePointPairs_ = number_of_vertices();

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

                // Pout<< "vert " << vert << endl;
                // Pout<< "    surfHit " << surfHit << endl;

                if (surfHit.hit())
                {
                    vit->setNearBoundary();

                    if (dualCellSurfaceAnyIntersection(vit))
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
        }

        Pout<< nl << "Initial conformation" << nl
            << "    Number of vertices " << number_of_vertices() << nl
            << "    Number of surface hits " << surfaceHits.size() << nl
            << "    Number of edge hits " << featureEdgeHits.size()
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

        initialTotalHits = surfaceHits.size() + featureEdgeHits.size();
    }

    // Build the parallel interface the initial surface conformation
    buildParallelInterface();

    label iterationNo = 0;

    label maxIterations =
        cvMeshControls().maxConformationIterations(reconfMode);

    scalar iterationToIntialHitRatioLimit =
        cvMeshControls().iterationToIntialHitRatioLimit(reconfMode);

    label hitLimit = label(iterationToIntialHitRatioLimit*initialTotalHits);

    Pout<< nl << "Stopping iterations when: " << nl
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
            // points can also generate protrusions and must be assessed too.

            if (vit->nearBoundary() || vit->ppMaster())
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
            else if (vit->ppSlave())
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

        Info<< nl << "Conformation iteration " << iterationNo << nl
            << "    Number of vertices " << number_of_vertices() << nl
            << "    Number of surface hits " << surfaceHits.size() << nl
            << "    Number of edge hits " << featureEdgeHits.size()
            << endl;

        totalHits = surfaceHits.size() + featureEdgeHits.size();

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

    // Info<< nl << "After iterations, check penetrations" << endl;

    // for
    // (
    //     Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
    //     vit != finite_vertices_end();
    //     vit++
    // )
    // {
    //     if (vit->internalOrBoundaryPoint())
    //     {
    //         Foam::point vert(topoint(vit->point()));
    //         pointIndexHit surfHit;
    //         label hitSurface;

    //         dualCellLargestSurfaceProtrusion(vit, surfHit, hitSurface);

    //         if (surfHit.hit())
    //         {
    //             Info<< nl << "Residual penetration: " << nl
    //                 << vit->index() << nl
    //                 << vit->type() << nl
    //                 << vit->ppMaster() << nl
    //                 << "nearFeaturePt "
    //                 << nearFeaturePt(surfHit.hitPoint()) << nl
    //                 << vert << nl
    //                 << surfHit.hitPoint()
    //                 << endl;
    //         }
    //     }
    // }

    // {
    //     // TEST - ASSESS CLOSE SURFACE POINTS

    //     setVertexSizeAndAlignment();

    //     for
    //     (
    //         Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
    //         vit != finite_vertices_end();
    //         vit++
    //     )
    //     {
    //         if
    //         (
    //             vit->index() >= startOfSurfacePointPairs_
    //          && vit->internalOrBoundaryPoint()
    //         )
    //         {
    //             std::list<Vertex_handle> adjacentVertices;

    //             adjacent_vertices(vit, std::back_inserter(adjacentVertices));

    //             Foam::point pt = topoint(vit->point());

    //             // Info<< nl << "vit: " << vit->index() << " "
    //             //     << topoint(vit->point())
    //             //     << endl;

    //             // Info<< adjacentVertices.size() << endl;

    //             for
    //             (
    //                 std::list<Vertex_handle>::iterator
    //                 avit = adjacentVertices.begin();
    //                 avit != adjacentVertices.end();
    //                 ++avit
    //             )
    //             {
    //                 Vertex_handle avh = *avit;

    //                 // The lower indexed vertex will perform the assessment
    //                 if
    //                 (
    //                     avh->index() >= startOfSurfacePointPairs_
    //                  && avh->internalOrBoundaryPoint()
    //                  && vit->index() < avh->index()
    //                  && vit->type() != avh->type()
    //                 )
    //                 {
    //                     scalar targetSize = 0.2*averageAnyCellSize(vit, avh);

    //                     // Info<< "diff " << mag(pt - topoint(avh->point()))
    //                     //     << " " << targetSize << endl;

    //                     if
    //                     (
    //                         magSqr(pt - topoint(avh->point()))
    //                       < sqr(targetSize)
    //                     )
    //                     {
    //                         Info<< nl << "vit: " << vit->index() << " "
    //                             << topoint(vit->point())
    //                             << endl;

    //                         Info<< "    adjacent too close: "
    //                             << avh->index() << " "
    //                             << topoint(avh->point())
    //                             << endl;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    Pout<< "NOT STORING SURFACE CONFORMATION" << endl;
    // storeSurfaceConformation();
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

        // Check for the edge passing through a surface
        if (geometryToConformTo_.findSurfaceAnyIntersection(dE0, dE1))
        {
            return true;
        }
    }

    return false;
}


void Foam::conformalVoronoiMesh::buildParallelInterface()
{
    if (!Pstream::parRun())
    {
        return;
    }

    label totalInterfaceVertices = 0;

    // Propagating vertices to other processors to form halo information

    List<DynamicList<label> > verticesToProc(number_of_vertices());

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            List<label> toProcs = parallelInterfaceIntersection(vit);

            label vIndex = vit->index();

            DynamicList<label>& vertexToProc = verticesToProc[vIndex];

            forAll(toProcs, tPI)
            {
                label toProc = toProcs[tPI];

                if (toProc > -1)
                {
                    if (findIndex(vertexToProc, toProc) == -1)
                    {
                        vertexToProc.append(toProc);
                    }

                    // Refer all incident vertices to neighbour
                    // processor too
                    std::list<Vertex_handle> incidentVertices;

                    incident_vertices
                    (
                        vit,
                        std::back_inserter(incidentVertices)
                    );

                    for
                    (
                        std::list<Vertex_handle>::iterator ivit =
                            incidentVertices.begin();
                        ivit != incidentVertices.end();
                        ++ivit
                    )
                    {
                        if (!(*ivit)->farPoint())
                        {
                            label ivIndex = (*ivit)->index();

                            DynamicList<label>& iVertexToProc =
                            verticesToProc[ivIndex];

                            if (findIndex(iVertexToProc, toProc) == -1)
                            {
                                iVertexToProc.append(toProc);
                            }
                        }
                    }
                }
            }
        }
    }

    // Store the vertices that have been received and added already so that
    // there is no attempt to add them more than once.
    HashSet<labelPair, labelPair ::Hash<> > receivedVertices;

    {
        DynamicList<Foam::point> parallelInterfacePoints;
        DynamicList<label> targetProcessor;
        DynamicList<label> parallelInterfaceIndices;

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (!vit->referred())
            {
                // If a vertex has been marked to go to another processor,
                // then send it

                label vIndex = vit->index();

                const DynamicList<label>& vertexToProc = verticesToProc[vIndex];

                if (!verticesToProc.empty())
                {
                    forAll(vertexToProc, vTPI)
                    {
                        parallelInterfacePoints.append
                        (
                            topoint(vit->point())
                        );

                        targetProcessor.append(vertexToProc[vTPI]);

                        if (vit->internalOrBoundaryPoint())
                        {
                            parallelInterfaceIndices.append(vit->index());
                        }
                        else
                        {
                            parallelInterfaceIndices.append(-vit->index());
                        }
                    }
                }
            }
        }

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelInterfacePointsToSend_initial.obj",
                parallelInterfacePoints
            );
        }

        mapDistribute pointMap = buildReferringMap(targetProcessor);

        totalInterfaceVertices = parallelInterfacePoints.size();

        reduce(totalInterfaceVertices, sumOp<label>());

        pointMap.distribute(parallelInterfacePoints);

        pointMap.distribute(parallelInterfaceIndices);

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelInterfacePointsReceived_initial.obj",
                parallelInterfacePoints
            );
        }

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& constructMap = pointMap.constructMap()[domain];

            if (constructMap.size())
            {
                forAll(constructMap, i)
                {
                    label origIndex = parallelInterfaceIndices[constructMap[i]];

                    labelPair vPair(domain, origIndex);

                    // For the initial referred vertices, the original processor
                    // is the one that is sending it.
                    label encodedDomain = -(domain + 1);

                    insertPoint
                    (
                        parallelInterfacePoints[constructMap[i]],
                        origIndex,
                        encodedDomain
                    );

                    receivedVertices.insert(vPair);
                }
            }
        }
    }

    Info<< "totalInterfaceVertices " << totalInterfaceVertices << endl;

    // Re-refer vertices

    // When a real vertex has referred vertices attached to it that originate
    // on several processors, then refer all vertices to all processors, i.e.:

    //     2x--|--0a--|--1x
    //     2y--|      |--1y

    // Vertex a on processor 0 (0a) has two vertices (x and y) from processors 1
    // and 2 referred to it that are connected/incident.  We therefore say
    // that 1x and 1y and need to be referred to processor 2, and vice versa.

    // There may be multiple processors trying to supply 1x and 1y to
    // processor 2.

    // Remember which vertices were re-referred from all iterations so only
    // updates are sent.
    List<HashSet<labelPair, labelPair ::Hash<> > > allReferralVertices
    (
        Pstream::nProcs()
    );

    label totalRereferredVertices = 0;

    label nRereferIters = 0;

    // Storage for which processors to refer to
    DynamicList<label> toProcs;

    do
    {
        // For each processor (outer list), maintains a HashSet of label pairs
        // [origProc, origIndex] that are to be sent.
        List<HashSet<labelPair, labelPair ::Hash<> > > referralVertices
        (
            Pstream::nProcs()
        );

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->internalOrBoundaryPoint() || vit->pairPoint())
            {
                toProcs.clear();

                std::list<Vertex_handle> incidentVertices;

                incident_vertices
                (
                    vit,
                    std::back_inserter(incidentVertices)
                );

                // Determine all of the processors that referred vertices that
                // are attached to
                for
                (
                    std::list<Vertex_handle>::iterator ivit =
                        incidentVertices.begin();
                    ivit != incidentVertices.end();
                    ++ivit
                )
                {
                    if ((*ivit)->referred())
                    {
                        label rPI = (*ivit)->procIndex();

                        if (rPI >= 0 && findIndex(toProcs, rPI) == -1)
                        {
                            toProcs.append(rPI);
                        }
                    }
                }

                // Re-refer all referred vertices to all other connected
                // processors
                for
                (
                    std::list<Vertex_handle>::iterator ivit =
                        incidentVertices.begin();
                    ivit != incidentVertices.end();
                    ++ivit
                )
                {
                    if ((*ivit)->referred())
                    {
                        forAll(toProcs, tPI)
                        {
                            label toProc = toProcs[tPI];

                            if (toProc == (*ivit)->procIndex())
                            {
                                // There is no point referring a vertex back to
                                // its original processor
                                continue;
                            }

                            labelPair vPair
                            (
                                label((*ivit)->procIndex()),
                                label((*ivit)->index())
                            );

                            if (!allReferralVertices[toProc].found(vPair))
                            {
                                referralVertices[toProc].insert(vPair);
                                allReferralVertices[toProc].insert(vPair);
                            }
                        }
                    }
                }
            }
        }

        // Re-declare the referring structure - cannot reuse the DynamicLists as
        // they get manipulated by the mapDistribute in a way that causes
        // segmentation faults.

        DynamicList<Foam::point> parallelInterfacePoints;
        DynamicList<label> targetProcessor;
        DynamicList<label> parallelInterfaceIndices;
        DynamicList<label> originalProcessor;

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->referred())
            {
                labelPair vPair(label(vit->procIndex()), label(vit->index()));

                forAll(referralVertices, proc)
                {
                    if(referralVertices[proc].found(vPair))
                    {
                        targetProcessor.append(proc);
                        parallelInterfacePoints.append(topoint(vit->point()));
                        parallelInterfaceIndices.append(vit->index());
                        originalProcessor.append(vit->procIndex());
                    }
                }
            }
        }

        if (cvMeshControls().objOutput())
        {
            writePoints
            (
                "parallelInterfacePointsToSend_reRefer_"
                  + name(nRereferIters) + ".obj",
                parallelInterfacePoints
            );
        }

        {
            mapDistribute pointMap = buildReferringMap(targetProcessor);

            totalRereferredVertices = 0;

            pointMap.distribute(parallelInterfacePoints);
            pointMap.distribute(parallelInterfaceIndices);
            pointMap.distribute(originalProcessor);

            if (cvMeshControls().objOutput())
            {
                writePoints
                (
                    "parallelInterfacePointsReceived_reRefer_"
                      + name(nRereferIters) + ".obj",
                    parallelInterfacePoints
                );
            }

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& constructMap = pointMap.constructMap()[domain];

                if (constructMap.size())
                {
                    forAll(constructMap, i)
                    {
                        label origProc = originalProcessor[constructMap[i]];

                        label origIndex =
                            parallelInterfaceIndices[constructMap[i]];

                        labelPair vPair(origProc, origIndex);

                        if (!receivedVertices.found(vPair))
                        {
                            label encodedDomain = -(origProc + 1);

                            insertPoint
                            (
                                parallelInterfacePoints[constructMap[i]],
                                origIndex,
                                encodedDomain
                            );

                            receivedVertices.insert(vPair);

                            totalRereferredVertices++;
                        }
                    }
                }
            }

            reduce(totalRereferredVertices, sumOp<label>());
        }

        Info<< "Rerefer iteration " << nRereferIters
            << " totalRereferredVertices " << totalRereferredVertices
            << endl;

    } while(++nRereferIters < 10 && totalRereferredVertices > 0);
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


Foam::List<Foam::label>
Foam::conformalVoronoiMesh::parallelInterfaceIntersection
(
    const Delaunay::Finite_vertices_iterator& vit
) const
{
    std::list<Facet> facets;
    incident_facets(vit, std::back_inserter(facets));

    DynamicList<label> procs;

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

        Foam::point boxPt(vector::one*GREAT);
        direction ptOnFace = -1;

        Foam::point start = dE0;
        Foam::point end   = dE1;

        bool intersects = geometryToConformTo_.bounds().intersects
        (
            start,
            end - start,
            start,
            end,
            boxPt,
            ptOnFace
        );

        if (intersects && ptOnFace == 0)
        {
            // If the box is intersected, but doesn't return the appropriate
            // bits, then this means that the start point was inside the box,
            // so reverse the direction of the query.

            start = dE1;
            end   = dE0;

            geometryToConformTo_.bounds().intersects
            (
                start,
                end - start,
                start,
                end,
                boxPt,
                ptOnFace
            );
        }

        // If an intersection with the box has been found, then it is the "end"
        // point that is inside the box, so test the line from end to
        // boxPt to see if the surface is in-between.
        if
        (
            ptOnFace > 0
         && !geometryToConformTo_.findSurfaceAnyIntersection(end, boxPt)
        )
        {
            // A surface penetration is not found, and a bounds penetration is
            // found, so the dual edge has crossed the parallel interface.

            label target = targetProc(ptOnFace);

            if
            (
                target > -1
             && findIndex(procs, target) == -1
            )
            {
                procs.append(target);
            }
        }
    }

    return procs;
}


Foam::label Foam::conformalVoronoiMesh::targetProc
(
     direction ptOnFace
) const
{
    // Hard coded to 8 proc, bound-box octant split:

    if (ptOnFace == 0)
    {
        return -1;
    }

    label faceIndex = 0;

    // Taking log2 of the direction to deduce the face, i.e. reverse the
    // encoding of faceId to faceBit in treeBoundBox.H
    // From: http://graphics.stanford.edu/~seander/bithacks.html
    while (ptOnFace >>= 1)
    {
        faceIndex++;
    }

    label procArray[8][6] =
    {
        {-1,  1, -1,  2, -1,  4},
        { 0, -1, -1,  3, -1,  5},
        {-1,  3,  0, -1, -1,  6},
        { 2, -1,  1, -1, -1,  7},
        {-1,  5, -1,  6,  0, -1},
        { 4, -1, -1,  7,  1, -1},
        {-1,  7,  4, -1,  2, -1},
        { 6, -1,  5, -1,  3, -1}
    };

    return initListList<List<label>, label, 8, 6>(procArray)
        [Pstream::myProcNo()][faceIndex];
}



void Foam::conformalVoronoiMesh::dualCellLargestSurfaceProtrusion
(
    const Delaunay::Finite_vertices_iterator& vit,
    pointIndexHit& surfHitLargest,
    label& hitSurfaceLargest
) const
{
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
}


void Foam::conformalVoronoiMesh::dualCellLargestSurfaceIncursion
(
    const Delaunay::Finite_vertices_iterator& vit,
    pointIndexHit& surfHitLargest,
    label& hitSurfaceLargest
) const
{
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
                // Cannot limit displacement, point  closer than tolerance
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
        geometryToConformTo_.globalBounds().extend(rndGen_, 1E-4)
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
        geometryToConformTo_.globalBounds().extend(rndGen_, 1E-4)
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
        number_of_vertices() - startOfSurfacePointPairs_
    );

    label surfPtI = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->index() >= startOfSurfacePointPairs_)
        {
            if (!vit->pairPoint())
            {
                FatalErrorIn("storeSurfaceConformation()")
                    << "Trying to store a vertex that is not a surface point"
                    << exit(FatalError);
            }

            surfaceConformationVertices_[surfPtI] = Vb(vit->point());

            surfaceConformationVertices_[surfPtI].index() =
                vit->index() - startOfSurfacePointPairs_;

            surfaceConformationVertices_[surfPtI].type() =
                vit->type() - startOfSurfacePointPairs_;

            surfPtI++;
        }
    }

    Info<< "    Stored " << surfaceConformationVertices_.size()
        << " vertices" << endl;
}


void Foam::conformalVoronoiMesh::reinsertSurfaceConformation()
{
    Info<< nl << "Reinserting stored surface conformation" << endl;

    startOfSurfacePointPairs_ = number_of_vertices();

    forAll(surfaceConformationVertices_, v)
    {
        insertVb(surfaceConformationVertices_[v], startOfSurfacePointPairs_);
    }

    Info<< "    Reinserted " << number_of_vertices() - startOfSurfacePointPairs_
        << " vertices" << endl;
}


// ************************************************************************* //
