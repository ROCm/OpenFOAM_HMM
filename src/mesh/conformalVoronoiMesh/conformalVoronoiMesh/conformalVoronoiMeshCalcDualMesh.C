/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::calcDualMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    timeCheck();

    setVertexSizeAndAlignment();

    timeCheck();

    // ------>  OLD, FOR REFERENCE
        // // ~~~~~~~~~~~ removing short edges by indexing dual vertices ~~~~~~~~~~~~~~

        // for
        // (
        //     Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        //     cit != finite_cells_end();
        //     ++cit
        // )
        // {
        //     cit->cellIndex() = -1;
        // }

        // points.setSize(number_of_cells());

        // // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // label dualVertI = 0;

        // // Scanning by number of short (dual) edges (nSE) attached to the
        // // circumcentre of each Delaunay tet.  A Delaunay tet may only have four
        // // dual edges emanating from its circumcentre, assigning positions and
        // // indices to those with 4 short edges attached first, then >= 3, then >= 2
        // // etc.
        // for (label nSE = 4; nSE >= 0; nSE--)
        // {
        //     Info<< nl << "Scanning for dual vertices with >= "
        //         << nSE
        //         << " short edges attached." << endl;

        //     for
        //     (
        //         Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        //         cit != finite_cells_end();
        //         ++cit
        //     )
        //     {
        //         // If the Delaunay tet has an index already then it has either
        //         // evaluated itself and taken action or has had its index dictated
        //         // by a neighbouring tet with more short edges attached.

        //         if (cit->cellIndex() == -1)
        //         {
        //             point dualVertex = topoint(dual(cit));

        //             label shortEdges = 0;

        //             List<bool> edgeIsShort(4, false);

        //             List<bool> neighbourAlreadyIndexed(4, false);

        //             // Loop over the four facets of the Delaunay tet
        //             for (label f = 0; f < 4; f++)
        //             {
        //                 // Check that at least one of the vertices of the facet is
        //                 // an internal or boundary point
        //                 if
        //                 (
        //                     cit->vertex(vertex_triple_index(f, 0))->
        //                         internalOrBoundaryPoint()
        //                  || cit->vertex(vertex_triple_index(f, 1))->
        //                         internalOrBoundaryPoint()
        //                  || cit->vertex(vertex_triple_index(f, 2))->
        //                         internalOrBoundaryPoint()
        //                 )
        //                 {
        //                     point neighDualVertex;

        //                     label cNI = cit->neighbor(f)->cellIndex();

        //                     if (cNI == -1)
        //                     {
        //                         neighDualVertex = topoint(dual(cit->neighbor(f)));
        //                     }
        //                     else
        //                     {
        //                         neighDualVertex = points[cNI];
        //                     }

        //                     if
        //                     (
        //                         magSqr(dualVertex - neighDualVertex)
        //                       < sqr
        //                         (
        //                             minimumEdgeLength
        //                             (
        //                                 0.5*(dualVertex + neighDualVertex)
        //                             )
        //                         )
        //                     )
        //                     {
        //                         edgeIsShort[f] = true;

        //                         if (cNI > -1)
        //                         {
        //                             neighbourAlreadyIndexed[f] = true;
        //                         }

        //                         shortEdges++;
        //                     }
        //                 }
        //             }

        //             if (nSE == 0 && shortEdges == 0)
        //             {
        //                 // Final iteration and no short edges are found, index
        //                 // remaining dual vertices.

        //                 if
        //                 (
        //                     cit->vertex(0)->internalOrBoundaryPoint()
        //                  || cit->vertex(1)->internalOrBoundaryPoint()
        //                  || cit->vertex(2)->internalOrBoundaryPoint()
        //                  || cit->vertex(3)->internalOrBoundaryPoint()
        //                 )
        //                 {
        //                     cit->cellIndex() = dualVertI;
        //                     points[dualVertI] = dualVertex;
        //                     dualVertI++;
        //                 }
        //             }
        //             else if
        //             (
        //                 shortEdges >= nSE
        //             )
        //             {
        //                 // Info<< neighbourAlreadyIndexed << ' '
        //                 //     << edgeIsShort << endl;

        //                 label numUnindexedNeighbours = 1;

        //                 for (label f = 0; f < 4; f++)
        //                 {
        //                     if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
        //                     {
        //                         dualVertex += topoint(dual(cit->neighbor(f)));

        //                         numUnindexedNeighbours++;
        //                     }
        //                 }

        //                 dualVertex /= numUnindexedNeighbours;

        //                 label nearestExistingIndex = -1;

        //                 point nearestIndexedNeighbourPos = vector::zero;

        //                 scalar minDistSqrToNearestIndexedNeighbour = VGREAT;

        //                 for (label f = 0; f < 4; f++)
        //                 {
        //                     if (edgeIsShort[f] && neighbourAlreadyIndexed[f])
        //                     {
        //                         label cNI = cit->neighbor(f)->cellIndex();

        //                         point indexedNeighbourPos = points[cNI];

        //                         if
        //                         (
        //                             magSqr(indexedNeighbourPos - dualVertex)
        //                           < minDistSqrToNearestIndexedNeighbour
        //                         )
        //                         {
        //                             nearestExistingIndex = cNI;

        //                             nearestIndexedNeighbourPos =
        //                             indexedNeighbourPos;

        //                             minDistSqrToNearestIndexedNeighbour =
        //                             magSqr(indexedNeighbourPos - dualVertex);
        //                         }
        //                     }
        //                 }

        //                 if
        //                 (
        //                     nearestExistingIndex > -1
        //                  && minDistSqrToNearestIndexedNeighbour
        //                   < sqr
        //                     (
        //                         minimumEdgeLength
        //                         (
        //                             0.5*(nearestIndexedNeighbourPos + dualVertex)
        //                         )
        //                     )
        //                 )
        //                 {
        //                     points[nearestExistingIndex] =
        //                     0.5*(dualVertex + nearestIndexedNeighbourPos);

        //                     for (label f = 0; f < 4; f++)
        //                     {
        //                         if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
        //                         {
        //                             cit->neighbor(f)->cellIndex() =
        //                             nearestExistingIndex;
        //                         }
        //                     }

        //                     cit->cellIndex() = nearestExistingIndex;
        //                 }
        //                 else
        //                 {
        //                     for (label f = 0; f < 4; f++)
        //                     {
        //                         if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
        //                         {
        //                             cit->neighbor(f)->cellIndex() = dualVertI;
        //                         }
        //                     }

        //                     cit->cellIndex() = dualVertI;

        //                     points[dualVertI] = dualVertex;

        //                     dualVertI++;
        //                 }
        //             }
        //         }
        //     }
        // }

        // points.setSize(dualVertI);

    // <------ OLD, FOR REFERENCE

    // Dual cell indexing

    // Assign an index to the Delaunay vertices which will be the dual cell
    // index used for owner neighbour assignment.

    // The indices of the points are reset *which **destroys** the point-pair
    // matching*, so the type of each vertex are reset to avoid any ambiguity.

    label dualCelli = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            vit->type() = Vb::ptInternalPoint;
            vit->index() = dualCelli;
            dualCelli++;
        }
        else
        {
            vit->type() = Vb::ptFarPoint;
            vit->index() = -1;
        }
    }

    // Dual face filtering

    // Indexing Delaunay cells, which are the dual vertices

    label dualVertI = 0;

    points.setSize(number_of_cells());

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
        {
            cit->cellIndex() = dualVertI;
            points[dualVertI] = topoint(dual(cit));
            dualVertI++;
        }
        else
        {
            cit->cellIndex() = -1;
        }
    }

    points.setSize(dualVertI);

    // Merge close points

    Info<< nl << "    Merging close points" << endl;

    mergeCloseDualVertices(points);

    // Smooth the surface of the mesh

    Info<< nl << "    Smoothing surface" << endl;

    smoothSurface(points);

    // Collapse faces throughout the mesh

    Info<< nl << "    Collapsing unnecessary faces" << endl;

    collapseFaces(points);

    // Final dual face and owner neighbour construction

    timeCheck();

    createFacesOwnerNeighbourAndPatches
    (
        points,
        faces,
        owner,
        neighbour,
        patchNames,
        patchSizes,
        patchStarts,
        false
    );
}


void Foam::conformalVoronoiMesh::calcTetMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    labelList vertexMap(number_of_vertices());

    label vertI = 0;

    points.setSize(number_of_vertices());

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() || vit->pairPoint())
        {
            vertexMap[vit->index()] = vertI;
            points[vertI] =  topoint(vit->point());
            vertI++;
        }
    }

    points.setSize(vertI);

    label cellI = 0;

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
        {
             cit->cellIndex() = cellI++;
        }
        else
        {
            cit->cellIndex() = -1;
        }
    }

    patchNames = geometryToConformTo_.patchNames();

    patchNames.setSize(patchNames.size() + 1);

    patchNames[patchNames.size() - 1] = "cvMesh_defaultPatch";

    label nPatches = patchNames.size();

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_facets());

    owner.setSize(number_of_facets());

    neighbour.setSize(number_of_facets());

    label faceI = 0;

    labelList verticesOnTriFace(3, -1);

    face newFace(verticesOnTriFace);

    for
    (
        Triangulation::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const int oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        label c1I = c1->cellIndex();
        label c2I = c2->cellIndex();

        label ownerCell = -1;
        label neighbourCell = -1;

        if (c1I == -1 && c2I == -1)
        {
            // Both tets are outside, skip
            continue;
        }

        for (label i = 0; i < 3; i++)
        {
            verticesOnTriFace[i] = vertexMap
            [
                c1->vertex(vertex_triple_index(oppositeVertex, i))->index()
            ];
        }

        newFace = face(verticesOnTriFace);

        if (c1I == -1 || c2I == -1)
        {
            // Boundary face...
            if (c1I == -1)
            {
                //... with c1 outside
                ownerCell = c2I;
            }
            else
            {
                // ... with c2 outside
                ownerCell = c1I;

                reverse(newFace);
            }

            label patchIndex = geometryToConformTo_.findPatch
            (
                newFace.centre(points)
            );

            if (patchIndex == -1)
            {
                patchIndex = patchNames.size() - 1;

                WarningIn("Foam::conformalVoronoiMesh::calcTetMesh")
                    << "Tet face centre at  " << nl
                    << "    " << newFace.centre(points) << nl
                    << "    did not find a surface patch. Adding to "
                    << patchNames[patchIndex]
                    << endl;
            }

            patchFaces[patchIndex].append(newFace);
            patchOwners[patchIndex].append(ownerCell);
        }
        else
        {
            // Internal face...
            if (c1I < c2I)
            {
                // ...with c1 as the ownerCell
                ownerCell = c1I;
                neighbourCell = c2I;

                reverse(newFace);
            }
            else
            {
                // ...with c2 as the ownerCell
                ownerCell = c2I;
                neighbourCell = c1I;
            }

            faces[faceI] = newFace;
            owner[faceI] = ownerCell;
            neighbour[faceI] = neighbourCell;
            faceI++;
        }
    }

    label nInternalFaces = faceI;

    faces.setSize(nInternalFaces);
    owner.setSize(nInternalFaces);
    neighbour.setSize(nInternalFaces);

    sortFaces(faces, owner, neighbour);

    addPatches
    (
        nInternalFaces,
        faces,
        owner,
        patchNames,
        patchSizes,
        patchStarts,
        patchFaces,
        patchOwners,
        false
    );
}


void Foam::conformalVoronoiMesh::mergeCloseDualVertices(const pointField& pts)
{
    // Assess close points to be merged

    label nPtsMerged = 0;

    do
    {
        Map<label> dualPtIndexMap;

        nPtsMerged = mergeCloseDualVertices(pts, dualPtIndexMap);

        if (nPtsMerged > 0)
        {
            Info<< "        Merged " << nPtsMerged << " points "
                << "closenessTolerance HARDCODED " << endl;
        }

        reindexDualVertices(dualPtIndexMap);

    } while (nPtsMerged > 0);
}


Foam::label Foam::conformalVoronoiMesh::mergeCloseDualVertices
(
    const pointField& pts,
    Map<label>& dualPtIndexMap
) const
{
    label nPtsMerged = 0;

    scalar closenessTolerance = 1e-6;

    for
    (
        Triangulation::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const int oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        label& c1I = c1->cellIndex();
        label& c2I = c2->cellIndex();

        if (dualPtIndexMap.found(c1I) || dualPtIndexMap.found(c2I))
        {
            // One of the points of this edge has already been
            // merged this sweep, leave for next sweep

            continue;
        }

        if(c1I != -1 && c2I != -1 && (c1I != c2I))
        {
            if
            (
                magSqr(pts[c1I] - pts[c2I])
              < sqr(averageAnyCellSize(fit)*closenessTolerance)
            )
            {
                dualPtIndexMap.insert(c1I, c1I);
                dualPtIndexMap.insert(c2I, c1I);

                nPtsMerged++;
            }
        }
    }

    return nPtsMerged;
}


void Foam::conformalVoronoiMesh::smoothSurface(pointField& pts)
{
    label nCollapsedFaces = 0;

    do
    {
        Map<label> dualPtIndexMap;

        nCollapsedFaces = smoothSurfaceDualFaces(pts, dualPtIndexMap);

        if (nCollapsedFaces > 0)
        {
            Info<< "        Collapsed " << nCollapsedFaces << " boundary faces"
                << endl;
        }

        reindexDualVertices(dualPtIndexMap);

        mergeCloseDualVertices(pts);

    } while (nCollapsedFaces > 0);

    // Force all points of the face to be on the surface

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        label ptI = cit->cellIndex();

        // Only cells with indices > -1 are valid
        if (ptI > -1)
        {
            // Test if this is a boundary dual vertex - if it is, at
            // least one of the Delaunay vertices of the Delaunay cell
            // will be outside

            if
            (
                !cit->vertex(0)->internalOrBoundaryPoint()
             || !cit->vertex(1)->internalOrBoundaryPoint()
             || !cit->vertex(2)->internalOrBoundaryPoint()
             || !cit->vertex(3)->internalOrBoundaryPoint()
            )
            {
                point& pt = pts[ptI];

                pointIndexHit surfHit;
                label hitSurface;

                geometryToConformTo_.findSurfaceNearest
                (
                    pt,
                    cvMeshControls().spanSqr(),
                    surfHit,
                    hitSurface
                );

                if (surfHit.hit())
                {
                    pt = surfHit.hitPoint();
                }
            }
        }

    }

    mergeCloseDualVertices(pts);
}


Foam::label Foam::conformalVoronoiMesh::smoothSurfaceDualFaces
(
    pointField& pts,
    Map<label>& dualPtIndexMap
) const
{
    label nCollapsedFaces = 0;

    const scalar cosPerpendicularToleranceAngle = cos(degToRad(75));

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_circulator ccStart = incident_cells(*eit);
        Cell_circulator cc = ccStart;

        do
        {
            if (dualPtIndexMap.found(cc->cellIndex()))
            {
                // One of the points of this face has already been
                // collapsed this sweep, leave for next sweep
                continue;
            }

        } while (++cc != ccStart);

        if (isBoundaryDualFace(eit))
        {
            face dualFace = buildDualFace(eit);

            if (dualFace.size() < 3)
            {
                // This face has been collapsed already
                continue;
            }

            pointIndexHit surfHit;
            label hitSurface;

            geometryToConformTo_.findSurfaceNearest
            (
                dualFace.centre(pts),
                cvMeshControls().spanSqr(),
                surfHit,
                hitSurface
            );

            vectorField norm(1);

            allGeometry_[hitSurface].getNormal
            (
                List<pointIndexHit>(1, surfHit),
                norm
            );

            const vector& surfaceNormal = norm[0];

            // Orient the face correctly before calculating the normal

            Cell_handle c = eit->first;
            Vertex_handle vA = c->vertex(eit->second);
            Vertex_handle vB = c->vertex(eit->third);

            if (!vA->internalOrBoundaryPoint())
            {
                reverse(dualFace);
            }

            vector faceNormal = dualFace.normal(pts);

            faceNormal /= mag(faceNormal);

            if ((faceNormal & surfaceNormal) < cosPerpendicularToleranceAngle)
            {
                scalar targetFaceSize = averageAnyCellSize(vA, vB);

                // Selecting faces to collapse based on angle to
                // surface, so set collapseSizeLimitCoeff to GREAT to
                // allow collapse of all faces

                if
                (
                    collapseFaceToEdge
                    (
                        dualFace,
                        pts,
                        dualPtIndexMap,
                        targetFaceSize,
                        GREAT
                    )
                )
                {
                    nCollapsedFaces++;
                }
            }
        }
    }

    return nCollapsedFaces;
}


void Foam::conformalVoronoiMesh::collapseFaces(pointField& pts)
{
    label nCollapsedFaces = 0;

    do
    {
        Map<label> dualPtIndexMap;

        nCollapsedFaces = collapseFaces(pts, dualPtIndexMap);

        reindexDualVertices(dualPtIndexMap);

        mergeCloseDualVertices(pts);

        if (nCollapsedFaces > 0)
        {
            Info<< "        Collapsed " << nCollapsedFaces << " faces" << endl;
        }

    } while (nCollapsedFaces > 0);
}


Foam::label Foam::conformalVoronoiMesh::collapseFaces
(
    pointField& pts,
    Map<label>& dualPtIndexMap
) const
{
    label nCollapsedFaces = 0;

    scalar collapseSizeLimitCoeff = cvMeshControls().minimumEdgeLengthCoeff();

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_circulator ccStart = incident_cells(*eit);
        Cell_circulator cc = ccStart;

        do
        {
            if (dualPtIndexMap.found(cc->cellIndex()))
            {
                // One of the points of this face has already been
                // collapsed this sweep, leave for next sweep
                continue;
            }

        } while (++cc != ccStart);

        Cell_handle c = eit->first;
        Vertex_handle vA = c->vertex(eit->second);
        Vertex_handle vB = c->vertex(eit->third);

        if
        (
            vA->internalOrBoundaryPoint()
         || vB->internalOrBoundaryPoint()
        )
        {
            face dualFace = buildDualFace(eit);

            if (dualFace.size() < 3)
            {
                // This face has been collapsed already
                continue;
            }

            scalar targetFaceSize = averageAnyCellSize(vA, vB);

            if
            (
                collapseFaceToEdge
                (
                    dualFace,
                    pts,
                    dualPtIndexMap,
                    targetFaceSize,
                    collapseSizeLimitCoeff
                )
            )
            {
                nCollapsedFaces++;
            }
        }
    }

    // ------>  OLD, FOR REFERENCE

    // scalar smallEdgeLengthCoeff = 1e-3;
    // scalar smallFaceAreaCoeff = sqr(smallEdgeLengthCoeff);
    // scalar collapseToEdgeCoeff = 0.02;
    // scalar longestEdgeLengthRatio = 0.35;

    // for
    // (
    //     Triangulation::Finite_edges_iterator eit = finite_edges_begin();
    //     eit != finite_edges_end();
    //     ++eit
    // )
    // {
    //     Cell_circulator ccStart = incident_cells(*eit);
    //     Cell_circulator cc = ccStart;

    //     do
    //     {
    //         if (dualPtIndexMap.found(cc->cellIndex()))
    //         {
    //             // One of the points of this face has already been
    //             // collapsed this sweep, leave for next sweep
    //             continue;
    //         }

    //     } while (++cc != ccStart);

    //     Cell_handle c = eit->first;
    //     Vertex_handle vA = c->vertex(eit->second);
    //     Vertex_handle vB = c->vertex(eit->third);

    //     if
    //     (
    //         vA->internalOrBoundaryPoint()
    //      || vB->internalOrBoundaryPoint()
    //     )
    //     {
    //         face dualFace = buildDualFace(eit);

    //         if (dualFace.size() < 3)
    //         {
    //             // This face has been collapsed already
    //             continue;
    //         }

    //         scalar area = dualFace.mag(pts);

    //         scalar targetFaceSize = averageAnyCellSize(vA, vB);
    //         scalar targetArea = sqr(targetFaceSize);

    //         if (area < smallFaceAreaCoeff*targetArea)
    //         {
    //             // Collapse the dual face

    //             // Determine if the face should be collapsed to a line or a
    //             // point

    //             const edgeList& eds = dualFace.edges();

    //             label longestEdgeI = -1;

    //             scalar longestEdgeLength = -SMALL;

    //             scalar perimeter = 0.0;

    //             forAll(eds, edI)
    //             {
    //                 scalar edgeLength = eds[edI].mag(pts);

    //                 perimeter += edgeLength;

    //                 if (edgeLength > longestEdgeLength)
    //                 {
    //                     longestEdgeI = edI;

    //                     longestEdgeLength = edgeLength;
    //                 }
    //             }

    //             if
    //             (
    //                 longestEdgeLength > collapseToEdgeCoeff*targetFaceSize
    //              && longestEdgeLength/perimeter > longestEdgeLengthRatio
    //             )
    //             {
    //                 // Collapse to edge

    //                 // Start at either end of the longest edge and consume the
    //                 // rest of the points of the face

    //                 const edge& longestEd = eds[longestEdgeI];

    //                 label longestEdStartPtI = longestEd.start();
    //                 label longestEdEndPtI = longestEd.end();

    //                 label revEdI = longestEdgeI;
    //                 label fwdEdI = longestEdgeI;

    //                 point revPt = pts[longestEdStartPtI];
    //                 point fwdPt = pts[longestEdEndPtI];

    //                 dualPtIndexMap.insert(longestEdStartPtI, longestEdStartPtI);
    //                 dualPtIndexMap.insert(longestEdEndPtI, longestEdEndPtI);

    //                 for (label fcI = 1; fcI <= label(eds.size()/2); fcI++)
    //                 {
    //                     revEdI = eds.rcIndex(revEdI);
    //                     fwdEdI = eds.fcIndex(fwdEdI);

    //                     const edge& revEd = eds[revEdI];
    //                     const edge& fwdEd = eds[fwdEdI];

    //                     if (fcI < label(eds.size()/2))
    //                     {
    //                         revPt += pts[revEd.start()];
    //                         fwdPt += pts[fwdEd.end()];

    //                         dualPtIndexMap.insert
    //                         (
    //                             revEd.start(),
    //                             longestEdStartPtI
    //                         );

    //                         dualPtIndexMap.insert
    //                         (
    //                             fwdEd.end(),
    //                             longestEdEndPtI
    //                         );
    //                     }
    //                     else
    //                     {
    //                         // Final circulation

    //                         if
    //                         (
    //                             eds.size() % 2 == 1
    //                             && revEd.start() == fwdEd.end()
    //                         )
    //                         {
    //                             // Odd number of edges, give final point to
    //                             // the edge direction that has the shorter
    //                             // final edge

    //                             if (fwdEd.mag(pts) < revEd.mag(pts))
    //                             {
    //                                 fwdPt += pts[fwdEd.end()];

    //                                 dualPtIndexMap.insert
    //                                 (
    //                                     fwdEd.end(),
    //                                     longestEdEndPtI
    //                                 );

    //                                 revPt /= fcI;
    //                                 fwdPt /= (fcI + 1);
    //                             }
    //                             else
    //                             {
    //                                 revPt += pts[revEd.start()];

    //                                 dualPtIndexMap.insert
    //                                 (
    //                                     revEd.start(),
    //                                     longestEdStartPtI
    //                                 );

    //                                 revPt /= (fcI + 1);
    //                                 fwdPt /= fcI;
    //                             }
    //                         }
    //                         else if
    //                         (
    //                             eds.size() % 2 == 0
    //                          && revEd.start() == fwdEd.start()
    //                          && revEd.end() == fwdEd.end()
    //                         )
    //                         {
    //                             // Even number of edges

    //                             revPt /= fcI;
    //                             fwdPt /= fcI;
    //                         }
    //                         else
    //                         {
    //                             FatalErrorIn
    //                             (
    //                                 "Foam::conformalVoronoiMesh::collapseFace"
    //                             )
    //                                 << "Face circulation failed for face "
    //                                 << dualFace << nl
    //                                 << exit(FatalError);
    //                         }
    //                     }
    //                 }

    //                 // Move the position of the accumulated points
    //                 pts[longestEdStartPtI] = revPt;
    //                 pts[longestEdEndPtI] = fwdPt;

    //                 nCollapsedFaces++;
    //             }
    //             else if
    //             (
    //                 longestEdgeLength <= collapseToEdgeCoeff*targetFaceSize
    //             )
    //             {
    //                 // Collapse to point

    //                 point resultantPt = vector::zero;

    //                 label collapseToPtI = dualFace[0];

    //                 forAll(dualFace, fPtI)
    //                 {
    //                     label ptI = dualFace[fPtI];

    //                     resultantPt += pts[ptI];

    //                     dualPtIndexMap.insert(ptI, collapseToPtI);
    //                 }

    //                 resultantPt /= dualFace.size();

    //                 pts[collapseToPtI] = resultantPt;

    //                 nCollapsedFaces++;
    //             }
    //         }
    //     }
    // }

    // <------ OLD, FOR REFERENCE

    return nCollapsedFaces;
}


bool Foam::conformalVoronoiMesh::collapseFaceToEdge
(
    const face& f,
    pointField& pts,
    Map<label>& dualPtIndexMap,
    scalar targetFaceSize,
    scalar collapseSizeLimitCoeff
) const
{
    scalar guardFraction = 0.5;

    const vector fC = f.centre(pts);

    vector fN = f.normal(pts);

    const scalar fA = mag(fN);

    tensor J = f.inertia(pts, fC);

    // Find the dominant collapse direction by finding the eigenvector
    // that corresponds to the normal direction, discarding it.  The
    // eigenvector corresponding to the smaller of the two remaining
    // eigenvalues is the dominant axis in a high aspect ratio face.

    // Normalise inertia tensor to remove problems with small values

    scalar magJ = mag(J);

    if (magJ > VSMALL)
    {
        J /= mag(J);
        // J /= cmptMax(J);
        // J /= max(eigenValues(J).x(), SMALL);
    }
    else
    {
        WarningIn
        (
            "Foam::conformalVoronoiMesh::collapseFaceToEdge"
            "("
            "const face& f,"
            "pointField& pts,"
            "Map<label>& dualPtIndexMap"
            ") const"
        )
            << "Inertia tensor magnitude too small, not collapsing." << nl
            << J << nl << "mag = " << magJ
            << endl;
    }

    // ***************************************************************
    // The maximum eigenvalue (z()) must be the direction of the
    // normal, as it has the greatest value.  The minimum eigenvalue
    // is the dominant collapse axis for high aspect ratio faces.
    // These assumption can be tested with the code below; they are
    // only likely to be invalid for very warped faces.

    // // Face normal is now a unit vector
    // fN /= fA;

    // direction normalMatchEigenValueCmpt = -1;

    // scalar normalMatchMaxMagDotProd = -SMALL;

    // for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    // {
    //     vector eVec = eigenVector(J, eVals.component(cmpt));

    //     scalar magDotProd = mag(eVec & fN);

    //     if (magDotProd > normalMatchMaxMagDotProd)
    //     {
    //         normalMatchMaxMagDotProd = magDotProd;

    //         normalMatchEigenValueCmpt = cmpt;
    //     }
    // }

    // direction faceAxisACmpt = (normalMatchEigenValueCmpt + 1) % 3;
    // vector faceAxisA = eigenVector(J, eVals.component(faceAxisACmpt));

    // direction faceAxisBCmpt = (normalMatchEigenValueCmpt + 2) % 3;
    // vector faceAxisB = eigenVector(J, eVals.component(faceAxisBCmpt));

    // Info<< nl << "# Normal match " << normalMatchEigenValueCmpt << nl
    //     << "# " << eVals.component(faceAxisACmpt) << nl
    //     << "# " << eVals.component(faceAxisBCmpt) << nl
    //     << "# Aspect ratio = " << sqrt
    //     (
    //         eVals.component(faceAxisBCmpt)
    //        /eVals.component(faceAxisACmpt)
    //     )
    //     << endl;

    // scalar scale = 2.0*mag(fC - pts[f[0]])/eVals.component(faceAxisACmpt);

    // meshTools::writeOBJ(Info, fC);
    // meshTools::writeOBJ
    // (
    //     Info,
    //     fC + scale*eVals.component(faceAxisACmpt)*faceAxisA
    // );
    // meshTools::writeOBJ
    // (
    //     Info,
    //     fC + scale*eVals.component(faceAxisBCmpt)*faceAxisB
    // );

    // Info<< "f 1 2" << nl << "f 1 3" << endl;

    // forAll(f, fPtI)
    // {
    //     meshTools::writeOBJ(Info, pts[f[fPtI]]);
    // }

    // Info<< "f";

    // forAll(f, fPtI)
    // {
    //     Info << " " << fPtI + 4;
    // }

    // Info << endl;

    // ***************************************************************

    vector collapseAxis = vector::zero;

    scalar aspectRatio = 1;

    scalar detJ = det(J);

    if (detJ < 1e-5)
    {
        const edgeList& eds = f.edges();

        label longestEdgeI = -1;

        scalar longestEdgeLength = -SMALL;

        forAll(eds, edI)
        {
            scalar edgeLength = eds[edI].mag(pts);

            if (edgeLength > longestEdgeLength)
            {
                longestEdgeI = edI;

                longestEdgeLength = edgeLength;
            }
        }

        collapseAxis = eds[longestEdgeI].vec(pts);

        if (mag(collapseAxis) < VSMALL)
        {
            Info<< "if (mag(collapseAxis) < VSMALL) " << collapseAxis << endl;
        }

        collapseAxis /= mag(collapseAxis);

        // Empirical correlation for high aspect ratio faces

        if (detJ < VSMALL)
        {
            Info<< "if (detJ < VSMALL) " << detJ << endl;
        }

        aspectRatio = sqrt(0.35/detJ);

        // Info<< "# Longest edge determined collapseAxis" << endl;
    }
    else
    {
        vector eVals = eigenValues(J);

        // The inertia calculation describes the mass distribution as a
        // function of distance squared to the axis, so the square root of
        // the ratio of face-plane moments gives a good indication of the
        // aspect ratio.

        aspectRatio = sqrt(eVals.y()/max(eVals.x(), SMALL));

        collapseAxis = eigenVector(J, eVals.x());
    }

    if (magSqr(collapseAxis) < VSMALL)
    {
        WarningIn
        (
            "Foam::conformalVoronoiMesh::collapseFaceToEdge"
            "("
                "const face& f,"
                "pointField& pts,"
                "Map<label>& dualPtIndexMap"
            ") const"
        )
            << "No collapse axis found for face, not collapsing."
            << endl;

        // // Output face and collapse axis for visualisation

        // Info<< nl << "# Aspect ratio = " << aspectRatio << nl
        //     << "# collapseAxis = " << collapseAxis << nl
        //     << "# eigenvalues = " << eVals << endl;

        // scalar scale = 2.0*mag(fC - pts[f[0]]);

        // meshTools::writeOBJ(Info, fC);
        // meshTools::writeOBJ(Info, fC + scale*collapseAxis);

        // Info<< "f 1 2" << endl;

        // forAll(f, fPtI)
        // {
        //     meshTools::writeOBJ(Info, pts[f[fPtI]]);
        // }

        // Info<< "f";

        // forAll(f, fPtI)
        // {
        //     Info << " " << fPtI + 3;
        // }

        // Info<< nl << endl;

        return false;
    }

    // The signed distance along the collapse axis passing through the
    // face centre that each vertex projects to.

    Field<scalar> d(f.size());

    forAll(f, fPtI)
    {
        const point& pt = pts[f[fPtI]];

        d[fPtI] = (collapseAxis & (pt - fC));
    }

    // Sort the projected distances and the corresponding vertex
    // indices along the collapse axis

    labelList facePts(f);

    labelList oldToNew;

    sortedOrder(d, oldToNew);

    oldToNew = invert(oldToNew.size(), oldToNew);

    inplaceReorder(oldToNew, d);

    inplaceReorder(oldToNew, facePts);

    // Shift the points so that they are relative to the centre of the
    // collapse line.

    scalar dShift = -0.5*(d.first() + d.last());

    d += dShift;

//     // Output face and collapse axis for visualisation

//     Info<< "# Aspect ratio = " << aspectRatio << nl
//         << "# determinant = " << detJ << nl
//         << "# collapseAxis = " << collapseAxis << nl
// //        << "# eigenvalues = " << eVals
//         << endl;

//     scalar scale = 2.0*mag(fC - pts[f[0]]);

//     meshTools::writeOBJ(Info, fC);
//     meshTools::writeOBJ(Info, fC + scale*collapseAxis);

//     Info<< "f 1 2" << endl;

//     forAll(f, fPtI)
//     {
//         meshTools::writeOBJ(Info, pts[f[fPtI]]);
//     }

//     Info<< "f";

//     forAll(f, fPtI)
//     {
//         Info << " " << fPtI + 3;
//     }

//     Info<< nl << "# " << d << endl;

//     Info<< "# " << d.first() << " " << d.last() << endl;

//     forAll(d, dI)
//     {
//         meshTools::writeOBJ(Info, fC + (d[dI] - dShift)*collapseAxis);
//     }

//     Info<< endl;

    // Form two lists, one for each half of the set of points
    // projected along the collapse axis.

    // Middle value, index of first entry in the second half
    label middle = -1;

    forAll(d, dI)
    {
        if (d[dI] > 0)
        {
            middle = dI;

            break;
        }
    }

    // Negative Half
    SubList<scalar> dNeg(d, middle, 0);
    SubList<label> facePtsNeg(facePts, middle, 0);

    // Positive Half
    SubList<scalar> dPos(d, d.size() - middle, middle);
    SubList<label> facePtsPos(facePts, d.size() - middle, middle);

    // Defining how close to the midpoint (M) of the projected
    // vertices line a projected vertex (X) can be before making this
    // an invalid edge collapse
    //
    // X---X-g----------------M----X-----------g----X--X
    //
    // Only allow a collapse if all projected vertices are outwith
    // guardFraction (g) of the distance form the face centre to the
    // furthest vertex in the considered direction

    if (dNeg.size() == 0 || dPos.size() == 0)
    {
        WarningIn
        (
            "Foam::conformalVoronoiMesh::collapseFaceToEdge"
            "("
                "const face& f,"
                "pointField& pts,"
                "Map<label>& dualPtIndexMap"
            ") const"
        )
            << "All points on one side of face centre, not collapsing."
            << endl;
    }

    bool validCollapse = false;

    if
    (
        (dNeg.last() < guardFraction*dNeg.first())
     && (dPos.first() > guardFraction*dPos.last())
     && (fA < aspectRatio*sqr(targetFaceSize*collapseSizeLimitCoeff))
     && f.size() <= 4
    )
    {
        validCollapse = true;
    }

    if (validCollapse)
    {
        // Arbitrarily choosing the most distant point as the index to
        // collapse to.

        label collapseToPtI = facePtsNeg.first();

        forAll(facePtsNeg, fPtI)
        {
            dualPtIndexMap.insert(facePtsNeg[fPtI], collapseToPtI);
        }

        pts[collapseToPtI] = collapseAxis*(sum(dNeg)/dNeg.size() - dShift) + fC;

        collapseToPtI = facePtsPos.last();

        forAll(facePtsPos, fPtI)
        {
            dualPtIndexMap.insert(facePtsPos[fPtI], collapseToPtI);
        }

        pts[collapseToPtI] = collapseAxis*(sum(dPos)/dPos.size() - dShift) + fC;
    }
    else
    {
        // If the face can't be collapsed to a line, and it is small
        // and low aspect ratio enough, collapse it to a point.

        if
        (
            (d.last() - d.first()) < targetFaceSize*0.35
         && fA < aspectRatio*sqr(targetFaceSize*collapseSizeLimitCoeff)
         && f.size() <= 4
        )
        {
            // Arbitrarily choosing the first point as the index to
            // collapse to.  Collapse to the face center.

            label collapseToPtI = facePts.first();

            forAll(facePts, fPtI)
            {
                dualPtIndexMap.insert(facePts[fPtI], collapseToPtI);
            }

            pts[collapseToPtI] = fC;

            validCollapse = true;
        }

        // Alternatively, do not topologically collapse face, but push
        // all points onto a line, so that the face area is zero and
        // either:
        //   + do not create it when dualising.  This may damage the edge
        //     addressing of the mesh;
        //   + split the face into two (or more?) edges later, sacrificing
        //     topological consistency with the Delaunay.
    }

    return validCollapse;
}


void Foam::conformalVoronoiMesh::reindexDualVertices
(
    const Map<label>& dualPtIndexMap
)
{
    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (dualPtIndexMap.found(cit->cellIndex()))
        {
            cit->cellIndex() = dualPtIndexMap[cit->cellIndex()];
        }
    }
}


void Foam::conformalVoronoiMesh::createFacesOwnerNeighbourAndPatches
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    bool includeEmptyPatches
) const
{
    patchNames = geometryToConformTo_.patchNames();

    patchNames.setSize(patchNames.size() + 1);

    patchNames[patchNames.size() - 1] = "cvMesh_defaultPatch";

    label nPatches = patchNames.size();

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_edges());

    owner.setSize(number_of_edges());

    neighbour.setSize(number_of_edges());

    label dualFaceI = 0;

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_handle c = eit->first;
        Vertex_handle vA = c->vertex(eit->second);
        Vertex_handle vB = c->vertex(eit->third);

        if
        (
            vA->internalOrBoundaryPoint()
         || vB->internalOrBoundaryPoint()
        )
        {
            face newDualFace = buildDualFace(eit);

            if (newDualFace.size() >= 3)
            {
                label dcA = vA->index();

                if (!vA->internalOrBoundaryPoint())
                {
                    dcA = -1;
                }

                label dcB = vB->index();

                if (!vB->internalOrBoundaryPoint())
                {
                    dcB = -1;
                }

                label dcOwn = -1;
                label dcNei = -1;

                if (dcA == -1 && dcB == -1)
                {
                    FatalErrorIn("calcDualMesh")
                        << "Attempting to create a face joining "
                        << "two unindexed dual cells "
                        << exit(FatalError);
                }
                else if (dcA == -1 || dcB == -1)
                {
                    // boundary face, find which is the owner

                    if (dcA == -1)
                    {
                        dcOwn = dcB;

                        // reverse face order to correctly orientate normal
                        reverse(newDualFace);
                    }
                    else
                    {
                        dcOwn = dcA;
                    }

                    // Find which patch this face is on by finding the
                    // intersection with the surface of the Delaunay edge
                    // generating the face and identify the region of the
                    // intersection.

                    point ptA = topoint(vA->point());

                    point ptB = topoint(vB->point());

                    label patchIndex = geometryToConformTo_.findPatch(ptA, ptB);

                    if (patchIndex == -1)
                    {
                        patchIndex = patchNames.size() - 1;

                        WarningIn("Foam::conformalVoronoiMesh::calcDualMesh")
                            << "Dual face found between Dv pair " << nl
                            << "    " << ptA << nl
                            << "    " << ptB << nl
                            << "    that is not on a surface patch. Adding to "
                            << patchNames[patchIndex]
                            << endl;
                    }

                    patchFaces[patchIndex].append(newDualFace);
                    patchOwners[patchIndex].append(dcOwn);
                }
                else
                {
                    // internal face, find the lower cell to be the owner

                    if (dcB > dcA)
                    {
                        dcOwn = dcA;
                        dcNei = dcB;
                    }
                    else
                    {
                        dcOwn = dcB;
                        dcNei = dcA;

                        // reverse face order to correctly orientate normal
                        reverse(newDualFace);
                    }

                    faces[dualFaceI] = newDualFace;

                    owner[dualFaceI] = dcOwn;

                    neighbour[dualFaceI] = dcNei;

                    dualFaceI++;
                }
            }
        }
    }

    label nInternalFaces = dualFaceI;

    faces.setSize(nInternalFaces);
    owner.setSize(nInternalFaces);
    neighbour.setSize(nInternalFaces);

    timeCheck();

    sortFaces(faces, owner, neighbour);

    timeCheck();

    addPatches
    (
        nInternalFaces,
        faces,
        owner,
        patchNames,
        patchSizes,
        patchStarts,
        patchFaces,
        patchOwners,
        false
    );

    removeUnusedPoints(faces, points);

    timeCheck();
}


void Foam::conformalVoronoiMesh::sortFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour
) const
{
    // Upper triangular order:
    // + owner is sorted in ascending cell order
    // + within each block of equal value for owner, neighbour is sorted in
    //   ascending cell order.
    // + faces sorted to correspond
    // e.g.
    // owner | neighbour
    // 0     | 2
    // 0     | 23
    // 0     | 71
    // 1     | 23
    // 1     | 24
    // 1     | 91

    // Two stage sort:
    // 1) sort by owner

    Info<< nl
        << "    Sorting faces, owner and neighbour into upper triangular order"
        << endl;

    labelList oldToNew;

    sortedOrder(owner, oldToNew);

    oldToNew = invert(oldToNew.size(), oldToNew);

    inplaceReorder(oldToNew, faces);
    inplaceReorder(oldToNew, owner);
    inplaceReorder(oldToNew, neighbour);

    // 2) in each block of owners sort by neighbour

    // Reset map.  Elements that are not sorted will retain their -1
    // value, which will mean that they are ignored by inplaceReorder

    oldToNew = -1;

    label ownerBlockStart = 0;

    for (label o = 1; o < owner.size(); o++)
    {
        label blockLength = -1;

        if (owner[o] > owner[o-1])
        {
            blockLength = o - ownerBlockStart;
        }
        else if (o == owner.size() - 1)
        {
            // If the last element is not a jump in owner, then it
            // needs to trigger a sort of the last block, but with a
            // block length that is one element longer so that it
            // sorts itself.

            // If it is a jump in owner, then it will form a block of
            // length one, and so will not need sorted.

            blockLength = o - ownerBlockStart + 1;
        }

        if (blockLength >= 1)
        {
            labelList blockIndices =
                identity(blockLength) + ownerBlockStart;

            SubList<label> neighbourBlock
            (
                neighbour,
                blockLength,
                ownerBlockStart
            );

            sortedOrder(neighbourBlock, blockIndices);

            blockIndices = invert(blockIndices.size(), blockIndices);

            forAll(blockIndices, b)
            {
                oldToNew[ownerBlockStart + b] =
                blockIndices[b] + ownerBlockStart;
            }

            ownerBlockStart = o;
        }
    }

    // owner does not need re-sorted
    inplaceReorder(oldToNew, faces);
    inplaceReorder(oldToNew, neighbour);
}


void Foam::conformalVoronoiMesh::addPatches
(
    const label nInternalFaces,
    faceList& faces,
    labelList& owner,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    List<DynamicList<face> >& patchFaces,
    List<DynamicList<label> >& patchOwners,
    bool includeEmptyPatches
) const
{
    label nTotalPatches = patchNames.size();

    label nValidPatches = 0;

    PackedBoolList validPatch(nTotalPatches, false);

    wordList allPatchNames = patchNames;

    patchSizes.setSize(nTotalPatches);
    patchStarts.setSize(nTotalPatches);

    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        // Check if the patch has any faces.  Never create an empty
        // default patch.

        if
        (
            patchFaces[p].size()
         || (includeEmptyPatches && (p != nTotalPatches - 1))
        )
        {
            patchNames[nValidPatches] = allPatchNames[p];
            patchSizes[nValidPatches] = patchFaces[p].size();
            patchStarts[nValidPatches] = nInternalFaces + nBoundaryFaces;

            nBoundaryFaces += patchSizes[p];

            nValidPatches++;

            validPatch[p] = 1;
        }
        else
        {
            // Warn if a patch is empty and includeEmptyPatches is
            // false, unless it is the default patch.

            if (p != nTotalPatches - 1)
            {
                WarningIn("void addPatches")
                    << "Patch " << patchNames[p]
                    << " has no faces, not creating." << endl;
            }
        }
    }

    patchNames.setSize(nValidPatches);
    patchSizes.setSize(nValidPatches);
    patchStarts.setSize(nValidPatches);

    faces.setSize(nInternalFaces + nBoundaryFaces);
    owner.setSize(nInternalFaces + nBoundaryFaces);

    label faceI = nInternalFaces;

    forAll(patchFaces, p)
    {
        if (validPatch[p])
        {
            forAll(patchFaces[p], f)
            {
                faces[faceI] = patchFaces[p][f];
                owner[faceI] = patchOwners[p][f];

                faceI++;
            }
        }
    }
}


void Foam::conformalVoronoiMesh::removeUnusedPoints
(
    faceList& faces,
    pointField& pts
) const
{
    Info<< nl << "    Removing unused points after filtering" << endl;

    PackedBoolList ptUsed(pts.size(), false);

    // Scan all faces to find all of the points that are used

    forAll(faces, fI)
    {
        const face& f = faces[fI];

        forAll(f, fPtI)
        {
            ptUsed[f[fPtI]] = true;
        }
    }

    label pointI = 0;
    labelList oldToNew(pts.size(), -1);

    // Move all of the used faces to the start of the pointField and
    // truncate it

    forAll(ptUsed, ptUI)
    {
        if (ptUsed[ptUI] == true)
        {
            oldToNew[ptUI] = pointI++;
        }
    }

    inplaceReorder(oldToNew, pts);

    Info<< "        Removing "
        << pts.size() - pointI
        << " unused points" << endl;

    pts.setSize(pointI);

    // Renumber the faces to use the new point numbers

    forAll(faces, fI)
    {
        inplaceRenumber(oldToNew, faces[fI]);
    }
}


// ************************************************************************* //
