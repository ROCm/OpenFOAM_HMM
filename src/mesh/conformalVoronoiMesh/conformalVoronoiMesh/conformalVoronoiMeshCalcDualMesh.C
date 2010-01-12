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
#include "motionSmoother.H"

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

    // Dual cell indexing

    // Assign an index to the Delaunay vertices which will be the dual cell
    // index used for owner neighbour assignment.

    // The indices of the points are reset *which **destroys** the point-pair
    // matching*, so the type of each vertex is reset to avoid any ambiguity.

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

    HashSet<labelPair, labelPair::Hash<> > deferredCollapseFaces;

    collapseFaces(points, deferredCollapseFaces);

    label nBadQualityFaces = checkPolyMeshQuality(points);

    Info<< "Found " << nBadQualityFaces << " bad quality faces" << endl;

    // Final dual face and owner neighbour construction

    timeCheck();

    createFacesOwnerNeighbourAndPatches
    (
        faces,
        owner,
        neighbour,
        patchNames,
        patchSizes,
        patchStarts,
        false
    );

    deferredCollapseFaceSet(owner, neighbour, deferredCollapseFaces);

    removeUnusedCells(owner, neighbour);

    removeUnusedPoints(faces, points);

    timeCheck();
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

    scalar closenessTolerance = cvMeshControls().mergeClosenessCoeff();

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

    const scalar cosPerpendicularToleranceAngle = cos
    (
        degToRad(cvMeshControls().surfaceStepFaceAngle())
    );

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_circulator ccStart = incident_cells(*eit);
        Cell_circulator cc = ccStart;

        bool skipFace = false;

        do
        {
            if (dualPtIndexMap.found(cc->cellIndex()))
            {
                // One of the points of this face has already been
                // collapsed this sweep, leave for next sweep

                skipFace = true;

                break;
            }

        } while (++cc != ccStart);

        if (skipFace)
        {
            continue;
        }

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

            if (mag(faceNormal) < VSMALL)
            {
                // If the face is essentially zero area, then force it
                // to be collapsed by making the dot product result -1
                faceNormal = -surfaceNormal;
            }
            else
            {
                faceNormal /= mag(faceNormal);
            }

            if ((faceNormal & surfaceNormal) < cosPerpendicularToleranceAngle)
            {
                scalar targetFaceSize = averageAnyCellSize(vA, vB);

                // Selecting faces to collapse based on angle to
                // surface, so set collapseSizeLimitCoeff to GREAT to
                // allow collapse of all faces

                faceCollapseMode mode = collapseFace
                (
                    dualFace,
                    pts,
                    dualPtIndexMap,
                    targetFaceSize,
                    GREAT
                );

                if (mode == fcmPoint || mode == fcmEdge)
                {
                    nCollapsedFaces++;
                }
            }
        }
    }

    return nCollapsedFaces;
}


void Foam::conformalVoronoiMesh::collapseFaces
(
    pointField& pts,
    HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
)
{
    label nCollapsedFaces = 0;

    do
    {
        Map<label> dualPtIndexMap;

        deferredCollapseFaces.clear();

        nCollapsedFaces = collapseFaces
        (
            pts,
            dualPtIndexMap,
            deferredCollapseFaces
        );

        reindexDualVertices(dualPtIndexMap);

        mergeCloseDualVertices(pts);

        if (nCollapsedFaces > 0)
        {
            Info<< "        Collapsed " << nCollapsedFaces << " faces" << endl;
            // Info<< "dualPtIndexMap" << nl << dualPtIndexMap << endl;
        }

    } while (nCollapsedFaces > 0);
}


Foam::label Foam::conformalVoronoiMesh::collapseFaces
(
    pointField& pts,
    Map<label>& dualPtIndexMap,
    HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
) const
{
    label nCollapsedFaces = 0;

    scalar collapseSizeLimitCoeff = cvMeshControls().filterSizeCoeff();

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_circulator ccStart = incident_cells(*eit);
        Cell_circulator cc = ccStart;

        bool skipFace = false;

        do
        {
            if (dualPtIndexMap.found(cc->cellIndex()))
            {
                // One of the points of this face has already been
                // collapsed this sweep, leave for next sweep

                skipFace = true;

                break;
            }

        } while (++cc != ccStart);

        if (skipFace)
        {
            continue;
        }

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

            faceCollapseMode mode = collapseFace
            (
                dualFace,
                pts,
                dualPtIndexMap,
                targetFaceSize,
                collapseSizeLimitCoeff
            );

            if (mode != fcmNone)
            {
                if (mode == fcmDeferredMultiEdge)
                {
                    // Determine the owner and neighbour labels

                    Pair<label> ownAndNei(-1, -1);

                    ownerAndNeighbour
                    (
                        vA,
                        vB,
                        ownAndNei.first(),
                        ownAndNei.second()
                    );

                    // Record the owner and neighbour of this face for a
                    // deferredMultiEdge collapse

                    deferredCollapseFaces.insert(ownAndNei);
                }
                else
                {
                    nCollapsedFaces++;
                }
            }
        }
    }

    return nCollapsedFaces;
}


Foam::conformalVoronoiMesh::faceCollapseMode
Foam::conformalVoronoiMesh::collapseFace
(
    const face& f,
    pointField& pts,
    Map<label>& dualPtIndexMap,
    scalar targetFaceSize,
    scalar collapseSizeLimitCoeff
) const
{
    bool limitToQuadsOrTris = true;

    bool allowEarlyCollapseToPoint = true;

    const vector fC = f.centre(pts);

    vector fN = f.normal(pts);

    const scalar fA = mag(fN);

    tensor J = f.inertia(pts, fC);

    // Find the dominant collapse direction by finding the eigenvector
    // that corresponds to the normal direction, discarding it.  The
    // eigenvector corresponding to the smaller of the two remaining
    // eigenvalues is the dominant axis in a high aspect ratio face.

    scalar magJ = mag(J);

    scalar detJ = SMALL;

    if (magJ > VSMALL)
    {
        // Normalise inertia tensor to remove problems with small values

        J /= mag(J);
        // J /= cmptMax(J);
        // J /= max(eigenValues(J).x(), SMALL);

        // Calculating determinant, including stabilisation for zero or
        // small negative values

        detJ = max(det(J), SMALL);
    }

    vector collapseAxis = vector::zero;

    scalar aspectRatio = 1;

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

        aspectRatio = sqrt(0.35/detJ);

        // Info<< "# Longest edge determined collapseAxis" << endl;
    }
    else
    {
        vector eVals = eigenValues(J);

        // The maximum eigenvalue (z()) must be the direction of the
        // normal, as it has the greatest value.  The minimum eigenvalue
        // is the dominant collapse axis for high aspect ratio faces.

        collapseAxis = eigenVector(J, eVals.x());

        // The inertia calculation describes the mass distribution as a
        // function of distance squared to the axis, so the square root of
        // the ratio of face-plane moments gives a good indication of the
        // aspect ratio.

        aspectRatio = sqrt(eVals.y()/max(eVals.x(), SMALL));
    }

    if (magSqr(collapseAxis) < VSMALL)
    {
        WarningIn
        (
            "Foam::conformalVoronoiMesh::collapseFace"
            "("
                "const face& f,"
                "pointField& pts,"
                "Map<label>& dualPtIndexMap"
            ") const"
        )
            << "No collapse axis found for face, not collapsing."
            << endl;

        return fcmNone;
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
            "Foam::conformalVoronoiMesh::collapseFace"
            "("
                "const face& f,"
                "pointField& pts,"
                "Map<label>& dualPtIndexMap"
            ") const"
        )
            << "All points on one side of face centre, not collapsing."
            << endl;
    }

    faceCollapseMode mode = fcmNone;

    if
    (
        (fA < aspectRatio*sqr(targetFaceSize*collapseSizeLimitCoeff))
     && (!limitToQuadsOrTris || f.size() <= 4)
    )
    {
        scalar guardFraction = cvMeshControls().edgeCollapseGuardFraction();

        cvMeshControls().maxCollapseFaceToPointSideLengthCoeff();

        if
        (
            allowEarlyCollapseToPoint
         && (d.last() - d.first())
          < targetFaceSize
           *0.2*cvMeshControls().maxCollapseFaceToPointSideLengthCoeff()
        )
        {
            mode = fcmPoint;
        }
        else if
        (
            (dNeg.last() < guardFraction*dNeg.first())
         && (dPos.first() > guardFraction*dPos.last())
        )
        {
            mode = fcmEdge;
        }
        else if
        (
            (d.last() - d.first())
          < targetFaceSize
           *cvMeshControls().maxCollapseFaceToPointSideLengthCoeff()
        )
        {
            // If the face can't be collapsed to an edge, and it has a
            // small enough span, collapse it to a point.

            mode = fcmPoint;
        }
        else
        {
            // Alternatively, do not topologically collapse face here,
            // but push all points onto a line, so that the face area
            // is zero and then collapse to a string of edges later.
            // The fcmDeferredMultiEdge collapse must be performed at
            // the polyMesh stage as this type of collapse can't be
            // performed and still maintain topological dual
            // consistency with the Delaunay structure

            mode = fcmDeferredMultiEdge;
        }
    }

    switch (mode)
    {
        case fcmEdge:
        {
            // Arbitrarily choosing the most distant point as the index to
            // collapse to.

            label collapseToPtI = facePtsNeg.first();

            forAll(facePtsNeg, fPtI)
            {
                dualPtIndexMap.insert(facePtsNeg[fPtI], collapseToPtI);
            }

            pts[collapseToPtI] =
                 collapseAxis*(sum(dNeg)/dNeg.size() - dShift) + fC;

            collapseToPtI = facePtsPos.last();

            forAll(facePtsPos, fPtI)
            {
                dualPtIndexMap.insert(facePtsPos[fPtI], collapseToPtI);
            }

            pts[collapseToPtI] =
                collapseAxis*(sum(dPos)/dPos.size() - dShift) + fC;

            break;
        }

        case fcmPoint:
        {
            // Arbitrarily choosing the first point as the index to
            // collapse to.  Collapse to the face centre.

            label collapseToPtI = facePts.first();

            forAll(facePts, fPtI)
            {
                dualPtIndexMap.insert(facePts[fPtI], collapseToPtI);
            }

            pts[collapseToPtI] = fC;

            break;
        }

        case fcmDeferredMultiEdge:
        {
            // forAll(facePts, fPtI)
            // {
            //     label ptI = facePts[fPtI];

            //     pts[ptI] = collapseAxis*(d[fPtI] - dShift) + fC;

            //     dualPtIndexMap.insert(ptI, ptI);
            // }

            break;
        }

        case fcmNone:
        {
            break;
        }
    }

    // if (mode == fcmDeferredMultiEdge)
    // if (mode != fcmNone)
    // {
    //     // Output face and collapse axis for visualisation

    //     Info<< "# Aspect ratio = " << aspectRatio << nl
    //         << "# determinant = " << detJ << nl
    //         << "# collapseAxis = " << collapseAxis << nl
    //         << "# mode = " << mode << nl
    //         << "# facePts = " << facePts << nl
    //     //        << "# eigenvalues = " << eVals
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
    // }

    return mode;
}


void Foam::conformalVoronoiMesh::deferredCollapseFaceSet
(
    labelList& owner,
    labelList& neighbour,
    const HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
) const
{
    DynamicList<label> faceLabels;

    forAll (neighbour, nI)
    {
        if (deferredCollapseFaces.found(Pair<label>(owner[nI], neighbour[nI])))
        {
            faceLabels.append(nI);
        }
    }

    Info<< "facesToCollapse" << nl << faceLabels << endl;
}


Foam::label Foam::conformalVoronoiMesh::checkPolyMeshQuality
(
    const pointField& pts
) const
{
    faceList faces;
    labelList owner;
    labelList neighbour;
    wordList patchNames;
    labelList patchSizes;
    labelList patchStarts;

    timeCheck();

    Info<< "    Creating polyMesh to assess quality" << endl;

    createFacesOwnerNeighbourAndPatches
    (
        faces,
        owner,
        neighbour,
        patchNames,
        patchSizes,
        patchStarts,
        false
    );

    removeUnusedCells(owner, neighbour);

    IOobject io
    (
        "cvMesh_temporary",
        runTime_.timeName(),
        runTime_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    polyMesh pMesh
    (
        io,
        xferCopy(pts),
        xferMove(faces),
        xferMove(owner),
        xferMove(neighbour)
    );

    List<polyPatch*> patches(patchStarts.size());

    forAll (patches, p)
    {
        patches[p] = new polyPatch
        (
            patchNames[p],
            patchSizes[p],
            patchStarts[p],
            p,
            pMesh.boundaryMesh()
        );
    }

    pMesh.addPatches(patches);

    timeCheck();

    labelHashSet wrongFaces(pMesh.nFaces()/100);

    motionSmoother::checkMesh
    (
        false,
        pMesh,
        cvMeshControls().cvMeshDict().subDict("meshQualityControls"),
        wrongFaces
    );

    const cellList& cells = pMesh.cells();

    forAll(cells, cI)
    {
        if (cells[cI].size() < 4)
        {
            Info<< "cell " << cI
                << " has " << cells[cI].size() << " faces."
                << endl;
        }
    }

    // forAllConstIter(labelHashSet, wrongFaces, iter)
    // {
    //     label faceI = iter.key();

    //     Info<< faceI << " " << pMesh.faces()[faceI] << endl;
    // }

    return wrongFaces.size();

    // For parallel running:
    // return returnReduce
    // (
    //     wrongFaces.size(),
    //     sumOp<label>()
    // );
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
                label own = -1;
                label nei = -1;

                if (ownerAndNeighbour(vA, vB, own, nei))
                {
                    reverse(newDualFace);
                }

                if (nei == -1)
                {
                    // boundary face

                    point ptA = topoint(vA->point());
                    point ptB = topoint(vB->point());

                    label patchIndex = geometryToConformTo_.findPatch(ptA, ptB);

                    if (patchIndex == -1)
                    {
                        patchIndex = patchNames.size() - 1;

                        // WarningIn("Foam::conformalVoronoiMesh::calcDualMesh")
                        //     << "Dual face found between Dv pair " << nl
                        //     << "    " << ptA << nl
                        //     << "    " << ptB << nl
                        //     << "    that is not on a surface patch. Adding to "
                        //     << patchNames[patchIndex]
                        //     << endl;
                    }

                    patchFaces[patchIndex].append(newDualFace);
                    patchOwners[patchIndex].append(own);
                }
                else
                {
                    // internal face

                    faces[dualFaceI] = newDualFace;
                    owner[dualFaceI] = own;
                    neighbour[dualFaceI] = nei;

                    dualFaceI++;
                }

                // --> OLD, FOR REFERENCE
                // label dualCellIndexA = vA->index();

                // if (!vA->internalOrBoundaryPoint())
                // {
                //     dualCellIndexA = -1;
                // }

                // label dualCellIndexB = vB->index();

                // if (!vB->internalOrBoundaryPoint())
                // {
                //     dualCellIndexB = -1;
                // }

                // label own = -1;
                // label nei = -1;

                // if (dualCellIndexA == -1 && dualCellIndexB == -1)
                // {
                //     FatalErrorIn
                //     (
                //         "void Foam::conformalVoronoiMesh::"
                //         "createFacesOwnerNeighbourAndPatches"
                //         "("
                //             "faceList& faces,"
                //             "labelList& owner,"
                //             "labelList& neighbour,"
                //             "wordList& patchNames,"
                //             "labelList& patchSizes,"
                //             "labelList& patchStarts,"
                //             "bool includeEmptyPatches"
                //         ") const"
                //     )
                //         << "Attempting to create a face joining "
                //         << "two unindexed dual cells "
                //         << exit(FatalError);
                // }
                // else if (dualCellIndexA == -1 || dualCellIndexB == -1)
                // {
                //     // boundary face, find which is the owner

                //     if (dualCellIndexA == -1)
                //     {
                //         own = dualCellIndexB;

                //         // reverse face order to correctly orientate normal
                //         reverse(newDualFace);
                //     }
                //     else
                //     {
                //         own = dualCellIndexA;
                //     }

                //     // Find which patch this face is on by finding the
                //     // intersection with the surface of the Delaunay edge
                //     // generating the face and identify the region of the
                //     // intersection.

                //     point ptA = topoint(vA->point());

                //     point ptB = topoint(vB->point());

                //     label patchIndex = geometryToConformTo_.findPatch(ptA, ptB);

                //     if (patchIndex == -1)
                //     {
                //         patchIndex = patchNames.size() - 1;

                //         // WarningIn("Foam::conformalVoronoiMesh::calcDualMesh")
                //         //     << "Dual face found between Dv pair " << nl
                //         //     << "    " << ptA << nl
                //         //     << "    " << ptB << nl
                //         //     << "    that is not on a surface patch. Adding to "
                //         //     << patchNames[patchIndex]
                //         //     << endl;
                //     }

                //     patchFaces[patchIndex].append(newDualFace);
                //     patchOwners[patchIndex].append(own);
                // }
                // else
                // {
                //     // internal face, find the lower cell to be the owner

                //     if (dualCellIndexB > dualCellIndexA)
                //     {
                //         own = dualCellIndexA;
                //         nei = dualCellIndexB;
                //     }
                //     else
                //     {
                //         own = dualCellIndexB;
                //         nei = dualCellIndexA;

                //         // reverse face order to correctly orientate normal
                //         reverse(newDualFace);
                //     }

                //     faces[dualFaceI] = newDualFace;

                //     owner[dualFaceI] = own;

                //     neighbour[dualFaceI] = nei;

                //     dualFaceI++;
                // }
                // <-- OLD, FOR REFERENCE
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

    Info<< "        Removing " << pts.size() - pointI << " unused points"
        << endl;

    pts.setSize(pointI);

    // Renumber the faces to use the new point numbers

    forAll(faces, fI)
    {
        inplaceRenumber(oldToNew, faces[fI]);
    }
}


void Foam::conformalVoronoiMesh::removeUnusedCells
(
    labelList& owner,
    labelList& neighbour
) const
{
    Info<< nl << "    Removing unused cells after filtering" << endl;

    PackedBoolList cellUsed(max(max(owner), max(neighbour)), false);

    // Scan all faces to find all of the cells that are used

    forAll(owner, oI)
    {
        cellUsed[owner[oI]] = true;
    }

    forAll(neighbour, nI)
    {
        cellUsed[neighbour[nI]] = true;
    }

    // Find all of the unused cells, create a list of them, then
    // subtract one from each owner and neighbour entry for each of
    // the unused cell indices that it is above.

    DynamicList<label> unusedCells;

    forAll(cellUsed, cUI)
    {
        if (cellUsed[cUI] == false)
        {
            unusedCells.append(cUI);
        }
    }

    if (unusedCells.size() > 0)
    {
        Info<< "        Removing " << unusedCells.size() <<  " unused cells "
            << endl;

        forAll(owner, oI)
        {
            label& o = owner[oI];

            o -= findLower(unusedCells, o) + 1;
        }

        forAll(neighbour, nI)
        {
            label& n = neighbour[nI];

            n -= findLower(unusedCells, n) + 1;
        }
    }
}


// ************************************************************************* //
