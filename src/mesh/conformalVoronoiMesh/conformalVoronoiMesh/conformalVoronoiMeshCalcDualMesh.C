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
#include "motionSmoother.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::calcDualMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchTypes,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    labelList& procNeighbours,
    pointField& cellCentres,
    bool filterFaces
)
{
    timeCheck("Start calcDualMesh");

    setVertexSizeAndAlignment();

    timeCheck("After setVertexSizeAndAlignment");

    // Make all filterCount values zero

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        cit->filterCount() = 0;
    }

    PackedBoolList boundaryPts;

    boundaryPts.reserve(number_of_cells());

    indexDualVertices(points, boundaryPts);

    {
        // Ideally requires a no-risk face filtering to get rid of zero area
        // faces and establish if the mesh can be produced at all to the
        // specified criteria

        Info<< nl << "Merging close points" << endl;

        // There is no guarantee that a merge of close points is no-risk
        mergeCloseDualVertices(points, boundaryPts);
    }

    timeCheck("After initial close point merge");

    if (filterFaces)
    {
        label nInitialBadQualityFaces = checkPolyMeshQuality(points).size();

        Info<< nl << "Initial check before face collapse, found "
            << nInitialBadQualityFaces << " bad quality faces"
            << endl;

        HashSet<labelPair, labelPair::Hash<> > deferredCollapseFaces;

        if (nInitialBadQualityFaces > 0)
        {
            Info<< nl
                << "A mesh could not be produced to satisfy the specified "
                << "quality criteria." << nl
                << "The quality and the surface conformation controls "
                << "can be altered and the " << nl
                << "internalDelaunayVertices read in to try again, or more "
                << "cell size resolution " << nl
                << "and motion iterations can be applied in areas where "
                << "problems are occurring."
                << endl;
        }

        if
        (
            nInitialBadQualityFaces == 0
         || cvMeshControls().continueFilteringOnBadInitialPolyMesh()
        )
        {
            label nBadQualityFaces = 0;

            labelHashSet lastWrongFaces;

            label nConsecutiveEqualFaceSets = 0;

            do
            {
                // Reindexing the Delaunay cells and regenerating the
                // points resets the mesh to the starting condition.

                indexDualVertices(points, boundaryPts);

                {
                    Info<< nl << "Merging close points" << endl;

                    mergeCloseDualVertices(points, boundaryPts);
                }

                {
                    // Risky and undo-able face filtering to reduce
                    // the face count as much as possible, staying
                    // within the specified criteria

                    Info<< nl << "Smoothing surface" << endl;

                    smoothSurface(points, boundaryPts);

                    Info<< nl << "Collapsing unnecessary faces" << endl;

                    collapseFaces(points, boundaryPts, deferredCollapseFaces);
                }

                labelHashSet wrongFaces = checkPolyMeshQuality(points);

                nBadQualityFaces = wrongFaces.size();

                Info<< nl << "Found " << nBadQualityFaces
                    << " bad quality faces" << endl;

                if (lastWrongFaces == wrongFaces)
                {
                    Info<< nl << "Consecutive iterations found the same set "
                        << "of bad quality faces." << endl;

                    if
                    (
                        ++nConsecutiveEqualFaceSets
                     >= cvMeshControls().maxConsecutiveEqualFaceSets()
                    )
                    {
                        Info<< nl << nConsecutiveEqualFaceSets
                            << " consecutive iterations produced the same "
                            << "bad quality faceSet, stopping filtering"
                            << endl;

                        break;
                    }
                }
                else
                {
                    nConsecutiveEqualFaceSets = 0;

                    lastWrongFaces = wrongFaces;
                }

                timeCheck("End of filtering iteration");

            } while (nBadQualityFaces > nInitialBadQualityFaces);
        }
    }

    // Final dual face and owner neighbour construction

    timeCheck("Before createFacesOwnerNeighbourAndPatches");

    createFacesOwnerNeighbourAndPatches
    (
        faces,
        owner,
        neighbour,
        patchTypes,
        patchNames,
        patchSizes,
        patchStarts,
        procNeighbours,
        false
    );

    // deferredCollapseFaceSet(owner, neighbour, deferredCollapseFaces);

    createCellCentres(cellCentres);

    removeUnusedCells(owner, neighbour, cellCentres);

    removeUnusedPoints(faces, points);

    timeCheck("End of calcDualMesh");
}


// void Foam::conformalVoronoiMesh::calcTetMesh
// (
//     pointField& points,
//     faceList& faces,
//     labelList& owner,
//     labelList& neighbour,
//     wordList& patchNames,
//     labelList& patchSizes,
//     labelList& patchStarts
// )
// {
//     labelList vertexMap(number_of_vertices());

//     label vertI = 0;

//     points.setSize(number_of_vertices());

//     for
//     (
//         Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
//         vit != finite_vertices_end();
//         ++vit
//     )
//     {
//         if (vit->internalPoint() || vit->pairPoint())
//         {
//             vertexMap[vit->index()] = vertI;
//             points[vertI] = topoint(vit->point());
//             vertI++;
//         }
//     }

//     points.setSize(vertI);

//     label cellI = 0;

//     for
//     (
//         Delaunay::Finite_cells_iterator cit = finite_cells_begin();
//         cit != finite_cells_end();
//         ++cit
//     )
//     {
//         if
//         (
//             cit->vertex(0)->internalOrBoundaryPoint()
//          || cit->vertex(1)->internalOrBoundaryPoint()
//          || cit->vertex(2)->internalOrBoundaryPoint()
//          || cit->vertex(3)->internalOrBoundaryPoint()
//         )
//         {
//              cit->cellIndex() = cellI++;
//         }
//         else
//         {
//             cit->cellIndex() = -1;
//         }
//     }

//     patchNames = geometryToConformTo_.patchNames();

//     patchNames.setSize(patchNames.size() + 1);

//     patchNames[patchNames.size() - 1] = "cvMesh_defaultPatch";

//     label nPatches = patchNames.size();

//     List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

//     List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

//     faces.setSize(number_of_facets());

//     owner.setSize(number_of_facets());

//     neighbour.setSize(number_of_facets());

//     label faceI = 0;

//     labelList verticesOnTriFace(3, -1);

//     face newFace(verticesOnTriFace);

//     for
//     (
//         Delaunay::Finite_facets_iterator fit = finite_facets_begin();
//         fit != finite_facets_end();
//         ++fit
//     )
//     {
//         const Cell_handle c1(fit->first);
//         const int oppositeVertex = fit->second;
//         const Cell_handle c2(c1->neighbor(oppositeVertex));

//         label c1I = c1->cellIndex();
//         label c2I = c2->cellIndex();

//         label ownerCell = -1;
//         label neighbourCell = -1;

//         if (c1I == -1 && c2I == -1)
//         {
//             // Both tets are outside, skip
//             continue;
//         }

//         for (label i = 0; i < 3; i++)
//         {
//             verticesOnTriFace[i] = vertexMap
//             [
//                 c1->vertex(vertex_triple_index(oppositeVertex, i))->index()
//             ];
//         }

//         newFace = face(verticesOnTriFace);

//         if (c1I == -1 || c2I == -1)
//         {
//             // Boundary face...
//             if (c1I == -1)
//             {
//                 //... with c1 outside
//                 ownerCell = c2I;
//             }
//             else
//             {
//                 // ... with c2 outside
//                 ownerCell = c1I;

//                 reverse(newFace);
//             }

//             label patchIndex = geometryToConformTo_.findPatch
//             (
//                 newFace.centre(points)
//             );

//             if (patchIndex == -1)
//             {
//                 patchIndex = patchNames.size() - 1;

//                 WarningIn("Foam::conformalVoronoiMesh::calcTetMesh")
//                     << "Tet face centre at  " << nl
//                     << newFace.centre(points) << nl
//                     << "did not find a surface patch. Adding to "
//                     << patchNames[patchIndex]
//                     << endl;
//             }

//             patchFaces[patchIndex].append(newFace);
//             patchOwners[patchIndex].append(ownerCell);
//         }
//         else
//         {
//             // Internal face...
//             if (c1I < c2I)
//             {
//                 // ...with c1 as the ownerCell
//                 ownerCell = c1I;
//                 neighbourCell = c2I;

//                 reverse(newFace);
//             }
//             else
//             {
//                 // ...with c2 as the ownerCell
//                 ownerCell = c2I;
//                 neighbourCell = c1I;
//             }

//             faces[faceI] = newFace;
//             owner[faceI] = ownerCell;
//             neighbour[faceI] = neighbourCell;
//             faceI++;
//         }
//     }

//     label nInternalFaces = faceI;

//     faces.setSize(nInternalFaces);
//     owner.setSize(nInternalFaces);
//     neighbour.setSize(nInternalFaces);

//     sortFaces(faces, owner, neighbour);

//     addPatches
//     (
//         nInternalFaces,
//         faces,
//         owner,
//         patchNames,
//         patchSizes,
//         patchStarts,
//         patchFaces,
//         patchOwners,
//         false
//     );
// }


void Foam::conformalVoronoiMesh::mergeCloseDualVertices
(
    const pointField& pts,
    const PackedBoolList& boundaryPts
)
{
    // Assess close points to be merged

    label nPtsMerged = 0;

    do
    {
        Map<label> dualPtIndexMap;

        nPtsMerged = mergeCloseDualVertices(pts, boundaryPts, dualPtIndexMap);

        if (nPtsMerged > 0)
        {
            Pout<< "    Merged " << nPtsMerged << " points " << endl;
        }

        reindexDualVertices(dualPtIndexMap);

    } while (nPtsMerged > 0);
}


Foam::label Foam::conformalVoronoiMesh::mergeCloseDualVertices
(
    const pointField& pts,
    const PackedBoolList& boundaryPts,
    Map<label>& dualPtIndexMap
) const
{
    label nPtsMerged = 0;

    scalar closenessTolerance = cvMeshControls().mergeClosenessCoeff();

    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
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
                if (boundaryPts[c2I] == true)
                {
                    // If c2I is a boundary point, then it is kept.
                    // If both are boundary points then c2I is chosen
                    // arbitrarily to be kept.

                    dualPtIndexMap.insert(c1I, c2I);
                    dualPtIndexMap.insert(c2I, c2I);
                }
                else
                {
                    dualPtIndexMap.insert(c1I, c1I);
                    dualPtIndexMap.insert(c2I, c1I);
                }

                nPtsMerged++;
            }
        }
    }

    return nPtsMerged;
}


void Foam::conformalVoronoiMesh::smoothSurface
(
    pointField& pts,
    const PackedBoolList& boundaryPts
)
{
    label nCollapsedFaces = 0;

    label iterI = 0;

    do
    {
        Map<label> dualPtIndexMap;

        nCollapsedFaces = smoothSurfaceDualFaces
        (
            pts,
            boundaryPts,
            dualPtIndexMap
        );

        reindexDualVertices(dualPtIndexMap);

        mergeCloseDualVertices(pts, boundaryPts);

        if (nCollapsedFaces > 0)
        {
            Info<< "    Collapsed " << nCollapsedFaces << " boundary faces"
                << endl;
        }

        if (++iterI > cvMeshControls().maxCollapseIterations())
        {
            Info<< "    maxCollapseIterations reached, stopping collapse"
                << endl;

            break;
        }

    } while (nCollapsedFaces > 0);

    // Force all points of boundary faces to be on the surface

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        label ptI = cit->cellIndex();

        label fC = cit->filterCount();

        if (fC > cvMeshControls().filterCountSkipThreshold())
        {
            // This vertex has been limited too many times, skip
            continue;
        }

        // Only cells with indices > -1 are valid
        if (ptI > -1)
        {
            if (boundaryPts[ptI] == true)
            {
                Foam::point& pt = pts[ptI];

                pointIndexHit surfHit;
                label hitSurface;

                geometryToConformTo_.findSurfaceNearest
                (
                    pt,
                    sqr(GREAT),
                    surfHit,
                    hitSurface
                );

                if (surfHit.hit())
                {
                    pt +=
                        (surfHit.hitPoint() - pt)
                       *pow(cvMeshControls().filterErrorReductionCoeff(), fC);
                }
            }
        }
    }

    mergeCloseDualVertices(pts, boundaryPts);
}


Foam::label Foam::conformalVoronoiMesh::smoothSurfaceDualFaces
(
    pointField& pts,
    const PackedBoolList& boundaryPts,
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
        Delaunay::Finite_edges_iterator eit = finite_edges_begin();
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

            label maxFC = maxFilterCount(eit);

            if (maxFC > cvMeshControls().filterCountSkipThreshold())
            {
                // A vertex on this face has been limited too many
                // times, skip
                continue;
            }

            pointIndexHit surfHit;
            label hitSurface;

            geometryToConformTo_.findSurfaceNearest
            (
                dualFace.centre(pts),
                sqr(GREAT),
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
                    boundaryPts,
                    dualPtIndexMap,
                    targetFaceSize,
                    GREAT,
                    maxFC
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
    const PackedBoolList& boundaryPts,
    HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
)
{
    label nCollapsedFaces = 0;

    label iterI = 0;

    do
    {
        Map<label> dualPtIndexMap;

        deferredCollapseFaces.clear();

        nCollapsedFaces = collapseFaces
        (
            pts,
            boundaryPts,
            dualPtIndexMap,
            deferredCollapseFaces
        );

        reindexDualVertices(dualPtIndexMap);

        mergeCloseDualVertices(pts, boundaryPts);

        if (nCollapsedFaces > 0)
        {
            Info<< "    Collapsed " << nCollapsedFaces << " faces" << endl;
            // Info<< "dualPtIndexMap" << nl << dualPtIndexMap << endl;
        }

        if (++iterI > cvMeshControls().maxCollapseIterations())
        {
            Info<< nl << "maxCollapseIterations reached, stopping collapse"
                << endl;

            break;
        }

    } while (nCollapsedFaces > 0);
}


Foam::label Foam::conformalVoronoiMesh::collapseFaces
(
    pointField& pts,
    const PackedBoolList& boundaryPts,
    Map<label>& dualPtIndexMap,
    HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
) const
{
    label nCollapsedFaces = 0;

    scalar collapseSizeLimitCoeff = cvMeshControls().filterSizeCoeff();

    for
    (
        Delaunay::Finite_edges_iterator eit = finite_edges_begin();
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

            label maxFC = maxFilterCount(eit);

            if (maxFC > cvMeshControls().filterCountSkipThreshold())
            {
                // A vertex on this face has been limited too many
                // times, skip
                continue;
            }

            scalar targetFaceSize = averageAnyCellSize(vA, vB);

            faceCollapseMode mode = collapseFace
            (
                dualFace,
                pts,
                boundaryPts,
                dualPtIndexMap,
                targetFaceSize,
                collapseSizeLimitCoeff,
                maxFC
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
    const PackedBoolList& boundaryPts,
    Map<label>& dualPtIndexMap,
    scalar targetFaceSize,
    scalar collapseSizeLimitCoeff,
    label maxFC
) const
{
    bool limitToQuadsOrTris = false;

    bool allowEarlyCollapseToPoint = true;

    // if (maxFC > cvMeshControls().filterCountSkipThreshold() - 3)
    // {
    //     limitToQuadsOrTris = true;

    //     allowEarlyCollapseToPoint = false;
    // }

    collapseSizeLimitCoeff *= pow
    (
        cvMeshControls().filterErrorReductionCoeff(),
        maxFC
    );

    labelList facePts(f);

    const Foam::point fC = f.centre(pts);

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

    scalar aspectRatio = 1.0;

    if (detJ < 1e-5)
    {
        collapseAxis = f.edges()[longestEdge(f, pts)].vec(pts);

        collapseAxis /= mag(collapseAxis);

        // Empirical correlation for high aspect ratio faces

        aspectRatio = sqrt(0.35/detJ);
    }
    else
    {
        vector eVals = eigenValues(J);

        if (mag(eVals.y() - eVals.x()) < 100*SMALL)
        {
            // First two eigenvalues are the same: i.e. a square face

            // Cannot necessarily determine linearly independent
            // eigenvectors, or any at all, use longest edge direction.

            collapseAxis = f.edges()[longestEdge(f, pts)].vec(pts);

            collapseAxis /= mag(collapseAxis);

            aspectRatio = 1.0;
        }
        else
        {
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
    }

    if (magSqr(collapseAxis) < VSMALL)
    {
        WarningIn
        (
            "Foam::conformalVoronoiMesh::collapseFace"
        )
            << "No collapse axis found for face, not collapsing."
            << endl;

        // Output face and collapse axis for visualisation

        Info<< "# Aspect ratio = " << aspectRatio << nl
            << "# inertia = " << J << nl
            << "# determinant = " << detJ << nl
            << "# eigenvalues = " << eigenValues(J) << nl
            << "# collapseAxis = " << collapseAxis << nl
            << "# facePts = " << facePts << nl
            << endl;

        forAll(f, fPtI)
        {
            meshTools::writeOBJ(Info, pts[f[fPtI]]);
        }

        Info<< "f";

        forAll(f, fPtI)
        {
            Info << " " << fPtI + 1;
        }

        Info<< endl;

        return fcmNone;
    }

    // The signed distance along the collapse axis passing through the
    // face centre that each vertex projects to.

    Field<scalar> d(f.size());

    forAll(f, fPtI)
    {
        const Foam::point& pt = pts[f[fPtI]];

        d[fPtI] = (collapseAxis & (pt - fC));
    }

    // Sort the projected distances and the corresponding vertex
    // indices along the collapse axis

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

    // Negative half
    SubList<scalar> dNeg(d, middle, 0);
    SubList<label> facePtsNeg(facePts, middle, 0);

    // Positive half
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
            // Negative half

            label collapseToPtI = facePtsNeg.first();

            Foam::point collapseToPt =
                collapseAxis*(sum(dNeg)/dNeg.size() - dShift) + fC;

            forAll(facePtsNeg, fPtI)
            {
                if (boundaryPts[facePtsNeg[fPtI]] == true)
                {
                    // If there is a point which is on the boundary,
                    // use it as the point to collapse others to, will
                    // use the first boundary point encountered if
                    // there are multiple boundary points.

                    collapseToPtI = facePtsNeg[fPtI];

                    collapseToPt = pts[collapseToPtI];

                    break;
                }
            }

            // ...otherwise arbitrarily choosing the most distant
            // point as the index to collapse to.

            forAll(facePtsNeg, fPtI)
            {
                dualPtIndexMap.insert(facePtsNeg[fPtI], collapseToPtI);
            }

            pts[collapseToPtI] = collapseToPt;

            // Positive half

            collapseToPtI = facePtsPos.last();

            collapseToPt = collapseAxis*(sum(dPos)/dPos.size() - dShift) + fC;

            forAll(facePtsPos, fPtI)
            {
                if (boundaryPts[facePtsPos[fPtI]] == true)
                {
                    // If there is a point which is on the boundary,
                    // use it as the point to collapse others to, will
                    // use the first boundary point encountered if
                    // there are multiple boundary points.

                    collapseToPtI = facePtsPos[fPtI];

                    collapseToPt = pts[collapseToPtI];

                    break;
                }
            }

            // ...otherwise arbitrarily choosing the most distant
            // point as the index to collapse to.

            forAll(facePtsPos, fPtI)
            {
                dualPtIndexMap.insert(facePtsPos[fPtI], collapseToPtI);
            }

            pts[collapseToPtI] = collapseToPt;

            break;
        }

        case fcmPoint:
        {
            label collapseToPtI = facePts.first();

            Foam::point collapseToPt = fC;

            forAll(facePts, fPtI)
            {
                if (boundaryPts[facePts[fPtI]] == true)
                {
                    // If there is a point which is on the boundary,
                    // use it as the point to collapse others to, will
                    // use the first boundary point encountered if
                    // there are multiple boundary points.

                    collapseToPtI = facePts[fPtI];

                    collapseToPt = pts[collapseToPtI];

                    break;
                }
            }

            // ...otherwise arbitrarily choosing the first point as
            // the index to collapse to.  Collapse to the face centre.

            forAll(facePts, fPtI)
            {
                dualPtIndexMap.insert(facePts[fPtI], collapseToPtI);
            }

            pts[collapseToPtI] = collapseToPt;

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


Foam::label Foam::conformalVoronoiMesh::longestEdge
(
    const face& f,
    const pointField& pts
) const
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

    return longestEdgeI;
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


Foam::labelHashSet Foam::conformalVoronoiMesh::checkPolyMeshQuality
(
    const pointField& pts
) const
{
    faceList faces;
    labelList owner;
    labelList neighbour;
    wordList patchTypes;
    wordList patchNames;
    labelList patchSizes;
    labelList patchStarts;
    labelList procNeighbours;
    pointField cellCentres;

    timeCheck("Start of checkPolyMeshQuality");

    Info<< nl << "Creating polyMesh to assess quality" << endl;

    createFacesOwnerNeighbourAndPatches
    (
        faces,
        owner,
        neighbour,
        patchTypes,
        patchNames,
        patchSizes,
        patchStarts,
        procNeighbours,
        false
    );

    createCellCentres(cellCentres);

    removeUnusedCells(owner, neighbour, cellCentres);

    polyMesh pMesh
    (
        IOobject
        (
            "cvMesh_temporary",
            runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        xferCopy(pts),
        xferMove(faces),
        xferMove(owner),
        xferMove(neighbour)
    );

    List<polyPatch*> patches(patchStarts.size());

    label nValidPatches = 0;

    forAll(patches, p)
    {
        if (patchTypes[p] == processorPolyPatch::typeName)
        {
            // Do not create empty processor patches

            if (patchSizes[p] > 0)
            {
                patches[nValidPatches] = new processorPolyPatch
                (
                    patchNames[p],
                    patchSizes[p],
                    patchStarts[p],
                    nValidPatches,
                    pMesh.boundaryMesh(),
                    Pstream::myProcNo(),
                    procNeighbours[p]
                );

                nValidPatches++;
            }
        }
        else
        {
            patches[nValidPatches] = polyPatch::New
            (
                patchTypes[p],
                patchNames[p],
                patchSizes[p],
                patchStarts[p],
                nValidPatches,
                pMesh.boundaryMesh()
            ).ptr();

            nValidPatches++;
        }
    }

    patches.setSize(nValidPatches);

    pMesh.addPatches(patches);

    // forAll(patches, p)
    // {
    //     patches[p] = new polyPatch
    //     (
    //         patchNames[p],
    //         patchSizes[p],
    //         patchStarts[p],
    //         p,
    //         pMesh.boundaryMesh()
    //     );
    // }

    // pMesh.addPatches(patches, false);

    // pMesh.overrideCellCentres(cellCentres);

    timeCheck("polyMesh created, checking quality");

    labelHashSet wrongFaces(pMesh.nFaces()/100);

    Info << endl;

    DynamicList<label> checkFaces(pMesh.nFaces());

    const vectorField& fAreas = pMesh.faceAreas();

    scalar faceAreaLimit = SMALL;

    forAll(fAreas, fI)
    {
        if (mag(fAreas[fI]) > faceAreaLimit)
        {
            checkFaces.append(fI);
        }
    }

    if (checkFaces.size() < fAreas.size())
    {
        Info<< "Excluding " << fAreas.size() - checkFaces.size()
            << " faces from check, < " << faceAreaLimit << " area" << endl;
    }

    motionSmoother::checkMesh
    (
        false,
        pMesh,
        cvMeshControls().cvMeshDict().subDict("meshQualityControls"),
        checkFaces,
        wrongFaces
    );

    {
        // Check for cells with more than 1 but fewer than 4 faces
        label nInvalidPolyhedra = 0;

        const cellList& cells = pMesh.cells();

        forAll(cells, cI)
        {
            if (cells[cI].size() < 4 && cells[cI].size() > 0)
            {
                // Info<< "cell " << cI << " " << cells[cI]
                //     << " has " << cells[cI].size() << " faces."
                //     << endl;

                nInvalidPolyhedra++;

                forAll(cells[cI], cFI)
                {
                    wrongFaces.insert(cells[cI][cFI]);
                }
            }
        }

        Info<< "    cells with more than 1 but fewer than 4 faces          : "
            << nInvalidPolyhedra << endl;

        // Check for cells with one internal face only

        labelList nInternalFaces(pMesh.nCells(), 0);

        for (label fI = 0; fI < pMesh.nInternalFaces(); fI++)
        {
            nInternalFaces[pMesh.faceOwner()[fI]]++;
            nInternalFaces[pMesh.faceNeighbour()[fI]]++;
        }

        const polyBoundaryMesh& patches = pMesh.boundaryMesh();

        forAll(patches, patchI)
        {
            if (patches[patchI].coupled())
            {
                const labelUList& owners = patches[patchI].faceCells();

                forAll(owners, i)
                {
                    nInternalFaces[owners[i]]++;
                }
            }
        }

        label oneInternalFaceCells = 0;

        forAll(nInternalFaces, cI)
        {
            if (nInternalFaces[cI] <= 1)
            {
                oneInternalFaceCells++;

                forAll(cells[cI], cFI)
                {
                    wrongFaces.insert(cells[cI][cFI]);
                }
            }
        }

        Info<< "    cells with with zero or one non-boundary face          : "
            << oneInternalFaceCells << endl;
    }


    PackedBoolList ptToBeLimited(pts.size(), false);

    forAllConstIter(labelHashSet, wrongFaces, iter)
    {
        const face f = pMesh.faces()[iter.key()];

        forAll(f, fPtI)
        {
            ptToBeLimited[f[fPtI]] = true;
        }
    }

    // // Limit connected cells

    // labelHashSet limitCells(pMesh.nCells()/100);

    // const labelListList& ptCells = pMesh.pointCells();

    // forAllConstIter(labelHashSet, wrongFaces, iter)
    // {
    //     const face f = pMesh.faces()[iter.key()];

    //     forAll(f, fPtI)
    //     {
    //         label ptI = f[fPtI];

    //         const labelList& pC = ptCells[ptI];

    //         forAll(pC, pCI)
    //         {
    //             limitCells.insert(pC[pCI]);
    //         }
    //     }
    // }

    // const labelListList& cellPts = pMesh.cellPoints();

    // forAllConstIter(labelHashSet, limitCells, iter)
    // {
    //     label cellI = iter.key();

    //     const labelList& cP = cellPts[cellI];

    //     forAll(cP, cPI)
    //     {
    //         ptToBeLimited[cP[cPI]] = true;
    //     }
    // }


    // Apply Delaunay cell filterCounts and determine the maximum
    // overall filterCount

    label maxFilterCount = 0;

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        label cI = cit->cellIndex();

        if (cI >= 0)
        {
            if (ptToBeLimited[cI] == true)
            {
                cit->filterCount()++;
            }

            if (cit->filterCount() > maxFilterCount)
            {
                maxFilterCount = cit->filterCount();
            }
        }
    }

    Info<< nl << "Maximum number of filter limits applied: "
        << maxFilterCount << endl;

    return wrongFaces;
}


void Foam::conformalVoronoiMesh::indexDualVertices
(
    pointField& pts,
    PackedBoolList& boundaryPts
)
{
    // Indexing Delaunay cells, which are the dual vertices

    label dualVertI = 0;

    pts.setSize(number_of_cells());

    boundaryPts.clear();

    boundaryPts.setSize(number_of_cells(), false);

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
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

            pts[dualVertI] = topoint(dual(cit));

            if
            (
                !cit->vertex(0)->internalOrBoundaryPoint()
             || !cit->vertex(1)->internalOrBoundaryPoint()
             || !cit->vertex(2)->internalOrBoundaryPoint()
             || !cit->vertex(3)->internalOrBoundaryPoint()
            )
            {
                // This is a boundary dual vertex

                boundaryPts[dualVertI] = true;
            }

            dualVertI++;
        }
        else
        {
            cit->cellIndex() = -1;
        }
    }

    pts.setSize(dualVertI);

    boundaryPts.setSize(dualVertI);
}


void Foam::conformalVoronoiMesh::reindexDualVertices
(
    const Map<label>& dualPtIndexMap
)
{
    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
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
    wordList& patchTypes,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    labelList& procNeighbours,
    bool includeEmptyPatches
) const
{
    patchNames = geometryToConformTo_.patchNames();
    patchTypes.setSize(patchNames.size() + 1, wallPolyPatch::typeName);
    procNeighbours.setSize(patchNames.size() + 1, -1);

    patchNames.setSize(patchNames.size() + 1);

    label defaultPatchIndex = patchNames.size() - 1;

    patchTypes[defaultPatchIndex] = wallPolyPatch::typeName;
    procNeighbours[defaultPatchIndex] = -1;
    patchNames[defaultPatchIndex] = "cvMesh_defaultPatch";

    label nProcPatches = 0;

    if (Pstream::parRun())
    {
        boolList procUsed(Pstream::nProcs(), false);

        // Determine which processor patches are required
        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->referred())
            {
                procUsed[vit->procIndex()] = true;
            }
        }

        forAll(procUsed, pUI)
        {
            if (procUsed[pUI])
            {
                nProcPatches++;
            }
        }

        label nNonProcPatches = patchNames.size();

        patchTypes.setSize(nNonProcPatches + nProcPatches);
        procNeighbours.setSize(nNonProcPatches + nProcPatches, -1);
        patchNames.setSize(nNonProcPatches + nProcPatches);

        label procAddI = 0;

        forAll(procUsed, pUI)
        {
            if (procUsed[pUI])
            {
                patchTypes[nNonProcPatches + procAddI] =
                    processorPolyPatch::typeName;

                patchNames[nNonProcPatches + procAddI] =
                    "procBoundary"
                   + name(Pstream::myProcNo())
                   + "to"
                   + name(pUI);

                procNeighbours[nNonProcPatches + procAddI] = pUI;

                procAddI++;
            }
        }
    }

    // Pout<< patchTypes << " " << patchNames << endl;

    label nPatches = patchNames.size();

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_edges());

    owner.setSize(number_of_edges());

    neighbour.setSize(number_of_edges());

    List<Pair<DynamicList<label> > > procPatchSortingIndex(nPatches);

    label dualFaceI = 0;

    for
    (
        Delaunay::Finite_edges_iterator eit = finite_edges_begin();
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

                    Foam::point ptA = topoint(vA->point());
                    Foam::point ptB = topoint(vB->point());

                    label patchIndex = -1;

                    if (vA->referredInternal() || vB->referredInternal())
                    {
                        // One (and only one) of the points is an internal
                        // point from another processor

                        label procIndex = max(vA->procIndex(), vB->procIndex());

                        patchIndex = findIndex(procNeighbours, procIndex);

                        // The lower processor index is the owner of the
                        // two for the purposed of sorting the patch faces.

                        if (Pstream::myProcNo() < procIndex)
                        {
                            // Use this processor's vertex index as the master
                            // for sorting

                           Pair<DynamicList<label> >& sortingIndex =
                               procPatchSortingIndex[patchIndex];

                            if (vB->referredInternal())
                            {
                                sortingIndex.first().append(vA->index());
                                sortingIndex.second().append(vB->index());
                            }
                            else
                            {
                                sortingIndex.first().append(vB->index());
                                sortingIndex.second().append(vA->index());
                            }
                        }
                        else
                        {
                            // Use the other processor's vertex index as the
                            // master for sorting

                            Pair<DynamicList<label> >& sortingIndex =
                                procPatchSortingIndex[patchIndex];

                            if (vA->referredInternal())
                            {
                                sortingIndex.first().append(vA->index());
                                sortingIndex.second().append(vB->index());
                            }
                            else
                            {
                                sortingIndex.first().append(vB->index());
                                sortingIndex.second().append(vA->index());
                            }
                        }

                        // Pout<< ptA << " " << ptB
                        //     << " proc indices "
                        //     << vA->procIndex() << " " << vB->procIndex()
                        //     << " indices " << vA->index()
                        //     << " " << vB->index()
                        //     << " my proc " << Pstream::myProcNo()
                        //     << " addedIndex "
                        //     << procPatchSortingIndex[patchIndex].last()
                        //     << endl;
                    }
                    else
                    {
                        patchIndex = geometryToConformTo_.findPatch(ptA, ptB);
                    }

                    if (patchIndex == -1)
                    {
                        // Did not find a surface patch between
                        // between Dv pair, adding face to default
                        // patch

                        patchIndex = defaultPatchIndex;
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
            }
        }
    }

    if (!patchFaces[defaultPatchIndex].empty())
    {
        Pout<< nl << patchFaces[defaultPatchIndex].size()
            << " faces were not able to have their patch determined from "
            << "the surface. "
            << nl <<  "Adding to patch " << patchNames[defaultPatchIndex]
            << endl;
    }

    label nInternalFaces = dualFaceI;

    faces.setSize(nInternalFaces);
    owner.setSize(nInternalFaces);
    neighbour.setSize(nInternalFaces);

    timeCheck("polyMesh quality checked");

    sortFaces(faces, owner, neighbour);

    sortProcPatches
    (
        patchFaces,
        patchOwners,
        procPatchSortingIndex
    );

    timeCheck("faces, owner, neighbour sorted");

    addPatches
    (
        nInternalFaces,
        faces,
        owner,
        patchSizes,
        patchStarts,
        patchFaces,
        patchOwners
    );
}


void Foam::conformalVoronoiMesh::createCellCentres
(
    pointField& cellCentres
) const
{
    cellCentres.setSize(number_of_vertices());

    label vertI = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            cellCentres[vit->index()] = topoint(vit->point());

            vertI++;
        }
    }

    cellCentres.setSize(vertI);
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
        << "Sorting faces, owner and neighbour into upper triangular order"
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

    label blockStart = 0;

    for (label i = 1; i < owner.size(); i++)
    {
        label blockLength = -1;

        if (owner[i] > owner[i - 1])
        {
            blockLength = i - blockStart;
        }
        else if (i == owner.size() - 1)
        {
            // If the last element is not a jump in owner, then it
            // needs to trigger a sort of the last block, but with a
            // block length that is one element longer so that it
            // sorts itself.

            // If it is a jump in owner, then it will form a block of
            // length one, and so will not need sorted.

            blockLength = i - blockStart + 1;
        }

        if (blockLength >= 1)
        {
            labelList blockIndices = identity(blockLength) + blockStart;

            SubList<label> neighbourBlock
            (
                neighbour,
                blockLength,
                blockStart
            );

            sortedOrder(neighbourBlock, blockIndices);

            blockIndices = invert(blockIndices.size(), blockIndices);

            forAll(blockIndices, b)
            {
                oldToNew[blockStart + b] = blockIndices[b] + blockStart;
            }

            blockStart = i;
        }
    }

    // owner does not need re-sorted
    inplaceReorder(oldToNew, faces);
    inplaceReorder(oldToNew, neighbour);
}


void Foam::conformalVoronoiMesh::sortProcPatches
(
    List<DynamicList<face> >& patchFaces,
    List<DynamicList<label> >& patchOwners,
    List<Pair<DynamicList<label> > >& patchSortingIndices
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    forAll(patchSortingIndices, patchI)
    {
        faceList& faces = patchFaces[patchI];
        labelList& owner = patchOwners[patchI];

        Pair<DynamicList<label> >& sortingIndices = patchSortingIndices[patchI];

        List<label>& primary = sortingIndices.first();
        List<label>& secondary = sortingIndices.second();

        if (!primary.empty())
        {
            if
            (
                faces.size() != primary.size()
             || owner.size() != primary.size()
            )
            {
                FatalErrorIn
                (
                    "void Foam::conformalVoronoiMesh::sortProcPatches"
                    "("
                        "List<DynamicList<face> >& patchFaces, "
                        "List<DynamicList<label> >& patchOwners, "
                        "const List<DynamicList<label> >& patchSortingIndices"
                    ") const"
                )
                    << "patch size and size of sorting indices is inconsistent "
                    << " for patch " << patchI << nl
                    << " faces.size() " << faces.size() << nl
                    << " owner.size() " << owner.size() << nl
                    << " sortingIndices.first().size() "
                    << sortingIndices.first().size()
                    << exit(FatalError) << endl;
            }

            // Two stage sort:
            // 1) sort by primary

            labelList oldToNew;

            sortedOrder(primary, oldToNew);

            oldToNew = invert(oldToNew.size(), oldToNew);

            inplaceReorder(oldToNew, primary);
            inplaceReorder(oldToNew, secondary);
            inplaceReorder(oldToNew, faces);
            inplaceReorder(oldToNew, owner);

            // 2) in each block of primary sort by secondary

            // Reset map.  Elements that are not sorted will retain their -1
            // value, which will mean that they are ignored by inplaceReorder

            oldToNew = -1;

            label blockStart = 0;

            for (label i = 1; i < primary.size(); i++)
            {
                label blockLength = -1;

                if (primary[i] > primary[i - 1])
                {
                    blockLength = i - blockStart;
                }
                else if (i == primary.size() - 1)
                {
                    // If the last element is not a jump in index, then it
                    // needs to trigger a sort of the last block, but with a
                    // block length that is one element longer so that it
                    // sorts itself.

                    // If it is a jump in index, then it will form a block of
                    // length one, and so will not need sorted.

                    blockLength = i - blockStart + 1;
                }

                if (blockLength >= 1)
                {
                    labelList blockIndices = identity(blockLength) + blockStart;

                    SubList<label> secondaryBlock
                    (
                        secondary,
                        blockLength,
                        blockStart
                    );

                    sortedOrder(secondaryBlock, blockIndices);

                    blockIndices = invert(blockIndices.size(), blockIndices);

                    forAll(blockIndices, b)
                    {
                        oldToNew[blockStart + b] = blockIndices[b] + blockStart;
                    }

                    blockStart = i;
                }
            }

            inplaceReorder(oldToNew, faces);
            inplaceReorder(oldToNew, owner);
        }
    }
}


void Foam::conformalVoronoiMesh::addPatches
(
    const label nInternalFaces,
    faceList& faces,
    labelList& owner,
    labelList& patchSizes,
    labelList& patchStarts,
    List<DynamicList<face> >& patchFaces,
    List<DynamicList<label> >& patchOwners
) const
{
    label nPatches = patchFaces.size();

    patchSizes.setSize(nPatches, -1);
    patchStarts.setSize(nPatches, -1);

    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        patchSizes[p] = patchFaces[p].size();
        patchStarts[p] = nInternalFaces + nBoundaryFaces;

        nBoundaryFaces += patchSizes[p];
    }

    faces.setSize(nInternalFaces + nBoundaryFaces);
    owner.setSize(nInternalFaces + nBoundaryFaces);

    label faceI = nInternalFaces;

    forAll(patchFaces, p)
    {
        forAll(patchFaces[p], f)
        {
            faces[faceI] = patchFaces[p][f];
            owner[faceI] = patchOwners[p][f];

            faceI++;
        }
    }
}


void Foam::conformalVoronoiMesh::removeUnusedPoints
(
    faceList& faces,
    pointField& pts
) const
{
    Info<< nl << "Removing unused points" << endl;

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

    // Move all of the used points to the start of the pointField and
    // truncate it

    forAll(ptUsed, ptUI)
    {
        if (ptUsed[ptUI] == true)
        {
            oldToNew[ptUI] = pointI++;
        }
    }

    inplaceReorder(oldToNew, pts);

    Info<< "    Removing " << pts.size() - pointI << " unused points"
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
    labelList& neighbour,
    pointField& cellCentres
) const
{
    Info<< nl << "Removing unused cells" << endl;

    PackedBoolList cellUsed(number_of_vertices(), false);

    // Scan all faces to find all of the cells that are used

    forAll(owner, oI)
    {
        cellUsed[owner[oI]] = true;
    }

    forAll(neighbour, nI)
    {
        cellUsed[neighbour[nI]] = true;
    }

    label cellI = 0;

    labelList oldToNew(cellUsed.size(), -1);

    // Move all of the used cellCentres to the start of the pointField and
    // truncate it

    forAll(cellUsed, cellUI)
    {
        if (cellUsed[cellUI] == true)
        {
            oldToNew[cellUI] = cellI++;
        }
    }

    inplaceReorder(oldToNew, cellCentres);

    cellCentres.setSize(cellI);

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
        // Pout<< "Removing " << unusedCells.size() <<  " unused cell labels"
        //     << endl;

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
