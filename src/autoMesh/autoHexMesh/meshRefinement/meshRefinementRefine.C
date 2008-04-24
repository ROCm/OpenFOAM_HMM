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

Class
    meshRefinement

\*----------------------------------------------------------------------------*/

#include "meshRefinement.H"
#include "trackedParticle.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "faceSet.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "polyTopoChange.H"
#include "mapDistributePolyMesh.H"
#include "featureEdgeMesh.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Get faces (on the new mesh) that have in some way been affected by the
// mesh change. Picks up all faces but those that are between two
// unrefined faces. (Note that of an unchanged face the edge still might be
// split but that does not change any face centre or cell centre.
Foam::labelList Foam::meshRefinement::getChangedFaces
(
    const mapPolyMesh& map,
    const labelList& oldCellsToRefine
)
{
    const polyMesh& mesh = map.mesh();

    labelList changedFaces;
    {
        // Mark any face on a cell which has been added or changed
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Note that refining a face changes the face centre (for a warped face)
        // which changes the cell centre. This again changes the cellcentre-
        // cellcentre edge across all faces of the cell.
        // Note: this does not happen for unwarped faces but unfortunately
        // we don't have this information.

        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighbour = mesh.faceNeighbour();
        const cellList& cells = mesh.cells();
        const label nInternalFaces = mesh.nInternalFaces();

        // Mark refined cells on old mesh
        PackedList<1> oldRefineCell(map.nOldCells(), 0u);

        forAll(oldCellsToRefine, i)
        {
            oldRefineCell.set(oldCellsToRefine[i], 1u);
        }

        // Mark refined faces
        PackedList<1> refinedInternalFace(nInternalFaces, 0u);

        // 1. Internal faces

        for (label faceI = 0; faceI < nInternalFaces; faceI++)
        {
            label oldOwn = map.cellMap()[faceOwner[faceI]];
            label oldNei = map.cellMap()[faceNeighbour[faceI]];

            if
            (
                oldOwn >= 0
             && oldRefineCell.get(oldOwn) == 0u
             && oldNei >= 0
             && oldRefineCell.get(oldNei) == 0u
            )
            {
                // Unaffected face since both neighbours were not refined.
            }
            else
            {
                refinedInternalFace.set(faceI, 1u);
            }
        }


        // 2. Boundary faces

        boolList refinedBoundaryFace(mesh.nFaces()-nInternalFaces, false);

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];

            label faceI = pp.start();

            forAll(pp, i)
            {
                label oldOwn = map.cellMap()[faceOwner[faceI]];

                if (oldOwn >= 0 && oldRefineCell.get(oldOwn) == 0u)
                {
                    // owner did exist and wasn't refined.
                }
                else
                {
                    refinedBoundaryFace[faceI-nInternalFaces] = true;
                }
                faceI++;
            }
        }

        // Synchronise refined face status
        syncTools::syncBoundaryFaceList
        (
            mesh,
            refinedBoundaryFace,
            orEqOp<bool>(),
            false
        );


        // 3. Mark all faces affected by refinement. Refined faces are in
        //    - refinedInternalFace
        //    - refinedBoundaryFace
        boolList changedFace(mesh.nFaces(), false);

        forAll(refinedInternalFace, faceI)
        {
            if (refinedInternalFace.get(faceI) == 1u)
            {
                const cell& ownFaces = cells[faceOwner[faceI]];
                forAll(ownFaces, ownI)
                {
                    changedFace[ownFaces[ownI]] = true;
                }
                const cell& neiFaces = cells[faceNeighbour[faceI]];
                forAll(neiFaces, neiI)
                {
                    changedFace[neiFaces[neiI]] = true;
                }
            }
        }

        forAll(refinedBoundaryFace, i)
        {
            if (refinedBoundaryFace[i])
            {
                const cell& ownFaces = cells[faceOwner[i+nInternalFaces]];
                forAll(ownFaces, ownI)
                {
                    changedFace[ownFaces[ownI]] = true;
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh,
            changedFace,
            orEqOp<bool>(),
            false
        );


        // Now we have in changedFace marked all affected faces. Pack.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        label nChanged = 0;

        forAll(changedFace, faceI)
        {
            if (changedFace[faceI])
            {
                nChanged++;
            }
        }

        changedFaces.setSize(nChanged);
        nChanged = 0;

        forAll(changedFace, faceI)
        {
            if (changedFace[faceI])
            {
                changedFaces[nChanged++] = faceI;
            }
        }
    }

    label nChangedFaces = changedFaces.size();
    reduce(nChangedFaces, sumOp<label>());

    if (debug)
    {
        Pout<< "getChangedFaces : Detected "
            << " local:" << changedFaces.size()
            << " global:" << nChangedFaces
            << " changed faces out of " << mesh.globalData().nTotalFaces()
            << endl;

        faceSet changedFacesSet(mesh, "changedFaces", changedFaces);
        Pout<< "getChangedFaces : Writing " << changedFaces.size()
            << " changed faces to faceSet " << changedFacesSet.name()
            << endl;
        changedFacesSet.write();
    }

    return changedFaces;
}


// Mark cell for refinement (if not already marked). Return false if
// refinelimit hit. Keeps running count (in nRefine) of cells marked for
// refinement
bool Foam::meshRefinement::markForRefine
(
    const label markValue,
    const label nAllowRefine,

    label& cellValue,
    label& nRefine
)
{
    if (cellValue == -1)
    {
        cellValue = markValue;
        nRefine++;
    }

    return nRefine <= nAllowRefine;
}


// Calculates list of cells to refine based on intersection with feature edge.
Foam::label Foam::meshRefinement::markFeatureRefinement
(
    const point& keepPoint,
    const PtrList<featureEdgeMesh>& featureMeshes,
    const labelList& featureLevels,
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    // We want to refine all cells containing a feature edge.
    // - don't want to search over whole mesh
    // - don't want to build octree for whole mesh
    // - so use tracking to follow the feature edges
    //
    // 1. find non-manifold points on feature edges (i.e. start of feature edge
    //    or 'knot')
    // 2. seed particle starting at keepPoint going to this non-manifold point
    // 3. track particles to their non-manifold point
    //
    // 4. track particles across their connected feature edges, marking all
    //    visited cells with their level (through trackingData)
    // 5. do 4 until all edges have been visited.


    // Find all start cells of features. Is done by tracking from keepPoint.
    Cloud<trackedParticle> cloud(mesh_, IDLList<trackedParticle>());

    // Create particles on whichever processor holds the keepPoint.
    label cellI = mesh_.findCell(keepPoint);

    if (cellI != -1)
    {
        forAll(featureMeshes, featI)
        {
            const featureEdgeMesh& featureMesh = featureMeshes[featI];
            const labelListList& pointEdges = featureMesh.pointEdges();

            forAll(pointEdges, pointI)
            {
                if (pointEdges[pointI].size() != 2)
                {
                    if (debug)
                    {
                        Pout<< "Adding particle from point:" << pointI
                            << " coord:" << featureMesh.points()[pointI]
                            << " pEdges:" << pointEdges[pointI]
                            << endl;
                    }

                    // Non-manifold point. Create particle.
                    cloud.addParticle
                    (
                        new trackedParticle
                        (
                            cloud,
                            keepPoint,
                            cellI,
                            featureMesh.points()[pointI],   // endpos
                            featureLevels[featI],           // level
                            featI,                          // featureMesh
                            pointI                          // end point
                        )
                    );
                }
            }
        }
    }


    // Largest refinement level of any feature passed through
    labelList maxFeatureLevel(mesh_.nCells(), -1);

    // Database to pass into trackedParticle::move
    trackedParticle::trackData td(cloud, maxFeatureLevel);

    // Track all particles to their end position.
    cloud.move(td);

    // Reset level
    maxFeatureLevel = -1;

    // Whether edge has been visited.
    List<PackedList<1> > featureEdgeVisited(featureMeshes.size());

    forAll(featureMeshes, featI)
    {
        featureEdgeVisited[featI].setSize(featureMeshes[featI].edges().size());
        featureEdgeVisited[featI] = 0u;
    }

    while (true)
    {
        label nParticles = 0;

        // Make particle follow edge.
        forAllIter(Cloud<trackedParticle>, cloud, iter)
        {
            trackedParticle& tp = iter();

            label featI = tp.i();
            label pointI = tp.j();

            const featureEdgeMesh& featureMesh = featureMeshes[featI];
            const labelList& pEdges = featureMesh.pointEdges()[pointI];

            // Particle now at pointI. Check connected edges to see which one
            // we have to visit now.

            bool keepParticle = false;

            forAll(pEdges, i)
            {
                label edgeI = pEdges[i];

                if (featureEdgeVisited[featI].get(edgeI) == 0)
                {
                    // Unvisited edge. Make the particle go to the other point
                    // on the edge.

                    const edge& e = featureMesh.edges()[edgeI];
                    label otherPointI = e.otherVertex(pointI);

                    featureEdgeVisited[featI].set(edgeI, 1u);
                    tp.end() = featureMesh.points()[otherPointI];
                    tp.j() = otherPointI;
                    keepParticle = true;
                    break;
                }
            }

            if (!keepParticle)
            {
                // Particle at 'knot' where another particle already has been
                // seeded. Delete particle.
                cloud.deleteParticle(tp);
            }
            else
            {
                // Keep particle
                nParticles++;
            }
        }

        reduce(nParticles, sumOp<label>());
        if (nParticles == 0)
        {
            break;
        }

        // Track all particles to their end position.
        cloud.move(td);
    }


    // See which cells to refine. maxFeatureLevel will hold highest level
    // of any feature edge that passed through.

    const labelList& cellLevel = meshCutter_.cellLevel();

    label oldNRefine = nRefine;

    forAll(maxFeatureLevel, cellI)
    {
        if (maxFeatureLevel[cellI] > cellLevel[cellI])
        {
            // Mark
            if
            (
               !markForRefine
                (
                    0,                      // surface (n/a)
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                )
            )
            {
                // Reached limit
                break;
            }
        }
    }

    return returnReduce(nRefine-oldNRefine,  sumOp<label>());
}


// Mark cells for non-surface intersection based refinement.
Foam::label Foam::meshRefinement::markInternalRefinement
(
    const PtrList<searchableSurface>& shells,
    const labelList& shellLevels,
    const boolList& shellRefineInside,
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();


    label oldNRefine = nRefine;

    // Number of cells marked for refinement per shell.
    labelList nCellsPerShell(shells.size(), 0);

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            // Cell not marked for refinement. Check if inside any shell
            // with higher refinement level than cell currently has.

            bool reachedLimit = false;

            forAll(shells, shellI)
            {
                // Cached inside-outside from tree
                searchableSurface::volumeType t =
                    shells[shellI].getVolumeType(cellCentres[cellI]);

                //// Uncached inside-outside from treeData
                //label t = shells[shellI].shapes().getVolumeType
                //    (
                //        shells[shellI],
                //        cellCentres[cellI]
                //    );

                // Which side of shell is to be refined
                searchableSurface::volumeType refSide =
                (
                    shellRefineInside[shellI]
                  ? searchableSurface::INSIDE
                  : searchableSurface::OUTSIDE
                );

                if (t == refSide && cellLevel[cellI] < shellLevels[shellI])
                {
                    // Cell is inside shell with higher refinement level. Mark
                    // for refinement.

                    reachedLimit = !markForRefine
                    (
                        labelMax,
                        nAllowRefine,
                        refineCell[cellI],
                        nRefine
                    );

                    if (reachedLimit)
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining internal cells"
                                << " since reaching my cell limit of "
                                << mesh_.nCells()+7*nRefine << endl;
                        }
                        break;
                    }
                    else
                    {
                        // Cell successfully marked for refinement
                        nCellsPerShell[shellI]++;
                    }
                }
            }

            if (reachedLimit)
            {
                break;
            }
        }
    }

    Pstream::listCombineGather(nCellsPerShell, plusEqOp<label>());
    Pstream::listCombineScatter(nCellsPerShell);

    Info<< "Marked for refinement per shell :" << endl;
    forAll(nCellsPerShell, shellI)
    {
        Info<< "    shell:" << shellI << " nCells:" << nCellsPerShell[shellI]
            << nl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Mark cells for surface intersection based refinement.
Foam::label Foam::meshRefinement::markSurfaceRefinement
(
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // Use cached surfaceIndex_ to detect if any intersection. If so
    // re-intersect to determine level wanted.

    label faceI;
    for (faceI = 0; faceI < surfaceIndex_.size(); faceI++)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            label own = mesh_.faceOwner()[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                // Test if not both sides already marked for refinement.
                if (refineCell[own] == -1 || refineCell[nei] == -1)
                {
                    pointIndexHit hit;
                    label surfI = surfaces_.findHigherIntersection
                    (
                        cellCentres[own],
                        cellCentres[nei],
                        min(cellLevel[own], cellLevel[nei]),
                        hit
                    );

                    if (surfI != -1)
                    {
                        // Found intersection with surface with higher wanted
                        // refinement. Mark cell for refining.
                        // Note: could do optimization here and only refine the
                        // side of the face the intersection actually goes
                        // through.
                        // This would require us though to find all
                        // intersections
                        // first instead of now only the first higher one. Also
                        // would need the exact intersection of face with cc-cc
                        // connection?

                        label surfaceMinLevel =
                            surfaces_.minLevelField(surfI)[hit.index()];

                        if (surfaceMinLevel > cellLevel[own])
                        {
                            // Owner needs refining
                            if
                            (
                               !markForRefine
                                (
                                    surfI,
                                    nAllowRefine,
                                    refineCell[own],
                                    nRefine
                                )
                            )
                            {
                                break;
                            }
                        }

                        if (surfaceMinLevel > cellLevel[nei])
                        {
                            // Neighbour needs refining
                            if
                            (
                               !markForRefine
                                (
                                    surfI,
                                    nAllowRefine,
                                    refineCell[nei],
                                    nRefine
                                )
                            )
                            {
                                break;
                            }
                        }
                    }
                }
            }
            else if (refineCell[own] != -1)
            {
                // boundary face with unmarked owner

                label bFaceI = faceI - mesh_.nInternalFaces();

                pointIndexHit hit;
                label surfI = surfaces_.findHigherIntersection
                (
                    cellCentres[own],
                    neiCc[bFaceI],
                    min(cellLevel[own], neiLevel[bFaceI]),
                    hit
                );

                if (surfI != -1)
                {
                    if
                    (
                       !markForRefine
                        (
                            surfI,
                            nAllowRefine,
                            refineCell[own],
                            nRefine
                        )
                    )
                    {
                        break;
                    }
                }
            }
        }
    }

    if (faceI < surfaceIndex_.size())
    {
        if (debug)
        {
            Pout<< "Stopped refining since reaching my cell limit of "
                << mesh_.nCells()+7*nRefine << endl;
        }
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Given intersection of (face of) cell by newSurfI, newTriI, check whether
// it needs to be refined. Updates maxCellSurf, maxCellTri and
// refineCell,nRefine if it decides to refine the cell. Returns false
// if the nRefine limit has been reached, true otherwise.
bool Foam::meshRefinement::checkCurvature
(
    const labelList& globalToPatch,
    const scalar curvature,
    const bool markDifferingRegions,
    const label nAllowRefine,

    const label newSurfI,
    const label newTriI,

    const label cellI,

    label& maxCellSurfI,
    label& maxCellTriI,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    if (maxCellSurfI == -1)
    {
        // First visit of cell. Store
        maxCellSurfI = newSurfI;
        maxCellTriI = newTriI;
    }
    else
    {
        // Second or more visit.
        label cellRegion = surfaces_.triangleRegion(maxCellSurfI, maxCellTriI);
        label newRegion = surfaces_.triangleRegion(newSurfI, newTriI);

        // Update max
        label maxLevel = surfaces_.maxLevel()[cellRegion];

        if (surfaces_.maxLevel()[newRegion] > maxLevel)
        {
            maxCellSurfI = newSurfI;
            maxCellTriI = newTriI;
            maxLevel = surfaces_.maxLevel()[newRegion];
        }


        // Check if cell is candidate for refinement
        if (cellLevel[cellI] < maxLevel)
        {
            // Test 1: different regions
            if (markDifferingRegions && cellRegion != newRegion)
            {
                return markForRefine
                (
                    globalToPatch[newRegion],   // mark with non-neg number.
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );
            }

            // Test 2: different normals
            const vector& cellN =
                surfaces_[maxCellSurfI].faceNormals()[maxCellTriI];
            const vector& newN =
                surfaces_[newSurfI].faceNormals()[newTriI];

            if ((cellN & newN) < curvature)
            {
                return markForRefine
                (
                    globalToPatch[newRegion],   // mark with non-neg number.
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );
            }
        }
    }

    // Did not reach refinement limit.
    return true;
}


// Mark cells for surface curvature based refinement.
Foam::label Foam::meshRefinement::markSurfaceCurvatureRefinement
(
    const labelList& globalToPatch,
    const scalar curvature,
    const bool markDifferingRegions,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // 1. Any cell on more than one surface gets refined (if its current level
    // is <= max of the surface max level)

    // 2. Any cell on only one surface with a neighbour on a different surface
    // gets refined (if its current level etc.)



    // Index of surface with the max maxLevel and the actual triangle
    // on the surface. We store these and not e.g. global region so we can get
    // at:
    //  - global region
    //  - global patch
    //  - face normal
    labelList cellMaxSurface(mesh_.nCells(), -1);
    labelList cellMaxTriangle(mesh_.nCells(), -1);

    // 1.

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            label own = mesh_.faceOwner()[faceI];

            // There is an intersection. Do more accurate test to get all
            // intersections
            labelList surfaceIndices;
            List<pointIndexHit> hits;

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                surfaces_.findAllIntersections
                (
                    cellCentres[own],
                    cellCentres[nei],
                    surfaceIndices,
                    hits
                );
            }
            else
            {
                label bFaceI = faceI - mesh_.nInternalFaces();

                surfaces_.findAllIntersections
                (
                    cellCentres[own],
                    neiCc[bFaceI],
                    surfaceIndices,
                    hits
                );
            }


            // See if intersection adds any new information to the owner or
            // neighbour cell.

            forAll(surfaceIndices, i)
            {
                label surfI = surfaceIndices[i];
                label triI = hits[i].index();

                checkCurvature
                (
                    globalToPatch,
                    curvature,
                    markDifferingRegions,
                    nAllowRefine,

                    surfI,
                    triI,

                    own,

                    cellMaxSurface[own],
                    cellMaxTriangle[own],

                    refineCell,
                    nRefine
                );

                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];

                    checkCurvature
                    (
                        globalToPatch,
                        curvature,
                        markDifferingRegions,
                        nAllowRefine,

                        surfI,
                        triI,

                        nei,

                        cellMaxSurface[nei],
                        cellMaxTriangle[nei],

                        refineCell,
                        nRefine
                    );
                }
            }
        }

        if (nRefine > nAllowRefine)
        {
            if (debug)
            {
                Pout<< "Stopped refining since reaching my cell limit of "
                    << mesh_.nCells()+7*nRefine << endl;
            }
            break;
        }
    }

    // 2. Find out a measure of surface curvature and region edges.
    // Send over surface region and surface normal to neighbour cell.

    // global region
    labelList ownRegion(mesh_.nFaces(), -1);
    labelList neiRegion(mesh_.nFaces(), -1);
    // local normal at hit
    vectorField ownNormal(mesh_.nFaces(), vector::zero);
    vectorField neiNormal(mesh_.nFaces(), vector::zero);

    // Internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];

        if (cellMaxSurface[own] != -1)
        {
            label surfI = cellMaxSurface[own];
            label triI = cellMaxTriangle[own];

            ownRegion[faceI] = surfaces_.triangleRegion(surfI, triI);
            ownNormal[faceI] = surfaces_[surfI].faceNormals()[triI];
        }

        label nei = mesh_.faceNeighbour()[faceI];

        if (cellMaxSurface[nei] != -1)
        {
            label surfI = cellMaxSurface[nei];
            label triI = cellMaxTriangle[nei];

            neiRegion[faceI] = surfaces_.triangleRegion(surfI, triI);
            neiNormal[faceI] = surfaces_[surfI].faceNormals()[triI];
        }
    }
    // Boundary faces
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];

        if (cellMaxSurface[own] != -1)
        {
            label surfI = cellMaxSurface[own];
            label triI = cellMaxTriangle[own];

            ownRegion[faceI] = surfaces_.triangleRegion(surfI, triI);
            ownNormal[faceI] = surfaces_[surfI].faceNormals()[triI];

            neiRegion[faceI] = ownRegion[faceI];
            neiNormal[faceI] = ownNormal[faceI];
        }
    }

    syncTools::swapFaceList(mesh_, neiRegion, false);
    syncTools::swapFaceList(mesh_, neiNormal, false);

    for (label faceI = 0; faceI < mesh_.nFaces(); faceI++)
    {
        if (ownRegion[faceI] != -1 && neiRegion[faceI] != -1)
        {
            if
            (
                (markDifferingRegions && (ownRegion[faceI] != neiRegion[faceI]))
             || ((ownNormal[faceI] & neiNormal[faceI]) < curvature)
            )
            {
                label own = mesh_.faceOwner()[faceI];

                if (cellLevel[own] < surfaces_.maxLevel()[ownRegion[faceI]])
                {
                    if
                    (
                        !markForRefine
                        (
                            globalToPatch[ownRegion[faceI]],
                            nAllowRefine,
                            refineCell[own],
                            nRefine
                        )
                    )
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining since reaching my cell"
                                << " limit of " << mesh_.nCells()+7*nRefine
                                << endl;
                        }
                        break;
                    }

                }

                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];

                    if (cellLevel[nei] < surfaces_.maxLevel()[neiRegion[faceI]])
                    {
                        if
                        (
                            !markForRefine
                            (
                                globalToPatch[neiRegion[faceI]],
                                nAllowRefine,

                                refineCell[nei],
                                nRefine
                            )
                        )
                        {
                            if (debug)
                            {
                                Pout<< "Stopped refining since reaching my cell"
                                    << " limit of " << mesh_.nCells()+7*nRefine
                                    << endl;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }


    return returnReduce(nRefine-oldNRefine, sumOp<label>());
    label totNRefined = returnReduce(totNRefined, sumOp<label>());

    return totNRefined;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate list of cells to refine. Gets for any edge (start - end)
// whether it intersects the surface. Does more accurate test and checks
// the wanted level on the surface intersected.
// Does approximate precalculation of how many cells can be refined before
// hitting overall limit maxGlobalCells.
Foam::labelList Foam::meshRefinement::refineCandidates
(
    const point& keepPoint,
    const labelList& globalToPatch,
    const scalar curvature,

    const PtrList<featureEdgeMesh>& featureMeshes,
    const labelList& featureLevels,

    const PtrList<searchableSurface>& shells,
    const labelList& shellLevels,
    const boolList& shellRefineInside,

    const bool featureRefinement,
    const bool internalRefinement,
    const bool surfaceRefinement,
    const bool curvatureRefinement,
    const label maxGlobalCells,
    const label maxLocalCells
) const
{
    label totNCells = mesh_.globalData().nTotalCells();

    labelList cellsToRefine;

    if (totNCells < maxGlobalCells)
    {
        // Every cell I refine adds 7 cells. Estimate the number of cells
        // I am allowed to refine.
        // Assume perfect distribution so can only refine as much the fraction
        // of the mesh I hold. This prediction step prevents us having to do
        // lots of reduces to keep count of the total number of cells selected
        // for refinement.

        //scalar fraction = scalar(mesh_.nCells())/totNCells;
        //label myMaxCells = label(maxGlobalCells*fraction);
        //label nAllowRefine = (myMaxCells - mesh_.nCells())/7;
        ////label nAllowRefine = (maxLocalCells - mesh_.nCells())/7;
        //
        //Pout<< "refineCandidates:" << nl
        //    << "    total cells:" << totNCells << nl
        //    << "    local cells:" << mesh_.nCells() << nl
        //    << "    local fraction:" << fraction << nl
        //    << "    local allowable cells:" << myMaxCells << nl
        //    << "    local allowable refinement:" << nAllowRefine << nl
        //    << endl;

        //- Disable refinement shortcut
        label nAllowRefine = labelMax;

        // Marked for refinement (>= 0) or not (-1). Actual value is the
        // index of the surface it intersects.
        labelList refineCell(mesh_.nCells(), -1);
        label nRefine = 0;


        // Swap neighbouring cell centres and cell level
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
        calcNeighbourData(neiLevel, neiCc);



        // Cells pierced by feature lines
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (featureRefinement)
        {
            label nFeatures = markFeatureRefinement
            (
                keepPoint,
                featureMeshes,
                featureLevels,
                nAllowRefine,

                refineCell,
                nRefine
            );

            Info<< "Marked for refinement due to explicit features    : "
                << nFeatures << " cells."  << endl;
        }

        // Inside refinement shells
        // ~~~~~~~~~~~~~~~~~~~~~~~~

        if (internalRefinement)
        {
            label nShell = markInternalRefinement
            (
                shells,
                shellLevels,
                shellRefineInside,
                nAllowRefine,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to refinement shells    : "
                << nShell << " cells."  << endl;
        }

        // Refinement based on intersection of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (surfaceRefinement)
        {
            label nSurf = markSurfaceRefinement
            (
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to surface intersection : "
                << nSurf << " cells."  << endl;
        }

        // Refinement based on curvature of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (curvatureRefinement)
        {
            label nCurv = markSurfaceCurvatureRefinement
            (
                globalToPatch,
                curvature,
                false,              // do not refine at multiple regions
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to curvature/regions    : "
                << nCurv << " cells."  << endl;
        }

        // Pack cells-to-refine
        // ~~~~~~~~~~~~~~~~~~~~

        cellsToRefine.setSize(nRefine);
        nRefine = 0;

        forAll(refineCell, cellI)
        {
            if (refineCell[cellI] != -1)
            {
                cellsToRefine[nRefine++] = cellI;
            }
        }
    }

    return cellsToRefine;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(mesh_);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    // Update fields
    mesh_.updateMesh(map);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    // Update intersection info
    updateMesh(map, getChangedFaces(map, cellsToRefine));

    return map;
}


// Do refinement of consistent set of cells followed by truncation and
// load balancing.
Foam::autoPtr<Foam::mapDistributePolyMesh>
 Foam::meshRefinement::refineAndBalance
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& cellsToRefine
)
{
    // Do all refinement
    refine(cellsToRefine);

    if (debug)
    {
        Pout<< "Writing refined but unbalanced " << msg
            << " mesh to time " << mesh_.time().timeName() << endl;
        write
        (
            debug,
            mesh_.time().path()
           /mesh_.time().timeName()
        );
        Pout<< "Dumped debug data in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // test all is still synced across proc patches
        checkData();
    }

    Info<< "Refined mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s" << endl;
    printMeshInfo(debug, "After refinement " + msg);


    //// Remove cells which are inside closed surfaces
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //if (findIndex(surfaces.closed(), true) != -1)
    //{
    //    Info<< "Removing cells fully inside closed surfaces."
    //        << endl;
    //    removeInsideCells
    //    (
    //        msg,               // refinement iteration for statistics only
    //        exposedFacesPatch  // patch to use for exposed internal faces
    //    );
    //    Info<< "Removed inside cells in = "
    //        << mesh_.time().cpuTimeIncrement() << " s" << endl;
    //    if (debug)
    //    {
    //       // test all is still synced across proc patches
    //       checkData();
    //    }
    //}


    // Load balancing
    // ~~~~~~~~~~~~~~

    autoPtr<mapDistributePolyMesh> distMap;

    if (Pstream::nProcs() > 1)
    {
        labelList distribution(decomposer.decompose(mesh_.cellCentres()));
        // Get distribution such that baffle faces stay internal to the
        // processor.
        //labelList distribution(decomposePreserveBaffles(decomposer));

        if (debug)
        {
            Pout<< "Wanted distribution:"
                << distributor.countCells(distribution)
                << endl;
        }
        // Do actual sending/receiving of mesh
        distMap = distributor.distribute(distribution);

        // Update numbering of meshRefiner
        distribute(distMap);

        Info<< "Balanced mesh in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        printMeshInfo(debug, "After balancing " + msg);


        if (debug)
        {
            Pout<< "Writing " << msg
                << " mesh to time " << mesh_.time().timeName() << endl;
            write
            (
                debug,
                mesh_.time().path()
               /mesh_.time().timeName()
            );
            Pout<< "Dumped debug data in = "
                << mesh_.time().cpuTimeIncrement() << " s" << endl;

            // test all is still synced across proc patches
            checkData();
        }
    }

    return distMap;
}


// ************************************************************************* //
