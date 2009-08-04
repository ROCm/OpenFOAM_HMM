/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "meshRefinement.H"
#include "fvMesh.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "pointSet.H"
#include "faceSet.H"
#include "indirectPrimitivePatch.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyModifyFace.H"
#include "polyModifyCell.H"
#include "polyAddFace.H"
#include "polyRemoveFace.H"
#include "polyAddPoint.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "OFstream.H"
#include "regionSplit.H"
#include "removeCells.H"
#include "motionSmoother.H"
#include "polyMeshGeometry.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Repatches external face or creates baffle for internal face
// with user specified patches (might be different for both sides).
// Returns label of added face.
Foam::label Foam::meshRefinement::createBaffle
(
    const label faceI,
    const label ownPatch,
    const label neiPatch,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[faceI];
    label zoneID = mesh_.faceZones().whichZone(faceI);
    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];
        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            f,                          // modified face
            faceI,                      // label of face
            mesh_.faceOwner()[faceI],   // owner
            -1,                         // neighbour
            false,                      // face flip
            ownPatch,                   // patch for face
            false,                      // remove from zone
            zoneID,                     // zone for face
            zoneFlip,                   // face flip in zone
            -1                          // sub patch
        )
    );


    label dupFaceI = -1;

    if (mesh_.isInternalFace(faceI))
    {
        if (neiPatch == -1)
        {
            FatalErrorIn
            (
                "meshRefinement::createBaffle"
                "(const label, const label, const label, polyTopoChange&)"
            )   << "No neighbour patch for internal face " << faceI
                << " fc:" << mesh_.faceCentres()[faceI]
                << " ownPatch:" << ownPatch << abort(FatalError);
        }

        dupFaceI = meshMod.setAction
        (
            polyAddFace
            (
                f.reverseFace(),            // modified face
                mesh_.faceNeighbour()[faceI],// owner
                -1,                         // neighbour
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                faceI,                      // masterFaceID,
                true,                       // face flip
                neiPatch,                   // patch for face
                zoneID,                     // zone for face
                zoneFlip,                   // face flip in zone
                -1                          // sub patch
            )
        );
    }
    return dupFaceI;
}


// Get an estimate for the patch the internal face should get. Bit heuristic.
Foam::label Foam::meshRefinement::getBafflePatch
(
    const labelList& facePatch,
    const label faceI
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Loop over face points
    // for each point check all faces patch IDs
    // as soon as an ID >= 0 is found, break and assign that ID
    // to the current face.
    // Check first for real patch (so proper surface intersection and then
    // in facePatch array for patches to block off faces

    forAll(mesh_.faces()[faceI], fp)
    {
        label pointI = mesh_.faces()[faceI][fp];

        forAll(mesh_.pointFaces()[pointI], pf)
        {
            label pFaceI = mesh_.pointFaces()[pointI][pf];

            label patchI = patches.whichPatch(pFaceI);

            if (patchI != -1 && !patches[patchI].coupled())
            {
                return patchI;
            }
            else if (facePatch[pFaceI] != -1)
            {
                return facePatch[pFaceI];
            }
        }
    }

    // Loop over owner and neighbour cells, looking for the first face with a
    // valid patch number
    const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

    forAll(ownFaces, i)
    {
        label cFaceI = ownFaces[i];

        label patchI = patches.whichPatch(cFaceI);

        if (patchI != -1 && !patches[patchI].coupled())
        {
            return patchI;
        }
        else if (facePatch[cFaceI] != -1)
        {
            return facePatch[cFaceI];
        }
    }

    if (mesh_.isInternalFace(faceI))
    {
        const cell& neiFaces = mesh_.cells()[mesh_.faceNeighbour()[faceI]];

        forAll(neiFaces, i)
        {
            label cFaceI = neiFaces[i];

            label patchI = patches.whichPatch(cFaceI);

            if (patchI != -1 && !patches[patchI].coupled())
            {
                return patchI;
            }
            else if (facePatch[cFaceI] != -1)
            {
                return facePatch[cFaceI];
            }
        }
    }

    WarningIn
    (
        "meshRefinement::getBafflePatch(const labelList& const label)"
    )   << "Could not find boundary face neighbouring internal face "
        << faceI << " with face centre " << mesh_.faceCentres()[faceI]
        << nl
        << "Using arbitrary patch " << 0 << " instead." << endl;

    return 0;
}


// Determine patches for baffles.
void Foam::meshRefinement::getBafflePatches
(
    const labelList& globalToPatch,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& ownPatch,
    labelList& neiPatch
) const
{
    autoPtr<OFstream> str;
    label vertI = 0;
    if (debug&OBJINTERSECTIONS)
    {
        str.reset(new OFstream(mesh_.time().path()/"intersections.obj"));

        Pout<< "getBafflePatches : Writing surface intersections to file "
            << str().name() << nl << endl;
    }

    const pointField& cellCentres = mesh_.cellCentres();

    // Build list of surfaces that are not to be baffled.
    const wordList& cellZoneNames = surfaces_.cellZoneNames();

    labelList surfacesToBaffle(cellZoneNames.size());
    label baffleI = 0;
    forAll(cellZoneNames, surfI)
    {
        if (cellZoneNames[surfI].size() > 0)
        {
            if (debug)
            {
                Pout<< "getBafflePatches : Not baffling surface "
                    << surfaces_.names()[surfI] << endl;
            }
        }
        else
        {
            surfacesToBaffle[baffleI++] = surfI;
        }
    }
    surfacesToBaffle.setSize(baffleI);

    ownPatch.setSize(mesh_.nFaces());
    ownPatch = -1;
    neiPatch.setSize(mesh_.nFaces());
    neiPatch = -1;


    // Collect candidate faces
    // ~~~~~~~~~~~~~~~~~~~~~~~

    labelList testFaces(intersectedFaces());

    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        label own = mesh_.faceOwner()[faceI];

        if (mesh_.isInternalFace(faceI))
        {
            start[i] = cellCentres[own];
            end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
        }
        else
        {
            start[i] = cellCentres[own];
            end[i] = neiCc[faceI-mesh_.nInternalFaces()];
        }
    }


    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    surfaces_.findNearestIntersection
    (
        surfacesToBaffle,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
    );

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        if (hit1[i].hit() && hit2[i].hit())
        {
            if (str.valid())
            {
                meshTools::writeOBJ(str(), start[i]);
                vertI++;
                meshTools::writeOBJ(str(), hit1[i].rawPoint());
                vertI++;
                meshTools::writeOBJ(str(), hit2[i].rawPoint());
                vertI++;
                meshTools::writeOBJ(str(), end[i]);
                vertI++;
                str()<< "l " << vertI-3 << ' ' << vertI-2 << nl;
                str()<< "l " << vertI-2 << ' ' << vertI-1 << nl;
                str()<< "l " << vertI-1 << ' ' << vertI << nl;
            }

            // Pick up the patches
            ownPatch[faceI] = globalToPatch
            [
                surfaces_.globalRegion(surface1[i], region1[i])
            ];
            neiPatch[faceI] = globalToPatch
            [
                surfaces_.globalRegion(surface2[i], region2[i])
            ];

            if (ownPatch[faceI] == -1 || neiPatch[faceI] == -1)
            {
                FatalErrorIn("getBafflePatches(..)")
                    << "problem." << abort(FatalError);
            }
        }
    }

    // No need to parallel sync since intersection data (surfaceIndex_ etc.)
    // already guaranteed to be synced...
    // However:
    // - owncc and neicc are reversed on different procs so might pick
    //   up different regions reversed? No problem. Neighbour on one processor
    //   might not be owner on the other processor but the neighbour is
    //   not used when creating baffles from proc faces.
    // - tolerances issues occasionally crop up.
    syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    syncTools::syncFaceList(mesh_, neiPatch, maxEqOp<label>());
}


// Create baffle for every face where ownPatch != -1
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createBaffles
(
    const labelList& ownPatch,
    const labelList& neiPatch
)
{
    if
    (
        ownPatch.size() != mesh_.nFaces()
     || neiPatch.size() != mesh_.nFaces()
    )
    {
        FatalErrorIn
        (
            "meshRefinement::createBaffles"
            "(const labelList&, const labelList&)"
        )   << "Illegal size :"
            << " ownPatch:" << ownPatch.size()
            << " neiPatch:" << neiPatch.size()
            << ". Should be number of faces:" << mesh_.nFaces()
            << abort(FatalError);
    }

    if (debug)
    {
        labelList syncedOwnPatch(ownPatch);
        syncTools::syncFaceList(mesh_, syncedOwnPatch, maxEqOp<label>());
        labelList syncedNeiPatch(neiPatch);
        syncTools::syncFaceList(mesh_, syncedNeiPatch, maxEqOp<label>());

        forAll(syncedOwnPatch, faceI)
        {
            if
            (
                (ownPatch[faceI] == -1 && syncedOwnPatch[faceI] != -1)
             || (neiPatch[faceI] == -1 && syncedNeiPatch[faceI] != -1)
            )
            {
                FatalErrorIn
                (
                    "meshRefinement::createBaffles"
                    "(const labelList&, const labelList&)"
                )   << "Non synchronised at face:" << faceI
                    << " on patch:" << mesh_.boundaryMesh().whichPatch(faceI)
                    << " fc:" << mesh_.faceCentres()[faceI] << endl
                    << "ownPatch:" << ownPatch[faceI]
                    << " syncedOwnPatch:" << syncedOwnPatch[faceI]
                    << " neiPatch:" << neiPatch[faceI]
                    << " syncedNeiPatch:" << syncedNeiPatch[faceI]
                    << abort(FatalError);
            }
        }
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nBaffles = 0;

    forAll(ownPatch, faceI)
    {
        if (ownPatch[faceI] != -1)
        {
            // Create baffle or repatch face. Return label of inserted baffle
            // face.
            createBaffle
            (
                faceI,
                ownPatch[faceI],   // owner side patch
                neiPatch[faceI],   // neighbour side patch
                meshMod
            );
            nBaffles++;
        }
    }

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    //- Redo the intersections on the newly create baffle faces. Note that
    //  this changes also the cell centre positions.
    faceSet baffledFacesSet(mesh_, "baffledFacesSet", 2*nBaffles);

    const labelList& reverseFaceMap = map().reverseFaceMap();
    const labelList& faceMap = map().faceMap();

    // Pick up owner side of baffle
    forAll(ownPatch, oldFaceI)
    {
        label faceI = reverseFaceMap[oldFaceI];

        if (ownPatch[oldFaceI] != -1 && faceI >= 0)
        {
            const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

            forAll(ownFaces, i)
            {
                baffledFacesSet.insert(ownFaces[i]);
            }
        }
    }
    // Pick up neighbour side of baffle (added faces)
    forAll(faceMap, faceI)
    {
        label oldFaceI = faceMap[faceI];

        if (oldFaceI >= 0 && reverseFaceMap[oldFaceI] != faceI)
        {
            const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

            forAll(ownFaces, i)
            {
                baffledFacesSet.insert(ownFaces[i]);
            }
        }
    }
    baffledFacesSet.sync(mesh_);

    updateMesh(map, baffledFacesSet.toc());

    return map;
}


void Foam::meshRefinement::markBoundaryFace
(
    const label faceI,
    boolList& isBoundaryFace,
    boolList& isBoundaryEdge,
    boolList& isBoundaryPoint
) const
{
    isBoundaryFace[faceI] = true;

    const labelList& fEdges = mesh_.faceEdges(faceI);

    forAll(fEdges, fp)
    {
        isBoundaryEdge[fEdges[fp]] = true;
    }

    const face& f = mesh_.faces()[faceI];

    forAll(f, fp)
    {
        isBoundaryPoint[f[fp]] = true;
    }
}


// Returns list with for every internal face -1 or the patch they should
// be baffled into. Gets run after createBaffles so all the surface
// intersections have already been turned into baffles. Used to remove cells
// by baffling all their faces and have the splitMeshRegions chuck away non
// used regions.
Foam::labelList Foam::meshRefinement::markFacesOnProblemCells
(
    const labelList& globalToPatch
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    // Per internal face (boundary faces not used) the patch that the
    // baffle should get (or -1)
    labelList facePatch(mesh_.nFaces(), -1);

    // Mark all points and edges on baffle patches (so not on any inlets,
    // outlets etc.)
    boolList isBoundaryPoint(mesh_.nPoints(), false);
    boolList isBoundaryEdge(mesh_.nEdges(), false);
    boolList isBoundaryFace(mesh_.nFaces(), false);

    // Fill boundary data. All elements on meshed patches get marked.
    // Get the labels of added patches.
    labelList adaptPatchIDs(meshRefinement::addedPatches(globalToPatch));

    forAll(adaptPatchIDs, i)
    {
        label patchI = adaptPatchIDs[i];

        const polyPatch& pp = patches[patchI];

        label faceI = pp.start();

        forAll(pp, j)
        {
            markBoundaryFace
            (
                faceI,
                isBoundaryFace,
                isBoundaryEdge,
                isBoundaryPoint
            );

            faceI++;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>()
    );


    // For each cell count the number of anchor points that are on
    // the boundary:
    // 8 : check the number of (baffle) boundary faces. If 3 or more block
    //     off the cell since the cell would get squeezed down to a diamond
    //     (probably; if the 3 or more faces are unrefined (only use the
    //      anchor points))
    // 7 : store. Used to check later on whether there are points with
    //     3 or more of these cells. (note that on a flat surface a boundary
    //     point will only have 4 cells connected to it)

    // Does cell have exactly 7 of its 8 anchor points on the boundary?
    PackedList<1> hasSevenBoundaryAnchorPoints(mesh_.nCells(), 0u);
    // If so what is the remaining non-boundary anchor point?
    labelHashSet nonBoundaryAnchors(mesh_.nCells()/10000);

    // On-the-fly addressing storage.
    DynamicList<label> dynFEdges;
    DynamicList<label> dynCPoints;

    // Count of faces marked for baffling
    label nBaffleFaces = 0;

    forAll(cellLevel, cellI)
    {
        const labelList& cPoints = mesh_.cellPoints(cellI, dynCPoints);

        // Get number of anchor points (pointLevel == cellLevel)

        label nBoundaryAnchors = 0;
        label nNonAnchorBoundary = 0;
        label nonBoundaryAnchor = -1;

        forAll(cPoints, i)
        {
            label pointI = cPoints[i];

            if (pointLevel[pointI] <= cellLevel[cellI])
            {
                // Anchor point
                if (isBoundaryPoint[pointI])
                {
                    nBoundaryAnchors++;
                }
                else
                {
                    // Anchor point which is not on the surface
                    nonBoundaryAnchor = pointI;
                }
            }
            else if (isBoundaryPoint[pointI])
            {
                nNonAnchorBoundary++;
            }
        }

        if (nBoundaryAnchors == 8)
        {
            const cell& cFaces = mesh_.cells()[cellI];

            // Count boundary faces.
            label nBfaces = 0;

            forAll(cFaces, cFaceI)
            {
                if (isBoundaryFace[cFaces[cFaceI]])
                {
                    nBfaces++;
                }
            }

            // If nBfaces > 1 make all non-boundary non-baffle faces baffles.
            // We assume that this situation is where there is a single
            // cell sticking out which would get flattened.

            // Eugene: delete cell no matter what.
            //if (nBfaces > 1)
            {
                forAll(cFaces, cf)
                {
                    label faceI = cFaces[cf];

                    if (facePatch[faceI] == -1 && mesh_.isInternalFace(faceI))
                    {
                        facePatch[faceI] = getBafflePatch(facePatch, faceI);
                        nBaffleFaces++;

                        // Mark face as a 'boundary'
                        markBoundaryFace
                        (
                            faceI,
                            isBoundaryFace,
                            isBoundaryEdge,
                            isBoundaryPoint
                        );
                    }
                }
            }
        }
        else if (nBoundaryAnchors == 7)
        {
            // Mark the cell. Store the (single!) non-boundary anchor point.
            hasSevenBoundaryAnchorPoints.set(cellI, 1u);
            nonBoundaryAnchors.insert(nonBoundaryAnchor);
        }
    }


    // Loop over all points. If a point is connected to 4 or more cells
    // with 7 anchor points on the boundary set those cell's non-boundary faces
    // to baffles

    DynamicList<label> dynPCells;

    forAllConstIter(labelHashSet, nonBoundaryAnchors, iter)
    {
        label pointI = iter.key();

        const labelList& pCells = mesh_.pointCells(pointI, dynPCells);

        // Count number of 'hasSevenBoundaryAnchorPoints' cells.
        label n = 0;

        forAll(pCells, i)
        {
            if (hasSevenBoundaryAnchorPoints.get(pCells[i]) == 1u)
            {
                n++;
            }
        }

        if (n > 3)
        {
            // Point in danger of being what? Remove all 7-cells.
            forAll(pCells, i)
            {
                label cellI = pCells[i];

                if (hasSevenBoundaryAnchorPoints.get(cellI) == 1u)
                {
                    const cell& cFaces = mesh_.cells()[cellI];

                    forAll(cFaces, cf)
                    {
                        label faceI = cFaces[cf];

                        if
                        (
                            facePatch[faceI] == -1
                         && mesh_.isInternalFace(faceI)
                        )
                        {
                            facePatch[faceI] = getBafflePatch(facePatch, faceI);
                            nBaffleFaces++;

                            // Mark face as a 'boundary'
                            markBoundaryFace
                            (
                                faceI,
                                isBoundaryFace,
                                isBoundaryEdge,
                                isBoundaryPoint
                            );
                        }
                    }
                }
            }
        }
    }


    // Sync all. (note that pointdata and facedata not used anymore but sync
    // anyway)

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>()
    );


    // Find faces with all edges on the boundary and make them baffles
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (facePatch[faceI] == -1)
        {
            const labelList& fEdges = mesh_.faceEdges(faceI, dynFEdges);
            label nFaceBoundaryEdges = 0;

            forAll(fEdges, fe)
            {
                if (isBoundaryEdge[fEdges[fe]])
                {
                    nFaceBoundaryEdges++;
                }
            }

            if (nFaceBoundaryEdges == fEdges.size())
            {
                facePatch[faceI] = getBafflePatch(facePatch, faceI);
                nBaffleFaces++;

                // Do NOT update boundary data since this would grow blocked
                // faces across gaps.
            }
        }
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                if (facePatch[faceI] == -1)
                {
                    const labelList& fEdges = mesh_.faceEdges(faceI, dynFEdges);
                    label nFaceBoundaryEdges = 0;

                    forAll(fEdges, fe)
                    {
                        if (isBoundaryEdge[fEdges[fe]])
                        {
                            nFaceBoundaryEdges++;
                        }
                    }

                    if (nFaceBoundaryEdges == fEdges.size())
                    {
                        facePatch[faceI] = getBafflePatch(facePatch, faceI);
                        nBaffleFaces++;

                        // Do NOT update boundary data since this would grow
                        // blocked faces across gaps.
                    }
                }

                faceI++;
            }
        }
    }

    Info<< "markFacesOnProblemCells : marked "
        << returnReduce(nBaffleFaces, sumOp<label>())
        << " additional internal faces to be converted into baffles."
        << endl;

    return facePatch;
}
//XXXXXXXXXXXXXX
// Mark faces to be baffled to prevent snapping problems. Does
// test to find nearest surface and checks which faces would get squashed.
Foam::labelList Foam::meshRefinement::markFacesOnProblemCellsGeometric
(
    const dictionary& motionDict,
    const labelList& globalToPatch
) const
{
    // Get the labels of added patches.
    labelList adaptPatchIDs(meshRefinement::addedPatches(globalToPatch));

    // Construct addressing engine.
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh_,
            adaptPatchIDs
        )
    );
    const indirectPrimitivePatch& pp = ppPtr();
    const pointField& localPoints = pp.localPoints();
    const labelList& meshPoints = pp.meshPoints();

    // Find nearest (non-baffle) surface
    pointField newPoints(mesh_.points());
    {
        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces_.findNearest
        (
            surfaces_.getUnnamedSurfaces(),
            localPoints,
            scalarField(localPoints.size(), sqr(GREAT)),    // sqr of attraction
            hitSurface,
            hitInfo
        );
    
        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                //label pointI = meshPoints[i];
                //Pout<< "   " << pointI << " moved from "
                //    << mesh_.points()[pointI] << " by "
                //    << mag(hitInfo[i].hitPoint()-mesh_.points()[pointI])
                //    << endl;
                newPoints[meshPoints[i]] = hitInfo[i].hitPoint();
            }
        }
    }

    // Per face (internal or coupled!) the patch that the
    // baffle should get (or -1).
    labelList facePatch(mesh_.nFaces(), -1);
    // Count of baffled faces
    label nBaffleFaces = 0;


//    // Sync position? Or not since same face on both side so just sync
//    // result of baffle.
//
//    const scalar minArea(readScalar(motionDict.lookup("minArea")));
//
//    Pout<< "markFacesOnProblemCellsGeometric : Comparing to minArea:"
//        << minArea << endl;
//
//    pointField facePoints;
//    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
//    {
//        const face& f = mesh_.faces()[faceI];
//
//        bool usesPatchPoint = false;
//
//        facePoints.setSize(f.size());
//        forAll(f, fp)
//        {
//            Map<label>::const_iterator iter = pp.meshPointMap().find(f[fp]);
//
//            if (iter != pp.meshPointMap().end())
//            {
//                facePoints[fp] = newPosition[iter()];
//                usesPatchPoint = true;
//            }
//            else
//            {
//                facePoints[fp] = mesh_.points()[f[fp]];
//            }
//        }
//
//        if (usesPatchPoint)
//        {
//            // Check area of face wrt original area
//            face identFace(identity(f.size()));
//
//            if (identFace.mag(facePoints) < minArea)
//            {
//                facePatch[faceI] = getBafflePatch(facePatch, faceI);
//                nBaffleFaces++;
//            }
//        }
//    }
//
//
//    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
//    forAll(patches, patchI)
//    {
//        const polyPatch& pp = patches[patchI];
//
//        if (pp.coupled())
//        {
//            forAll(pp, i)
//            {
//                label faceI = pp.start()+i;
//
//                const face& f = mesh_.faces()[faceI];
//
//                bool usesPatchPoint = false;
//
//                facePoints.setSize(f.size());
//                forAll(f, fp)
//                {
//                    Map<label>::const_iterator iter =
//                        pp.meshPointMap().find(f[fp]);
//
//                    if (iter != pp.meshPointMap().end())
//                    {
//                        facePoints[fp] = newPosition[iter()];
//                        usesPatchPoint = true;
//                    }
//                    else
//                    {
//                        facePoints[fp] = mesh_.points()[f[fp]];
//                    }
//                }
//
//                if (usesPatchPoint)
//                {
//                    // Check area of face wrt original area
//                    face identFace(identity(f.size()));
//
//                    if (identFace.mag(facePoints) < minArea)
//                    {
//                        facePatch[faceI] = getBafflePatch(facePatch, faceI);
//                        nBaffleFaces++;
//                    }
//                }
//            }
//        }
//    }

    {
        pointField oldPoints(mesh_.points());
        mesh_.movePoints(newPoints);
        faceSet wrongFaces(mesh_, "wrongFaces", 100);
        {
            //motionSmoother::checkMesh(false, mesh_, motionDict, wrongFaces);

            // Just check the errors from squashing
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            const labelList allFaces(identity(mesh_.nFaces()));
            label nWrongFaces = 0;

            scalar minArea(readScalar(motionDict.lookup("minArea")));
            if (minArea > -SMALL)
            {
                polyMeshGeometry::checkFaceArea
                (
                    false,
                    minArea,
                    mesh_,
                    mesh_.faceAreas(),
                    allFaces,
                    &wrongFaces
                );

                label nNewWrongFaces = returnReduce
                (
                    wrongFaces.size(),
                    sumOp<label>()
                );

                Info<< "    faces with area < "
                    << setw(5) << minArea
                    << " m^2                            : "
                    << nNewWrongFaces-nWrongFaces << endl;

                nWrongFaces = nNewWrongFaces;
            }

//            scalar minDet(readScalar(motionDict.lookup("minDeterminant")));
            scalar minDet = 0.01;
            if (minDet > -1)
            {
                polyMeshGeometry::checkCellDeterminant
                (
                    false,
                    minDet,
                    mesh_,
                    mesh_.faceAreas(),
                    allFaces,
                    polyMeshGeometry::affectedCells(mesh_, allFaces),
                    &wrongFaces
                );

                label nNewWrongFaces = returnReduce
                (
                    wrongFaces.size(),
                    sumOp<label>()
                );

                Info<< "    faces on cells with determinant < "
                    << setw(5) << minDet << "                : "
                    << nNewWrongFaces-nWrongFaces << endl;

                nWrongFaces = nNewWrongFaces;
            }
        }


        forAllConstIter(faceSet, wrongFaces, iter)
        {
            label patchI = mesh_.boundaryMesh().whichPatch(iter.key());

            if (patchI == -1 || mesh_.boundaryMesh()[patchI].coupled())
            {
                facePatch[iter.key()] = getBafflePatch(facePatch, iter.key());
                nBaffleFaces++;

                //Pout<< "    " << iter.key()
                //    //<< " on patch " << mesh_.boundaryMesh()[patchI].name()
                //    << " is destined for patch " << facePatch[iter.key()]
                //    << endl;
            }
        }
        // Restore points.
        mesh_.movePoints(oldPoints);
    }


    Info<< "markFacesOnProblemCellsGeometric : marked "
        << returnReduce(nBaffleFaces, sumOp<label>())
        << " additional internal and coupled faces"
        << " to be converted into baffles." << endl;

    syncTools::syncFaceList
    (
        mesh_,
        facePatch,
        maxEqOp<label>()
    );

    return facePatch;
}
//XXXXXXXX


// Return a list of coupled face pairs, i.e. faces that use the same vertices.
// (this information is recalculated instead of maintained since would be too
// hard across splitMeshRegions).
Foam::List<Foam::labelPair> Foam::meshRefinement::getDuplicateFaces
(
    const labelList& testFaces
) const
{
    labelList duplicateFace
    (
        localPointRegion::findDuplicateFaces
        (
            mesh_,
            testFaces
        )
    );

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Convert into list of coupled face pairs (mesh face labels).
    List<labelPair> duplicateFaces(testFaces.size());
    label dupI = 0;

    forAll(duplicateFace, i)
    {
        label otherFaceI = duplicateFace[i];

        if (otherFaceI != -1 && i < otherFaceI)
        {
            label meshFace0 = testFaces[i];
            label patch0 = patches.whichPatch(meshFace0);
            label meshFace1 = testFaces[otherFaceI];
            label patch1 = patches.whichPatch(meshFace1);

            if
            (
                (patch0 != -1 && isA<processorPolyPatch>(patches[patch0]))
             || (patch1 != -1 && isA<processorPolyPatch>(patches[patch1]))
            )
            {
                FatalErrorIn
                (
                    "meshRefinement::getDuplicateFaces"
                    "(const bool, const labelList&)"
                )   << "One of two duplicate faces is on"
                    << " processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << meshFace0
                    << " is on patch:" << patches[patch0].name()
                    << nl
                    << "Face:" << meshFace1
                    << " is on patch:" << patches[patch1].name()
                    << abort(FatalError);
            }

            duplicateFaces[dupI++] = labelPair(meshFace0, meshFace1);
        }
    }
    duplicateFaces.setSize(dupI);

    Info<< "getDuplicateFaces : found " << returnReduce(dupI, sumOp<label>())
        << " pairs of duplicate faces." << nl << endl;


    if (debug)
    {
        faceSet duplicateFaceSet(mesh_, "duplicateFaces", 2*dupI);

        forAll(duplicateFaces, i)
        {
            duplicateFaceSet.insert(duplicateFaces[i][0]);
            duplicateFaceSet.insert(duplicateFaces[i][1]);
        }
        Pout<< "Writing duplicate faces (baffles) to faceSet "
            << duplicateFaceSet.name() << nl << endl;
        duplicateFaceSet.write();
    }

    return duplicateFaces;
}


// Extract those baffles (duplicate) faces that are on the edge of a baffle
// region. These are candidates for merging.
// Done by counting the number of baffles faces per mesh edge. If edge
// has 2 boundary faces and both are baffle faces it is the edge of a baffle
// region.
Foam::List<Foam::labelPair> Foam::meshRefinement::filterDuplicateFaces
(
    const List<labelPair>& couples
) const
{
    // Construct addressing engine for all duplicate faces (only one
    // for each pair)

    // Duplicate faces in mesh labels (first face of each pair only)
    // (reused later on to mark off filtered couples. see below)
    labelList duplicateFaces(couples.size());


    // All duplicate faces on edge of the patch are to be merged.
    // So we count for all edges of duplicate faces how many duplicate
    // faces use them.
    labelList nBafflesPerEdge(mesh_.nEdges(), 0);



    // Count number of boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Count number of boundary faces. Discard coupled boundary faces.
        if (!pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                const labelList& fEdges = mesh_.faceEdges(faceI);

                forAll(fEdges, fEdgeI)
                {
                    nBafflesPerEdge[fEdges[fEdgeI]]++;
                }
                faceI++;
            }
        }
    }


    DynamicList<label> fe0;
    DynamicList<label> fe1;


    // Count number of duplicate boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(couples, i)
    {
        const labelList& fEdges0 = mesh_.faceEdges(couples[i].first(), fe0);

        forAll(fEdges0, fEdgeI)
        {
            nBafflesPerEdge[fEdges0[fEdgeI]] += 1000000;
        }

        const labelList& fEdges1 = mesh_.faceEdges(couples[i].second(), fe1);

        forAll(fEdges1, fEdgeI)
        {
            nBafflesPerEdge[fEdges1[fEdgeI]] += 1000000;
        }
    }

    // Add nBaffles on shared edges
    syncTools::syncEdgeList
    (
        mesh_,
        nBafflesPerEdge,
        plusEqOp<label>(),  // in-place add
        0                   // initial value
    );


    // Baffles which are not next to other boundaries and baffles will have
    // value 2*1000000+2*1

    List<labelPair> filteredCouples(couples.size());
    label filterI = 0;

    forAll(couples, i)
    {
        const labelPair& couple = couples[i];

        if
        (
            patches.whichPatch(couple.first())
         == patches.whichPatch(couple.second())
        )
        {
            const labelList& fEdges = mesh_.faceEdges(couples[i].first());

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (nBafflesPerEdge[edgeI] == 2*1000000+2*1)
                {
                    filteredCouples[filterI++] = couples[i];
                    break;
                }
            }
        }
    }
    filteredCouples.setSize(filterI);

    //Info<< "filterDuplicateFaces : from "
    //    << returnReduce(couples.size(), sumOp<label>())
    //    << " down to "
    //    << returnReduce(filteredCouples.size(), sumOp<label>())
    //    << " baffles." << nl << endl;

    return filteredCouples;
}


// Merge baffles. Gets pairs of faces.
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeBaffles
(
    const List<labelPair>& couples
)
{
    // Mesh change engine
    polyTopoChange meshMod(mesh_);

    const faceList& faces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();
    const faceZoneMesh& faceZones = mesh_.faceZones();

    forAll(couples, i)
    {
        label face0 = couples[i].first();
        label face1 = couples[i].second();

        // face1 < 0 signals a coupled face that has been converted to baffle.

        label own0 = faceOwner[face0];
        label own1 = faceOwner[face1];

        if (face1 < 0 || own0 < own1)
        {
            // Use face0 as the new internal face.
            label zoneID = faceZones.whichZone(face0);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
            }

            label nei = (face1 < 0 ? -1 : own1);

            meshMod.setAction(polyRemoveFace(face1));
            meshMod.setAction
            (
                polyModifyFace
                (
                    faces[face0],           // modified face
                    face0,                  // label of face being modified
                    own0,                   // owner
                    nei,                    // neighbour
                    false,                  // face flip
                    -1,                     // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip,               // face flip in zone
                    -1                      // sub patch
                )
            );
        }
        else
        {
            // Use face1 as the new internal face.
            label zoneID = faceZones.whichZone(face1);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
            }

            meshMod.setAction(polyRemoveFace(face0));
            meshMod.setAction
            (
                polyModifyFace
                (
                    faces[face1],           // modified face
                    face1,                  // label of face being modified
                    own1,                   // owner
                    own0,                   // neighbour
                    false,                  // face flip
                    -1,                     // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip,               // face flip in zone
                    -1                      // sub patch
                )
            );
        }
    }

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    // Update intersections. Recalculate intersections on merged faces since
    // this seems to give problems? Note: should not be nessecary since
    // baffles preserve intersections from when they were created.
    labelList newExposedFaces(2*couples.size());
    label newI = 0;

    forAll(couples, i)
    {
        label newFace0 = map().reverseFaceMap()[couples[i].first()];
        if (newFace0 != -1)
        {
            newExposedFaces[newI++] = newFace0;
        }

        label newFace1 = map().reverseFaceMap()[couples[i].second()];
        if (newFace1 != -1)
        {
            newExposedFaces[newI++] = newFace1;
        }
    }
    newExposedFaces.setSize(newI);
    updateMesh(map, newExposedFaces);

    return map;
}


// Finds region per cell for cells inside closed named surfaces
void Foam::meshRefinement::findCellZoneGeometric
(
    const labelList& closedNamedSurfaces,   // indices of closed surfaces
    const labelList& namedSurfaceIndex,     // per face index of named surface
    const labelList& surfaceToCellZone,     // cell zone index per surface

    labelList& cellToZone
) const
{
    const pointField& cellCentres = mesh_.cellCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Check if cell centre is inside
    labelList insideSurfaces;
    surfaces_.findInside
    (
        closedNamedSurfaces,
        cellCentres,
        insideSurfaces
    );

    forAll(insideSurfaces, cellI)
    {
        if (cellToZone[cellI] == -2)
        {
            label surfI = insideSurfaces[cellI];

            if (surfI != -1)
            {
                cellToZone[cellI] = surfaceToCellZone[surfI];
            }
        }
    }


    // Some cells with cell centres close to surface might have
    // had been put into wrong surface. Recheck with perturbed cell centre.


    // 1. Collect points

    // Count points to test.
    label nCandidates = 0;
    forAll(namedSurfaceIndex, faceI)
    {
        label surfI = namedSurfaceIndex[faceI];

        if (surfI != -1)
        {
            if (mesh_.isInternalFace(faceI))
            {
                nCandidates += 2;
            }
            else
            {
                nCandidates += 1;
            }
        }
    }

    // Collect points.
    pointField candidatePoints(nCandidates);
    nCandidates = 0;
    forAll(namedSurfaceIndex, faceI)
    {
        label surfI = namedSurfaceIndex[faceI];

        if (surfI != -1)
        {
            label own = faceOwner[faceI];
            const point& ownCc = cellCentres[own];

            if (mesh_.isInternalFace(faceI))
            {
                label nei = faceNeighbour[faceI];
                const point& neiCc = cellCentres[nei];
                // Perturbed cc
                const vector d = 1E-4*(neiCc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
                candidatePoints[nCandidates++] = neiCc+d;
            }
            else
            {
                const point& neiFc = mesh_.faceCentres()[faceI];
                // Perturbed cc
                const vector d = 1E-4*(neiFc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
            }
        }
    }


    // 2. Test points for inside

    surfaces_.findInside
    (
        closedNamedSurfaces,
        candidatePoints,
        insideSurfaces
    );


    // 3. Update zone information

    nCandidates = 0;
    forAll(namedSurfaceIndex, faceI)
    {
        label surfI = namedSurfaceIndex[faceI];

        if (surfI != -1)
        {
            label own = faceOwner[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }

                label neiSurfI = insideSurfaces[nCandidates++];
                if (neiSurfI != -1)
                {
                    label nei = faceNeighbour[faceI];

                    cellToZone[nei] = surfaceToCellZone[neiSurfI];
                }
            }
            else
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }
            }
        }
    }
}


// Finds region per cell. Assumes:
// - region containing keepPoint does not go into a cellZone
// - all other regions can be found by crossing faces marked in
//   namedSurfaceIndex.
void Foam::meshRefinement::findCellZoneTopo
(
    const point& keepPoint,
    const labelList& namedSurfaceIndex,
    const labelList& surfaceToCellZone,
    labelList& cellToZone
) const
{
    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());

    forAll(namedSurfaceIndex, faceI)
    {
        if (namedSurfaceIndex[faceI] == -1)
        {
            blockedFace[faceI] = false;
        }
        else
        {
            blockedFace[faceI] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Per mesh region the zone the cell should be put in.
    // -2   : not analysed yet
    // -1   : keepPoint region. Not put into any cellzone.
    // >= 0 : index of cellZone
    labelList regionToCellZone(cellRegion.nRegions(), -2);

    // See which cells already are set in the cellToZone (from geometric
    // searching) and use these to take over their zones.
    // Note: could be improved to count number of cells per region.
    forAll(cellToZone, cellI)
    {
        if (cellToZone[cellI] != -2)
        {
            regionToCellZone[cellRegion[cellI]] = cellToZone[cellI];
        }
    }



    // Find the region containing the keepPoint
    label keepRegionI = -1;

    label cellI = mesh_.findCell(keepPoint);

    if (cellI != -1)
    {
        keepRegionI = cellRegion[cellI];
    }
    reduce(keepRegionI, maxOp<label>());

    Info<< "Found point " << keepPoint << " in cell " << cellI
        << " in global region " << keepRegionI
        << " out of " << cellRegion.nRegions() << " regions." << endl;

    if (keepRegionI == -1)
    {
        FatalErrorIn
        (
            "meshRefinement::findCellZoneTopo"
            "(const point&, const labelList&, const labelList&, labelList&)"
        )   << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh_.globalData().bb()
            << exit(FatalError);
    }

    // Mark default region with zone -1.
    if (regionToCellZone[keepRegionI] == -2)
    {
        regionToCellZone[keepRegionI] = -1;
    }


    // Find correspondence between cell zone and surface
    // by changing cell zone every time we cross a surface.
    while (true)
    {
        bool changed = false;

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label surfI = namedSurfaceIndex[faceI];

            if (surfI != -1)
            {
                // Get cell zone that surface cells are in
                label surfZoneI = surfaceToCellZone[surfI];

                // Check whether inbetween different regions
                label ownRegion = cellRegion[mesh_.faceOwner()[faceI]];
                label neiRegion = cellRegion[mesh_.faceNeighbour()[faceI]];

                if (ownRegion != neiRegion)
                {
                    // Jump. Change one of the sides to my type.

                    // 1. Interface between my type and unset region.
                    // Set region to keepRegion

                    if (regionToCellZone[ownRegion] == -2)
                    {
                        if (regionToCellZone[neiRegion] == surfZoneI)
                        {
                            // Face between unset and my region. Put unset
                            // region into keepRegion
                            regionToCellZone[ownRegion] = -1;
                            changed = true;
                        }
                        else if (regionToCellZone[neiRegion] != -2)
                        {
                            // Face between unset and other region.
                            // Put unset region into my region
                            regionToCellZone[ownRegion] = surfZoneI;
                            changed = true;
                        }
                    }
                    else if (regionToCellZone[neiRegion] == -2)
                    {
                        if (regionToCellZone[ownRegion] == surfZoneI)
                        {
                            // Face between unset and my region. Put unset
                            // region into keepRegion
                            regionToCellZone[neiRegion] = -1;
                            changed = true;
                        }
                        else if (regionToCellZone[ownRegion] != -2)
                        {
                            // Face between unset and other region.
                            // Put unset region into my region
                            regionToCellZone[neiRegion] = surfZoneI;
                            changed = true;
                        }
                    }
                }
            }
        }

        if (!returnReduce(changed, orOp<bool>()))
        {
            break;
        }
    }


    // Synchronise regionToCellZone.
    // Note:
    // - region numbers are identical on all processors
    // - keepRegion is identical ,,
    // - cellZones are identical ,,
    Pstream::listCombineGather(regionToCellZone, maxEqOp<label>());
    Pstream::listCombineScatter(regionToCellZone);


    forAll(regionToCellZone, regionI)
    {
        label zoneI = regionToCellZone[regionI];

        if (zoneI ==  -2)
        {
            FatalErrorIn
            (
                "meshRefinement::findCellZoneTopo"
                "(const point&, const labelList&, const labelList&, labelList&)"
            )   << "For region " << regionI << " haven't set cell zone."
                << exit(FatalError);
        }
    }

    if (debug)
    {
        forAll(regionToCellZone, regionI)
        {
            Pout<< "Region " << regionI
                << " becomes cellZone:" << regionToCellZone[regionI]
                << endl;
        }
    }

    // Rework into cellToZone
    forAll(cellToZone, cellI)
    {
        cellToZone[cellI] = regionToCellZone[cellRegion[cellI]];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Split off unreachable areas of mesh.
void Foam::meshRefinement::baffleAndSplitMesh
(
    const bool handleSnapProblems,
    const bool mergeFreeStanding,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToPatch,
    const point& keepPoint
)
{
    // Introduce baffles
    // ~~~~~~~~~~~~~~~~~

    // Split the mesh along internal faces wherever there is a pierce between
    // two cell centres.

    Info<< "Introducing baffles for "
        << returnReduce(countHits(), sumOp<label>())
        << " faces that are intersected by the surface." << nl << endl;

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList ownPatch, neiPatch;
    getBafflePatches
    (
        globalToPatch,
        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );

    if (debug)
    {
        runTime++;
    }

    createBaffles(ownPatch, neiPatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }

    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles");

    if (debug)
    {
        Pout<< "Writing baffled mesh to time " << mesh_.time().timeName()
            << endl;
        write(debug, runTime.path()/"baffles");
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }


    // Introduce baffles to delete problem cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Create some additional baffles where we want surface cells removed.

    if (handleSnapProblems)
    {
        Info<< nl
            << "Introducing baffles to block off problem cells" << nl
            << "----------------------------------------------" << nl
            << endl;

        labelList facePatch
        (
            markFacesOnProblemCells
            (
                globalToPatch
            )
            //markFacesOnProblemCellsGeometric
            //(
            //    motionDict,
            //    globalToPatch
            //)
        );
        Info<< "Analyzed problem cells in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

        if (debug)
        {
            // Dump all these faces to a faceSet.
            faceSet problemGeom(mesh_, "problemFacesGeom", 100);

            const labelList facePatchGeom
            (
                markFacesOnProblemCellsGeometric
                (
                    motionDict,
                    globalToPatch
                )
            );
            forAll(facePatchGeom, faceI)
            {
                if (facePatchGeom[faceI] != -1)
                {
                    problemGeom.insert(faceI);
                }
            }
            Pout<< "Dumping " << problemGeom.size()
                << " problem faces to " << problemGeom.objectPath() << endl;
            problemGeom.write();


            faceSet problemTopo(mesh_, "problemFacesTopo", 100);

            const labelList facePatchTopo
            (    
                markFacesOnProblemCells
                (
                    globalToPatch
                )
            );
            forAll(facePatchTopo, faceI)
            {
                if (facePatchTopo[faceI] != -1)
                {
                    problemTopo.insert(faceI);
                }
            }
            Pout<< "Dumping " << problemTopo.size()
                << " problem faces to " << problemTopo.objectPath() << endl;
            problemTopo.write();
        }

        Info<< "Introducing baffles to delete problem cells." << nl << endl;

        if (debug)
        {
            runTime++;
        }

        // Create baffles with same owner and neighbour for now.
        createBaffles(facePatch, facePatch);

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            checkData();
        }
        Info<< "Created baffles in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

        printMeshInfo(debug, "After introducing baffles");

        if (debug)
        {
            Pout<< "Writing extra baffled mesh to time "
                << mesh_.time().timeName() << endl;
            write(debug, runTime.path()/"extraBaffles");
            Pout<< "Dumped debug data in = "
                << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
        }
    }


    // Select part of mesh
    // ~~~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Remove unreachable sections of mesh" << nl
        << "-----------------------------------" << nl
        << endl;

    if (debug)
    {
        runTime++;
    }

    splitMeshRegions(keepPoint);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Split mesh in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After subsetting");

    if (debug)
    {
        Pout<< "Writing subsetted mesh to time " << mesh_.time().timeName()
            << endl;
        write(debug, runTime.path()/runTime.timeName());
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }


    // Merge baffles
    // ~~~~~~~~~~~~~

    if (mergeFreeStanding)
    {
        Info<< nl
            << "Merge free-standing baffles" << nl
            << "---------------------------" << nl
            << endl;

        if (debug)
        {
            runTime++;
        }

        // List of pairs of freestanding baffle faces.
        List<labelPair> couples
        (
            filterDuplicateFaces    // filter out freestanding baffles
            (
                getDuplicateFaces   // get all baffles
                (
                    identity(mesh_.nFaces()-mesh_.nInternalFaces())
                   +mesh_.nInternalFaces()
                )
            )
        );

        label nCouples = couples.size();
        reduce(nCouples, sumOp<label>());

        Info<< "Detected free-standing baffles : " << nCouples << endl;

        if (nCouples > 0)
        {
            // Actually merge baffles. Note: not exactly parallellized. Should
            // convert baffle faces into processor faces if they resulted
            // from them.
            mergeBaffles(couples);

            if (debug)
            {
                // Debug:test all is still synced across proc patches
                checkData();
            }
        }
        Info<< "Merged free-standing baffles in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


// Split off (with optional buffer layers) unreachable areas of mesh.
Foam::autoPtr<Foam::mapPolyMesh>  Foam::meshRefinement::splitMesh
(
    const label nBufferLayers,
    const labelList& globalToPatch,
    const point& keepPoint
)
{
    // Determine patches to put intersections into
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList ownPatch, neiPatch;
    getBafflePatches
    (
        globalToPatch,
        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces(), false);

    forAll(ownPatch, faceI)
    {
        if (ownPatch[faceI] != -1 || neiPatch[faceI] != -1)
        {
            blockedFace[faceI] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Find the region containing the keepPoint
    label keepRegionI = -1;

    label cellI = mesh_.findCell(keepPoint);

    if (cellI != -1)
    {
        keepRegionI = cellRegion[cellI];
    }
    reduce(keepRegionI, maxOp<label>());

    Info<< "Found point " << keepPoint << " in cell " << cellI
        << " in global region " << keepRegionI
        << " out of " << cellRegion.nRegions() << " regions." << endl;

    if (keepRegionI == -1)
    {
        FatalErrorIn
        (
            "meshRefinement::findCellZoneTopo"
            "(const point&, const labelList&, const labelList&, labelList&)"
        )   << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh_.globalData().bb()
            << exit(FatalError);
    }


    // Walk out nBufferlayers from region split
    // (modifies cellRegion, ownPatch)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Takes over face patch onto points and then back to faces and cells
    // (so cell-face-point walk)

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Patch for exposed faces for lack of anything sensible.
    label defaultPatch = 0;
    if (globalToPatch.size() > 0)
    {
        defaultPatch = globalToPatch[0];
    }

    for (label i = 0; i < nBufferLayers; i++)
    {
        // 1. From cells (via faces) to points

        labelList pointBaffle(mesh_.nPoints(), -1);

        forAll(faceNeighbour, faceI)
        {
            const face& f = mesh_.faces()[faceI];

            label ownRegion = cellRegion[faceOwner[faceI]];
            label neiRegion = cellRegion[faceNeighbour[faceI]];

            if (ownRegion == keepRegionI && neiRegion != keepRegionI)
            {
                // Note max(..) since possibly regionSplit might have split
                // off extra unreachable parts of mesh. Note: or can this only
                // happen for boundary faces?
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[faceI]);
                }
            }
            else if (ownRegion != keepRegionI && neiRegion == keepRegionI)
            {
                label newPatchI = neiPatch[faceI];
                if (newPatchI == -1)
                {
                    newPatchI = max(defaultPatch, ownPatch[faceI]);
                }
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = newPatchI;
                }
            }
        }
        for
        (
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            const face& f = mesh_.faces()[faceI];

            label ownRegion = cellRegion[faceOwner[faceI]];

            if (ownRegion == keepRegionI)
            {
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[faceI]);
                }
            }
        }

        // Sync
        syncTools::syncPointList
        (
            mesh_,
            pointBaffle,
            maxEqOp<label>(),
            -1                  // null value
        );


        // 2. From points back to faces

        const labelListList& pointFaces = mesh_.pointFaces();

        forAll(pointFaces, pointI)
        {
            if (pointBaffle[pointI] != -1)
            {
                const labelList& pFaces = pointFaces[pointI];

                forAll(pFaces, pFaceI)
                {
                    label faceI = pFaces[pFaceI];

                    if (ownPatch[faceI] == -1)
                    {
                        ownPatch[faceI] = pointBaffle[pointI];
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());


        // 3. From faces to cells (cellRegion) and back to faces (ownPatch)

        labelList newOwnPatch(ownPatch);

        forAll(ownPatch, faceI)
        {
            if (ownPatch[faceI] != -1)
            {
                label own = faceOwner[faceI];

                if (cellRegion[own] != keepRegionI)
                {
                    cellRegion[own] = keepRegionI;

                    const cell& ownFaces = mesh_.cells()[own];
                    forAll(ownFaces, j)
                    {
                        if (ownPatch[ownFaces[j]] == -1)
                        {
                            newOwnPatch[ownFaces[j]] = ownPatch[faceI];
                        }
                    }
                }
                if (mesh_.isInternalFace(faceI))
                {
                    label nei = faceNeighbour[faceI];

                    if (cellRegion[nei] != keepRegionI)
                    {
                        cellRegion[nei] = keepRegionI;

                        const cell& neiFaces = mesh_.cells()[nei];
                        forAll(neiFaces, j)
                        {
                            if (ownPatch[neiFaces[j]] == -1)
                            {
                                newOwnPatch[neiFaces[j]] = ownPatch[faceI];
                            }
                        }
                    }
                }
            }
        }

        ownPatch.transfer(newOwnPatch);

        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    }



    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, cellI)
    {
        if (cellRegion[cellI] != keepRegionI)
        {
            cellsToRemove.append(cellI);
        }
    }
    cellsToRemove.shrink();

    label nCellsToKeep = mesh_.nCells() - cellsToRemove.size();
    reduce(nCellsToKeep, sumOp<label>());

    Info<< "Keeping all cells in region " << keepRegionI
        << " containing point " << keepPoint << endl
        << "Selected for keeping : " << nCellsToKeep
        << " cells." << endl;


    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatches(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label faceI = exposedFaces[i];

        if (ownPatch[faceI] != -1)
        {
            exposedPatches[i] = ownPatch[faceI];
        }
        else
        {
            WarningIn("meshRefinement::splitMesh(..)")
                << "For exposed face " << faceI
                << " fc:" << mesh_.faceCentres()[faceI]
                << " found no patch." << endl
                << "    Taking patch " << defaultPatch
                << " instead." << endl;
            exposedPatches[i] = defaultPatch;
        }
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


// Find boundary points that connect to more than one cell region and
// split them.
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints()
{
    // Topochange container
    polyTopoChange meshMod(mesh_);


    // Analyse which points need to be duplicated
    localPointRegion regionSide(mesh_);

    label nNonManifPoints = returnReduce
    (
        regionSide.meshPointMap().size(),
        sumOp<label>()
    );

    Info<< "dupNonManifoldPoints : Found : " << nNonManifPoints
        << " non-manifold points (out of "
        << mesh_.globalData().nTotalPoints()
        << ')' << endl;

    // Topo change engine
    duplicatePoints pointDuplicator(mesh_);

    // Insert changes into meshMod
    pointDuplicator.setRefinement(regionSide, meshMod);

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    // Update intersections. Is mapping only (no faces created, positions stay
    // same) so no need to recalculate intersections.
    updateMesh(map, labelList(0));

    return map;
}


// Zoning
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::zonify
(
    const point& keepPoint
)
{
    const wordList& cellZoneNames = surfaces_.cellZoneNames();
    const wordList& faceZoneNames = surfaces_.faceZoneNames();

    labelList namedSurfaces(surfaces_.getNamedSurfaces());

    boolList isNamedSurface(cellZoneNames.size(), false);

    forAll(namedSurfaces, i)
    {
        label surfI = namedSurfaces[i];

        isNamedSurface[surfI] = true;

        Info<< "Surface : " << surfaces_.names()[surfI] << nl
            << "    faceZone : " << faceZoneNames[surfI] << nl
            << "    cellZone : " << cellZoneNames[surfI] << endl;
    }


    // Add zones to mesh

    labelList surfaceToFaceZone(faceZoneNames.size(), -1);
    {
        faceZoneMesh& faceZones = mesh_.faceZones();

        forAll(namedSurfaces, i)
        {
            label surfI = namedSurfaces[i];

            label zoneI = faceZones.findZoneID(faceZoneNames[surfI]);

            if (zoneI == -1)
            {
                zoneI = faceZones.size();
                faceZones.setSize(zoneI+1);
                faceZones.set
                (
                    zoneI,
                    new faceZone
                    (
                        faceZoneNames[surfI],   //name
                        labelList(0),           //addressing
                        boolList(0),            //flipmap
                        zoneI,                  //index
                        faceZones               //faceZoneMesh
                    )
                );
            }

            if (debug)
            {
                Pout<< "Faces on " << surfaces_.names()[surfI]
                    << " will go into faceZone " << zoneI << endl;
            }
            surfaceToFaceZone[surfI] = zoneI;
        }

        // Check they are synced
        List<wordList> allFaceZones(Pstream::nProcs());
        allFaceZones[Pstream::myProcNo()] = faceZones.names();
        Pstream::gatherList(allFaceZones);
        Pstream::scatterList(allFaceZones);

        for (label procI = 1; procI < allFaceZones.size(); procI++)
        {
            if (allFaceZones[procI] != allFaceZones[0])
            {
                FatalErrorIn
                (
                    "meshRefinement::zonify"
                    "(const label, const point&)"
                )   << "Zones not synchronised among processors." << nl
                    << " Processor0 has faceZones:" << allFaceZones[0]
                    << " , processor" << procI
                    << " has faceZones:" << allFaceZones[procI]
                    << exit(FatalError);
            }
        }
    }

    labelList surfaceToCellZone(cellZoneNames.size(), -1);
    {
        cellZoneMesh& cellZones = mesh_.cellZones();

        forAll(namedSurfaces, i)
        {
            label surfI = namedSurfaces[i];

            label zoneI = cellZones.findZoneID(cellZoneNames[surfI]);

            if (zoneI == -1)
            {
                zoneI = cellZones.size();
                cellZones.setSize(zoneI+1);
                cellZones.set
                (
                    zoneI,
                    new cellZone
                    (
                        cellZoneNames[surfI],   //name
                        labelList(0),           //addressing
                        zoneI,                  //index
                        cellZones               //cellZoneMesh
                    )
                );
            }

            if (debug)
            {
                Pout<< "Cells inside " << surfaces_.names()[surfI]
                    << " will go into cellZone " << zoneI << endl;
            }
            surfaceToCellZone[surfI] = zoneI;
        }

        // Check they are synced
        List<wordList> allCellZones(Pstream::nProcs());
        allCellZones[Pstream::myProcNo()] = cellZones.names();
        Pstream::gatherList(allCellZones);
        Pstream::scatterList(allCellZones);

        for (label procI = 1; procI < allCellZones.size(); procI++)
        {
            if (allCellZones[procI] != allCellZones[0])
            {
                FatalErrorIn
                (
                    "meshRefinement::zonify"
                    "(const label, const point&)"
                )   << "Zones not synchronised among processors." << nl
                    << " Processor0 has cellZones:" << allCellZones[0]
                    << " , processor" << procI
                    << " has cellZones:" << allCellZones[procI]
                    << exit(FatalError);
            }
        }
    }



    const pointField& cellCentres = mesh_.cellCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    // Mark faces intersecting zoned surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // Like surfaceIndex_ but only for named surfaces.
    labelList namedSurfaceIndex(mesh_.nFaces(), -1);

    {
        // Statistics: number of faces per faceZone
        labelList nSurfFaces(faceZoneNames.size(), 0);

        // Swap neighbouring cell centres and cell level
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
        calcNeighbourData(neiLevel, neiCc);

        // Note: for all internal faces? internal + coupled?
        // Since zonify is run after baffling the surfaceIndex_ on baffles is
        // not synchronised across both baffle faces. Fortunately we don't
        // do zonify baffle faces anyway (they are normal boundary faces).

        // Collect candidate faces
        // ~~~~~~~~~~~~~~~~~~~~~~~

        labelList testFaces(intersectedFaces());

        // Collect segments
        // ~~~~~~~~~~~~~~~~

        pointField start(testFaces.size());
        pointField end(testFaces.size());

        forAll(testFaces, i)
        {
            label faceI = testFaces[i];

            label own = mesh_.faceOwner()[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                start[i] = cellCentres[own];
                end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
            }
            else
            {
                start[i] = cellCentres[own];
                end[i] = neiCc[faceI-mesh_.nInternalFaces()];
            }
        }


        // Do test for intersections
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Note that we intersect all intersected faces again. Could reuse
        // the information already in surfaceIndex_.

        labelList surface1;
        labelList surface2;
        {
            List<pointIndexHit> hit1;
            labelList region1;
            List<pointIndexHit> hit2;
            labelList region2;
            surfaces_.findNearestIntersection
            (
                namedSurfaces,
                start,
                end,

                surface1,
                hit1,
                region1,
                surface2,
                hit2,
                region2
            );
        }

        forAll(testFaces, i)
        {
            label faceI = testFaces[i];

            if (surface1[i] != -1)
            {
                // If both hit should probably choose nearest. For later.
                namedSurfaceIndex[faceI] = surface1[i];
                nSurfFaces[surface1[i]]++;
            }
            else if (surface2[i] != -1)
            {
                namedSurfaceIndex[faceI] = surface2[i];
                nSurfFaces[surface2[i]]++;
            }
        }


        // surfaceIndex migh have different surfaces on both sides if
        // there happen to be a (obviously thin) surface with different
        // regions between the cell centres. If one is on a named surface
        // and the other is not this might give problems so sync.
        syncTools::syncFaceList
        (
            mesh_,
            namedSurfaceIndex,
            maxEqOp<label>()
        );

        // Print a bit
        if (debug)
        {
            forAll(nSurfFaces, surfI)
            {
                Pout<< "Surface:"
                    << surfaces_.names()[surfI]
                    << "  nZoneFaces:" << nSurfFaces[surfI] << nl;
            }
            Pout<< endl;
        }
    }


    // Topochange container
    polyTopoChange meshMod(mesh_);


    // Put the faces into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label surfI = namedSurfaceIndex[faceI];

        if (surfI != -1)
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[faceI],           // modified face
                    faceI,                          // label of face
                    faceOwner[faceI],               // owner
                    faceNeighbour[faceI],           // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    surfaceToFaceZone[surfI],       // zone for face
                    false,                          // face flip in zone
                    -1                              // sub patch
                )
            );
        }
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        label faceI = pp.start();

        forAll(pp, i)
        {
            label surfI = namedSurfaceIndex[faceI];

            if (surfI != -1)
            {
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        mesh_.faces()[faceI],           // modified face
                        faceI,                          // label of face
                        faceOwner[faceI],               // owner
                        -1,                             // neighbour
                        false,                          // face flip
                        patchI,                         // patch for face
                        false,                          // remove from zone
                        surfaceToFaceZone[surfI],       // zone for face
                        false,                          // face flip in zone
                        polyTopoChange::whichSubPatch(pp, faceI)    // sub patch
                    )
                );
            }
            faceI++;
        }
    }


    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Closed surfaces with cellZone specified.
    labelList closedNamedSurfaces(surfaces_.getClosedNamedSurfaces());

    // Zone per cell:
    // -2 : unset
    // -1 : not in any zone
    // >=0: zoneID
    labelList cellToZone(mesh_.nCells(), -2);


    // Set using geometric test
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    if (closedNamedSurfaces.size() > 0)
    {
        findCellZoneGeometric
        (
            closedNamedSurfaces,   // indices of closed surfaces
            namedSurfaceIndex,     // per face index of named surface
            surfaceToCellZone,     // cell zone index per surface
            cellToZone
        );
    }


    // Set using walking
    // ~~~~~~~~~~~~~~~~~

    //if (returnReduce(nSet, sumOp<label>()) < mesh_.globalData().nTotalCells())
    {
        // Topological walk
        findCellZoneTopo
        (
            keepPoint,
            namedSurfaceIndex,
            surfaceToCellZone,
            cellToZone
        );
    }


    // Actually move the cells to their zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(cellToZone, cellI)
    {
        label zoneI = cellToZone[cellI];

        if (zoneI >= 0)
        {
            meshMod.setAction
            (
                polyModifyCell
                (
                    cellI,
                    false,          // removeFromZone
                    zoneI
                )
            );
        }
    }

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    return map;
}


// ************************************************************************* //
