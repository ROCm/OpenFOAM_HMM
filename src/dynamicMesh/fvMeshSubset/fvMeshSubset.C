/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "fvMeshSubset.H"
#include "boolList.H"
#include "BitOps.H"
#include "pointIndList.H"
#include "Pstream.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "removeCells.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Mark faces/points (with 0) in labelList
inline void markUsed
(
    const Foam::labelUList& locations,
    Foam::labelList& map
)
{
    for (auto idx : locations)
    {
        map[idx] = 0;
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fvMeshSubset::exposedPatchName("oldInternalFaces");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fvMeshSubset::checkCellSubset() const
{
    if (!fvMeshSubsetPtr_)
    {
        FatalErrorInFunction
            << "setCellSubset()" << nl
            << "before attempting to access subset data"
            << abort(FatalError);

        return false;
    }

    return true;
}


void Foam::fvMeshSubset::calcFaceFlipMap() const
{
    const labelList& subToBaseFace = faceMap();
    const labelList& subToBaseCell = cellMap();

    faceFlipMapPtr_.reset(new labelList(subToBaseFace.size()));
    auto& faceFlipMap = *faceFlipMapPtr_;

    // Only exposed internal faces might be flipped (since we don't do
    // any cell renumbering, just compacting)
    const label subInt = subMesh().nInternalFaces();

    const labelList& subOwn = subMesh().faceOwner();
    const labelList& own = baseMesh_.faceOwner();

    for (label subFaceI = 0; subFaceI < subInt; ++subFaceI)
    {
        faceFlipMap[subFaceI] = subToBaseFace[subFaceI]+1;
    }
    for (label subFaceI = subInt; subFaceI < subOwn.size(); ++subFaceI)
    {
        const label faceI = subToBaseFace[subFaceI];
        if (subToBaseCell[subOwn[subFaceI]] == own[faceI])
        {
            faceFlipMap[subFaceI] = faceI+1;
        }
        else
        {
            faceFlipMap[subFaceI] = -faceI-1;
        }
    }
}


void Foam::fvMeshSubset::doCoupledPatches
(
    const bool syncPar,
    labelList& nCellsUsingFace
) const
{
    // Synchronize facesToSubset on both sides of coupled patches.
    // Marks faces that become 'uncoupled' with 3.

    const polyBoundaryMesh& oldPatches = baseMesh().boundaryMesh();

    label nUncoupled = 0;

    if (syncPar && Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send face usage across processor patches
        for (const polyPatch& pp : oldPatches)
        {
            const auto* procPatch = isA<processorPolyPatch>(pp);

            if (procPatch)
            {
                const label nbrProci = procPatch->neighbProcNo();

                UOPstream toNeighbour(nbrProci, pBufs);

                if (!nCellsUsingFace.empty())
                {
                    toNeighbour <<
                        SubList<label>(nCellsUsingFace, pp.size(), pp.start());
                }
                else
                {
                    toNeighbour << labelList();
                }
            }
        }

        pBufs.finishedSends();

        // Receive face usage count and check for faces that become uncoupled.
        for (const polyPatch& pp : oldPatches)
        {
            const auto* procPatch = isA<processorPolyPatch>(pp);

            if (procPatch)
            {
                const label nbrProci = procPatch->neighbProcNo();

                UIPstream fromNeighbour(nbrProci, pBufs);

                const labelList nbrList(fromNeighbour);

                // Combine with this side.

                if (!nCellsUsingFace.empty())
                {
                    const labelList& nbrCellsUsingFace(nbrList);

                    // Combine with this side.

                    forAll(pp, i)
                    {
                        if
                        (
                            nCellsUsingFace[pp.start()+i] == 1
                         && nbrCellsUsingFace[i] == 0
                        )
                        {
                            // Face's neighbour is no longer there. Mark face
                            // off as coupled
                            nCellsUsingFace[pp.start()+i] = 3;
                            ++nUncoupled;
                        }
                    }
                }
            }
        }
    }

    // Do same for cyclics.
    for (const polyPatch& pp : oldPatches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && !nCellsUsingFace.empty())
        {
            const auto& cycPatch = *cpp;

            forAll(cycPatch, i)
            {
                label thisFacei = cycPatch.start() + i;
                label otherFacei = cycPatch.transformGlobalFace(thisFacei);

                if
                (
                    nCellsUsingFace[thisFacei] == 1
                 && nCellsUsingFace[otherFacei] == 0
                )
                {
                    nCellsUsingFace[thisFacei] = 3;
                    ++nUncoupled;
                }
            }
        }
    }

    if (syncPar)
    {
        reduce(nUncoupled, sumOp<label>());
    }

    if (nUncoupled > 0)
    {
        Info<< "Uncoupled " << nUncoupled << " faces on coupled patches. "
            << "(processorPolyPatch, cyclicPolyPatch)" << endl;
    }
}


void Foam::fvMeshSubset::removeCellsImpl
(
    const bitSet& cellsToRemove,
    const labelList& exposedFaces,
    const labelList& patchIDs,
    const bool syncCouples
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(baseMesh());

    removeCells cellRemover(baseMesh(), syncCouples);

    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        patchIDs,
        meshMod
    );

    // Create mesh, return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        fvMeshSubsetPtr_,
        IOobject
        (
            baseMesh().name(),
            baseMesh().time().timeName(),
            baseMesh().time(),
            IOobject::READ_IF_PRESENT,  // read fv* if present
            IOobject::NO_WRITE
        ),
        baseMesh(),
        syncCouples
    );

    pointMap_ = map().pointMap();
    faceMap_ = map().faceMap();
    cellMap_ = map().cellMap();
    patchMap_ = identity(baseMesh().boundaryMesh().size());
}


Foam::labelList Foam::fvMeshSubset::subsetSubset
(
    const label nElems,
    const labelUList& selectedElements,
    const labelUList& subsetMap
)
{
    // Mark selected elements.
    const bitSet selected(nElems, selectedElements);

    // Count subset of selected elements
    label n = 0;
    forAll(subsetMap, i)
    {
        if (selected[subsetMap[i]])
        {
            ++n;
        }
    }

    // Collect selected elements
    labelList subsettedElements(n);
    n = 0;

    forAll(subsetMap, i)
    {
        if (selected[subsetMap[i]])
        {
            subsettedElements[n] = i;
            ++n;
        }
    }

    return subsettedElements;
}


void Foam::fvMeshSubset::subsetZones()
{
    // Keep all zones, even if zero size.

    // PointZones

    const pointZoneMesh& pointZones = baseMesh().pointZones();

    List<pointZone*> pZonePtrs(pointZones.size());

    forAll(pointZones, zonei)
    {
        const pointZone& pz = pointZones[zonei];

        pZonePtrs[zonei] = new pointZone
        (
            pz.name(),
            subsetSubset(baseMesh().nPoints(), pz, pointMap()),
            zonei,
            fvMeshSubsetPtr_().pointZones()
        );
    }


    // FaceZones
    // Do we need to remove zones where the side we're interested in
    // no longer exists? Guess not.

    const faceZoneMesh& faceZones = baseMesh().faceZones();

    List<faceZone*> fZonePtrs(faceZones.size());

    forAll(faceZones, zonei)
    {
        const faceZone& fz = faceZones[zonei];

        // Expand faceZone to full mesh
        // +1 : part of faceZone, flipped
        // -1 :    ,,           , unflipped
        //  0 : not part of faceZone
        labelList zone(baseMesh().nFaces(), Zero);
        forAll(fz, j)
        {
            if (fz.flipMap()[j])
            {
                zone[fz[j]] = 1;
            }
            else
            {
                zone[fz[j]] = -1;
            }
        }

        // Select faces
        label nSub = 0;
        forAll(faceMap(), j)
        {
            if (zone[faceMap()[j]] != 0)
            {
                ++nSub;
            }
        }
        labelList subAddressing(nSub);
        boolList subFlipStatus(nSub);
        nSub = 0;
        forAll(faceMap(), subFacei)
        {
            const label meshFacei = faceMap()[subFacei];
            if (zone[meshFacei] != 0)
            {
                subAddressing[nSub] = subFacei;
                const label subOwner = subMesh().faceOwner()[subFacei];
                const label baseOwner = baseMesh().faceOwner()[meshFacei];
                // If subowner is the same cell as the base keep the flip status
                const bool sameOwner = (cellMap()[subOwner] == baseOwner);
                const bool flip = (zone[meshFacei] == 1);
                subFlipStatus[nSub] = (sameOwner == flip);

                ++nSub;
            }
        }

        fZonePtrs[zonei] = new faceZone
        (
            fz.name(),
            subAddressing,
            subFlipStatus,
            zonei,
            fvMeshSubsetPtr_().faceZones()
        );
    }

    // Cell Zones

    const cellZoneMesh& cellZones = baseMesh().cellZones();

    List<cellZone*> cZonePtrs(cellZones.size());

    forAll(cellZones, zonei)
    {
        const cellZone& cz = cellZones[zonei];

        cZonePtrs[zonei] = new cellZone
        (
            cz.name(),
            subsetSubset(baseMesh().nCells(), cz, cellMap()),
            zonei,
            fvMeshSubsetPtr_().cellZones()
        );
    }


    // Add the zones
    fvMeshSubsetPtr_().addZones(pZonePtrs, fZonePtrs, cZonePtrs);
}


Foam::bitSet Foam::fvMeshSubset::getCellsToRemove
(
    const bitSet& selectedCells
) const
{
    // Work on a copy
    bitSet cellsToRemove(selectedCells);

    // Ensure we have the full range
    cellsToRemove.resize(baseMesh().nCells(), false);

    // Invert the selection
    cellsToRemove.flip();

    return cellsToRemove;
}


Foam::bitSet Foam::fvMeshSubset::getCellsToRemove
(
    const label regioni,
    const labelUList& regions
) const
{
    return BitSetOps::create
    (
        baseMesh().nCells(),
        regioni,
        regions,
        false  // on=false: invert return cells to remove
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshSubset::fvMeshSubset(const fvMesh& baseMesh)
:
    baseMesh_(baseMesh),
    fvMeshSubsetPtr_(nullptr),
    faceFlipMapPtr_(nullptr),
    pointMap_(),
    faceMap_(),
    cellMap_(),
    patchMap_()
{}


Foam::fvMeshSubset::fvMeshSubset
(
    const fvMesh& baseMesh,
    const bitSet& selectedCells,
    const label patchID,
    const bool syncCouples
)
:
    fvMeshSubset(baseMesh)
{
    setCellSubset(selectedCells, patchID, syncCouples);
}


Foam::fvMeshSubset::fvMeshSubset
(
    const fvMesh& baseMesh,
    const labelHashSet& selectedCells,
    const label patchID,
    const bool syncCouples
)
:
    fvMeshSubset(baseMesh)
{
    setCellSubset(selectedCells, patchID, syncCouples);
}


Foam::fvMeshSubset::fvMeshSubset
(
    const fvMesh& baseMesh,
    const labelUList& selectedCells,
    const label patchID,
    const bool syncCouples
)
:
    fvMeshSubset(baseMesh)
{
    setCellSubset(selectedCells, patchID, syncCouples);
}


Foam::fvMeshSubset::fvMeshSubset
(
    const fvMesh& baseMesh,
    const label regioni,
    const labelUList& regions,
    const label patchID,
    const bool syncCouples
)
:
    fvMeshSubset(baseMesh)
{
    setCellSubset(regioni, regions, patchID, syncCouples);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshSubset::clear()
{
    fvMeshSubsetPtr_.reset(nullptr);
    faceFlipMapPtr_.reset(nullptr);

    pointMap_.clear();
    faceMap_.clear();
    cellMap_.clear();
    patchMap_.clear();
}


void Foam::fvMeshSubset::setCellSubset
(
    const bitSet& selectedCells,
    const label patchID,
    const bool syncPar
)
{
    const cellList& oldCells = baseMesh().cells();
    const faceList& oldFaces = baseMesh().faces();
    const pointField& oldPoints = baseMesh().points();
    const labelList& oldOwner = baseMesh().faceOwner();
    const labelList& oldNeighbour = baseMesh().faceNeighbour();
    const polyBoundaryMesh& oldPatches = baseMesh().boundaryMesh();
    const label oldNInternalFaces = baseMesh().nInternalFaces();

    // Initial check on patches before doing anything time consuming.

    label wantedPatchID = patchID;

    if (wantedPatchID == -1)
    {
        // No explicit patch specified. Put in oldInternalFaces patch.
        // Check if patch with this name already exists.
        wantedPatchID = oldPatches.findPatchID(exposedPatchName);
    }
    else if (wantedPatchID < 0 || wantedPatchID >= oldPatches.size())
    {
        FatalErrorInFunction
            << "Non-existing patch index " << wantedPatchID << endl
            << "Should be between 0 and " << oldPatches.size()-1
            << abort(FatalError);
    }

    // Clear all old maps and pointers
    clear();

    // The selected cells - sorted in ascending order
    cellMap_ = selectedCells.sortedToc();

    // The selectedCells should normally be in the [0,nCells) range,
    // but it is better not to trust that.
    {
        label len = cellMap_.size();
        for
        (
            label i = len-1;
            i >= 0 && (cellMap_[i] >= oldCells.size());
            --i
        )
        {
            len = i;
        }
        cellMap_.resize(len);
    }


    // Mark all used faces. Count number of cells using them
    // 0: face not used anymore
    // 1: face used by one cell, face becomes/stays boundary face
    // 2: face still used and remains internal face
    // 3: face coupled and used by one cell only (so should become normal,
    //    non-coupled patch face)
    //
    // Note that this is not really necessary - but means we can size things
    // correctly. Also makes handling coupled faces much easier.

    labelList nCellsUsingFace(oldFaces.size(), Zero);

    label nFacesInSet = 0;
    forAll(oldFaces, oldFacei)
    {
        bool faceUsed = false;

        if (selectedCells.test(oldOwner[oldFacei]))
        {
            ++nCellsUsingFace[oldFacei];
            faceUsed = true;
        }

        if
        (
            baseMesh().isInternalFace(oldFacei)
         && selectedCells.test(oldNeighbour[oldFacei])
        )
        {
            ++nCellsUsingFace[oldFacei];
            faceUsed = true;
        }

        if (faceUsed)
        {
            ++nFacesInSet;
        }
    }
    faceMap_.setSize(nFacesInSet);

    // Handle coupled faces. Modifies patch faces to be uncoupled to 3.
    doCoupledPatches(syncPar, nCellsUsingFace);


    // See which patch to use for exposed internal faces.
    label oldInternalPatchID = 0;

    // Insert faces before which patch
    label nextPatchID = oldPatches.size();

    // old to new patches
    labelList globalPatchMap(oldPatches.size());

    // New patch size
    label nbSize = oldPatches.size();

    if (wantedPatchID == -1)
    {
        // Create 'oldInternalFaces' patch at the end (or before
        // processorPatches)
        // and put all exposed internal faces in there.

        forAll(oldPatches, patchi)
        {
            if (isA<processorPolyPatch>(oldPatches[patchi]))
            {
                nextPatchID = patchi;
                break;
            }
            ++oldInternalPatchID;
        }

        ++nbSize;

        // adapt old to new patches for inserted patch
        for (label oldPatchi = 0; oldPatchi < nextPatchID; oldPatchi++)
        {
            globalPatchMap[oldPatchi] = oldPatchi;
        }
        for
        (
            label oldPatchi = nextPatchID;
            oldPatchi < oldPatches.size();
            oldPatchi++
        )
        {
            globalPatchMap[oldPatchi] = oldPatchi+1;
        }
    }
    else
    {
        oldInternalPatchID = wantedPatchID;
        nextPatchID = wantedPatchID+1;

        // old to new patches
        globalPatchMap = identity(oldPatches.size());
    }

    labelList boundaryPatchSizes(nbSize, Zero);


    // Make a global-to-local point map
    labelList globalPointMap(oldPoints.size(), -1);
    labelList globalFaceMap(oldFaces.size(), -1);

    label facei = 0;

    // 1. Pick up all preserved internal faces.
    for (label oldFacei = 0; oldFacei < oldNInternalFaces; ++oldFacei)
    {
        if (nCellsUsingFace[oldFacei] == 2)
        {
            globalFaceMap[oldFacei] = facei;
            faceMap_[facei++] = oldFacei;

            // Mark all points from the face
            markUsed(oldFaces[oldFacei], globalPointMap);
        }
    }

    // These are all the internal faces in the mesh.
    const label nInternalFaces = facei;

    // 2. Boundary faces up to where we want to insert old internal faces
    for
    (
        label oldPatchi = 0;
        oldPatchi < oldPatches.size()
     && oldPatchi < nextPatchID;
        oldPatchi++
    )
    {
        const polyPatch& oldPatch = oldPatches[oldPatchi];

        label oldFacei = oldPatch.start();

        forAll(oldPatch, i)
        {
            if (nCellsUsingFace[oldFacei] == 1)
            {
                // Boundary face is kept.

                // Mark face and increment number of points in set
                globalFaceMap[oldFacei] = facei;
                faceMap_[facei++] = oldFacei;

                // Mark all points from the face
                markUsed(oldFaces[oldFacei], globalPointMap);

                // Increment number of patch faces
                ++boundaryPatchSizes[globalPatchMap[oldPatchi]];
            }
            ++oldFacei;
        }
    }

    // 3a. old internal faces that have become exposed.
    for (label oldFacei = 0; oldFacei < oldNInternalFaces; ++oldFacei)
    {
        if (nCellsUsingFace[oldFacei] == 1)
        {
            globalFaceMap[oldFacei] = facei;
            faceMap_[facei++] = oldFacei;

            // Mark all points from the face
            markUsed(oldFaces[oldFacei], globalPointMap);

            // Increment number of patch faces
            ++boundaryPatchSizes[oldInternalPatchID];
        }
    }

    // 3b. coupled patch faces that have become uncoupled.
    for
    (
        label oldFacei = oldNInternalFaces;
        oldFacei < oldFaces.size();
        oldFacei++
    )
    {
        if (nCellsUsingFace[oldFacei] == 3)
        {
            globalFaceMap[oldFacei] = facei;
            faceMap_[facei++] = oldFacei;

            // Mark all points from the face
            markUsed(oldFaces[oldFacei], globalPointMap);

            // Increment number of patch faces
            ++boundaryPatchSizes[oldInternalPatchID];
        }
    }

    // 4. Remaining boundary faces
    for
    (
        label oldPatchi = nextPatchID;
        oldPatchi < oldPatches.size();
        oldPatchi++
    )
    {
        const polyPatch& oldPatch = oldPatches[oldPatchi];

        label oldFacei = oldPatch.start();

        forAll(oldPatch, i)
        {
            if (nCellsUsingFace[oldFacei] == 1)
            {
                // Boundary face is kept.

                // Mark face and increment number of points in set
                globalFaceMap[oldFacei] = facei;
                faceMap_[facei++] = oldFacei;

                // Mark all points from the face
                markUsed(oldFaces[oldFacei], globalPointMap);

                // Increment number of patch faces
                ++boundaryPatchSizes[globalPatchMap[oldPatchi]];
            }
            ++oldFacei;
        }
    }

    if (facei != nFacesInSet)
    {
        FatalErrorInFunction
            << "Problem" << abort(FatalError);
    }


    // Grab the points map
    label nPointsInSet = 0;

    forAll(globalPointMap, pointi)
    {
        if (globalPointMap[pointi] != -1)
        {
            ++nPointsInSet;
        }
    }
    pointMap_.setSize(nPointsInSet);

    nPointsInSet = 0;

    forAll(globalPointMap, pointi)
    {
        if (globalPointMap[pointi] != -1)
        {
            pointMap_[nPointsInSet] = pointi;
            globalPointMap[pointi] = nPointsInSet;
            ++nPointsInSet;
        }
    }


    //Pout<< "Number of points,faces,cells in new mesh : "
    //    << pointMap_.size() << ' '
    //    << faceMap_.size() << ' '
    //    << cellMap_.size() << nl;


    // Make a new mesh

    //
    // Create new points
    //
    pointField newPoints(pointUIndList(oldPoints, pointMap_));


    //
    // Create new faces
    //

    faceList newFaces(faceMap_.size());
    {
        auto iter = newFaces.begin();
        const auto& renumbering = globalPointMap;

        // Internal faces
        for (label facei = 0; facei < nInternalFaces; ++facei)
        {
            face& newItem = *iter;
            ++iter;

            const face& oldItem = oldFaces[faceMap_[facei]];

            // Copy relabelled
            newItem.resize(oldItem.size());

            forAll(oldItem, i)
            {
                newItem[i] = renumbering[oldItem[i]];
            }
        }

        // Boundary faces - may need to be flipped
        for (label facei = nInternalFaces; facei < faceMap_.size(); ++facei)
        {
            const label oldFacei = faceMap_[facei];

            face& newItem = *iter;
            ++iter;

            const face oldItem =
            (
                (
                    baseMesh().isInternalFace(oldFacei)
                 && selectedCells.test(oldNeighbour[oldFacei])
                 && !selectedCells.test(oldOwner[oldFacei])
                )
                // Face flipped to point outwards
              ? oldFaces[oldFacei].reverseFace()
              : oldFaces[oldFacei]
            );

            // Copy relabelled
            newItem.resize(oldItem.size());

            forAll(oldItem, i)
            {
                newItem[i] = renumbering[oldItem[i]];
            }
        }
    }


    //
    // Create new cells
    //

    cellList newCells(cellMap_.size());
    {
        auto iter = newCells.begin();
        const auto& renumbering = globalFaceMap;

        for (const label oldCelli : cellMap_)
        {
            cell& newItem = *iter;
            ++iter;

            const labelList& oldItem = oldCells[oldCelli];

            // Copy relabelled
            newItem.resize(oldItem.size());

            forAll(oldItem, i)
            {
                newItem[i] = renumbering[oldItem[i]];
            }
        }
    }


    // Make a new mesh
    //
    // Note that mesh gets registered with same name as original mesh.
    // This is not proper but cannot be avoided since otherwise
    // surfaceInterpolation cannot find its fvSchemes.
    // It will try to read, for example.  "system/region0SubSet/fvSchemes"
    //
    fvMeshSubsetPtr_ = autoPtr<fvMesh>::New
    (
        IOobject
        (
            baseMesh().name(),
            baseMesh().time().timeName(),
            baseMesh().time(),
            IOobject::NO_READ,      // do not read any dictionaries
            IOobject::NO_WRITE
        ),
        baseMesh(),                 // get dictionaries from base mesh
        std::move(newPoints),
        std::move(newFaces),
        std::move(newCells),
        syncPar           // parallel synchronisation
    );


    // Add old patches
    List<polyPatch*> newBoundary(nbSize);
    patchMap_.setSize(nbSize);
    label nNewPatches = 0;
    label patchStart = nInternalFaces;


    // For parallel: only remove patch if none of the processors has it.
    // This only gets done for patches before the one being inserted
    // (so patches < nextPatchID)

    // Get sum of patch sizes. Zero if patch can be deleted.
    labelList globalPatchSizes(boundaryPatchSizes);
    globalPatchSizes.setSize(nextPatchID);

    if (syncPar && Pstream::parRun())
    {
        // Get patch names (up to nextPatchID)
        List<wordList> patchNames(Pstream::nProcs());
        patchNames[Pstream::myProcNo()] = oldPatches.names();
        patchNames[Pstream::myProcNo()].setSize(nextPatchID);
        Pstream::gatherList(patchNames);
        Pstream::scatterList(patchNames);

        // Get patch sizes (up to nextPatchID).
        // Note that up to nextPatchID the globalPatchMap is an identity so
        // no need to index through that.
        Pstream::listCombineGather(globalPatchSizes, plusEqOp<label>());
        Pstream::listCombineScatter(globalPatchSizes);

        // Now all processors have all the patchnames.
        // Decide: if all processors have the same patch names and size is zero
        // everywhere remove the patch.
        bool samePatches = true;

        for (label proci = 1; proci < patchNames.size(); ++proci)
        {
            if (patchNames[proci] != patchNames[0])
            {
                samePatches = false;
                break;
            }
        }

        if (!samePatches)
        {
            // Patchnames not sync on all processors so disable removal of
            // zero sized patches.
            globalPatchSizes = labelMax;
        }
    }


    // Old patches

    for
    (
        label oldPatchi = 0;
        oldPatchi < oldPatches.size()
     && oldPatchi < nextPatchID;
        oldPatchi++
    )
    {
        const label newSize = boundaryPatchSizes[globalPatchMap[oldPatchi]];

        if (oldInternalPatchID != oldPatchi)
        {
            // Pure subset of patch faces (no internal faces added to this
            // patch). Can use mapping.
            labelList map(newSize);
            for (label patchFacei = 0; patchFacei < newSize; patchFacei++)
            {
                const label facei = patchStart+patchFacei;
                const label oldFacei = faceMap_[facei];
                map[patchFacei] = oldPatches[oldPatchi].whichFace(oldFacei);
            }

            newBoundary[nNewPatches] = oldPatches[oldPatchi].clone
            (
                fvMeshSubsetPtr_().boundaryMesh(),
                nNewPatches,
                map,
                patchStart
            ).ptr();
        }
        else
        {
            // Clone (even if 0 size)
            newBoundary[nNewPatches] = oldPatches[oldPatchi].clone
            (
                fvMeshSubsetPtr_().boundaryMesh(),
                nNewPatches,
                newSize,
                patchStart
            ).ptr();
        }

        patchStart += newSize;
        patchMap_[nNewPatches] = oldPatchi;    // compact patchMap
        ++nNewPatches;
    }

    // Inserted patch

    if (wantedPatchID == -1)
    {
        label oldInternalSize = boundaryPatchSizes[oldInternalPatchID];

        if (syncPar)
        {
            reduce(oldInternalSize, sumOp<label>());
        }

        // Newly created patch so is at end. Check if any faces in it.
        if (oldInternalSize > 0)
        {
            newBoundary[nNewPatches] = new emptyPolyPatch
            (
                exposedPatchName,
                boundaryPatchSizes[oldInternalPatchID],
                patchStart,
                nNewPatches,
                fvMeshSubsetPtr_().boundaryMesh(),
                emptyPolyPatch::typeName
            );

            //Pout<< "    " << exposedPatchName << " : "
            //    << boundaryPatchSizes[oldInternalPatchID] << endl;

            // The index for the first patch is -1 as it originates from
            // the internal faces
            patchStart += boundaryPatchSizes[oldInternalPatchID];
            patchMap_[nNewPatches] = -1;
            ++nNewPatches;
        }
    }

    // Old patches

    for
    (
        label oldPatchi = nextPatchID;
        oldPatchi < oldPatches.size();
        oldPatchi++
    )
    {
        const label newSize = boundaryPatchSizes[globalPatchMap[oldPatchi]];

        if (oldInternalPatchID != oldPatchi)
        {
            // Pure subset of patch faces (no internal faces added to this
            // patch). Can use mapping.
            labelList map(newSize);
            for (label patchFacei = 0; patchFacei < newSize; patchFacei++)
            {
                const label facei = patchStart+patchFacei;
                const label oldFacei = faceMap_[facei];
                map[patchFacei] = oldPatches[oldPatchi].whichFace(oldFacei);
            }

            newBoundary[nNewPatches] = oldPatches[oldPatchi].clone
            (
                fvMeshSubsetPtr_().boundaryMesh(),
                nNewPatches,
                map,
                patchStart
            ).ptr();
        }
        else
        {
            // Patch still exists. Add it
            newBoundary[nNewPatches] = oldPatches[oldPatchi].clone
            (
                fvMeshSubsetPtr_().boundaryMesh(),
                nNewPatches,
                newSize,
                patchStart
            ).ptr();
        }

        //Pout<< "    " << oldPatches[oldPatchi].name() << " : "
        //    << newSize << endl;

        patchStart += newSize;
        patchMap_[nNewPatches] = oldPatchi;    // compact patchMap
        ++nNewPatches;
    }


    // Reset the patch lists
    newBoundary.setSize(nNewPatches);
    patchMap_.setSize(nNewPatches);


    // Add the fvPatches
    fvMeshSubsetPtr_().addFvPatches(newBoundary, syncPar);

    // Subset and add any zones
    subsetZones();
}


void Foam::fvMeshSubset::setCellSubset
(
    const labelUList& selectedCells,
    const label patchID,
    const bool syncPar
)
{
    setCellSubset
    (
        BitSetOps::create(baseMesh().nCells(), selectedCells),
        patchID,
        syncPar
    );
}


void Foam::fvMeshSubset::setCellSubset
(
    const labelHashSet& selectedCells,
    const label patchID,
    const bool syncPar
)
{
    setCellSubset
    (
        BitSetOps::create(baseMesh().nCells(), selectedCells),
        patchID,
        syncPar
    );
}


void Foam::fvMeshSubset::setCellSubset
(
    const label regioni,
    const labelUList& regions,
    const label patchID,
    const bool syncPar
)
{
    setCellSubset
    (
        BitSetOps::create(baseMesh().nCells(), regioni, regions),
        patchID,
        syncPar
    );
}


Foam::labelList Foam::fvMeshSubset::getExposedFaces
(
    const bitSet& selectedCells,
    const bool syncCouples
) const
{
    return
        removeCells(baseMesh(), syncCouples)
       .getExposedFaces(getCellsToRemove(selectedCells));
}


Foam::labelList Foam::fvMeshSubset::getExposedFaces
(
    const label regioni,
    const labelUList& regions,
    const bool syncCouples
) const
{
    return
        removeCells(baseMesh(), syncCouples)
       .getExposedFaces(getCellsToRemove(regioni, regions));
}


void Foam::fvMeshSubset::setCellSubset
(
    const bitSet& selectedCells,
    const labelList& exposedFaces,
    const labelList& patchIDs,
    const bool syncCouples
)
{
    removeCellsImpl
    (
        getCellsToRemove(selectedCells),
        exposedFaces,
        patchIDs,
        syncCouples
    );
}


void Foam::fvMeshSubset::setCellSubset
(
    const label selectRegion,
    const labelList& regions,
    const labelList& exposedFaces,
    const labelList& patchIDs,
    const bool syncCouples
)
{
    removeCellsImpl
    (
        getCellsToRemove(selectRegion, regions),
        exposedFaces,
        patchIDs,
        syncCouples
    );
}


// ************************************************************************* //
