/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "InteractionLists.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::InteractionLists<ParticleType>::findExtendedProcBbsInRange
(
    const treeBoundBox& procBb,
    const List<treeBoundBox>& allExtendedProcBbs,
    const globalIndexAndTransform& globalTransforms,
    List<treeBoundBox>& extendedProcBbsInRange,
    List<label>& extendedProcBbsTransformIndex,
    List<label>& extendedProcBbsOrigProc
)
{
    extendedProcBbsInRange.setSize(0);
    extendedProcBbsTransformIndex.setSize(0);
    extendedProcBbsOrigProc.setSize(0);

    DynamicList<treeBoundBox> tmpExtendedProcBbsInRange;
    DynamicList<label> tmpExtendedProcBbsTransformIndex;
    DynamicList<label> tmpExtendedProcBbsOrigProc;

    label nTrans = globalTransforms.nIndependentTransforms();

    forAll(allExtendedProcBbs, procI)
    {
        List<label> permutationIndices(nTrans, 0);

        vector s = vector::zero;

        if (nTrans == 0 && procI != Pstream::myProcNo())
        {
            treeBoundBox extendedReferredProcBb = allExtendedProcBbs[procI];

            if (procBb.overlaps(extendedReferredProcBb))
            {
                tmpExtendedProcBbsInRange.append
                (
                    extendedReferredProcBb
                );

                // Dummy index, there are no transforms, so there will
                // be no resultant transform when this is decoded.
                tmpExtendedProcBbsTransformIndex.append(0);

                tmpExtendedProcBbsOrigProc.append(procI);
            }
        }
        else if (nTrans == 3)
        {
            label& i = permutationIndices[0];
            label& j = permutationIndices[1];
            label& k = permutationIndices[2];

            for (i = -1; i <= 1; i++)
            {
                for (j = -1; j <= 1; j++)
                {
                    for (k = -1; k <= 1; k++)
                    {
                        if
                        (
                            i == 0
                         && j == 0
                         && k == 0
                         && procI == Pstream::myProcNo()
                        )
                        {
                            // Skip this processor's extended boundBox
                            // when it has no transformation
                            continue;
                        }

                        label transI = globalTransforms.encodeTransformIndex
                        (
                            permutationIndices
                        );

                        const vector& transform = globalTransforms.transform
                        (
                            transI
                        );

                        treeBoundBox extendedReferredProcBb
                        (
                            allExtendedProcBbs[procI].min() + transform,
                            allExtendedProcBbs[procI].max() + transform
                        );

                        if (procBb.overlaps(extendedReferredProcBb))
                        {
                            tmpExtendedProcBbsInRange.append
                            (
                                extendedReferredProcBb
                            );

                            tmpExtendedProcBbsTransformIndex.append(transI);

                            tmpExtendedProcBbsOrigProc.append(procI);
                        }
                    }
                }
            }
        }
        else if (nTrans == 2)
        {
            label& i = permutationIndices[0];
            label& j = permutationIndices[1];

            for (i = -1; i <= 1; i++)
            {
                for (j = -1; j <= 1; j++)
                {
                    if (i == 0 && j == 0 && procI == Pstream::myProcNo())
                    {
                        // Skip this processor's extended boundBox
                        // when it has no transformation
                        continue;
                    }

                    label transI = globalTransforms.encodeTransformIndex
                    (
                        permutationIndices
                    );

                    const vector& transform = globalTransforms.transform
                    (
                        transI
                    );

                    treeBoundBox extendedReferredProcBb
                    (
                        allExtendedProcBbs[procI].min() + transform,
                        allExtendedProcBbs[procI].max() + transform
                    );

                    if (procBb.overlaps(extendedReferredProcBb))
                    {
                        tmpExtendedProcBbsInRange.append
                        (
                            extendedReferredProcBb
                        );

                        tmpExtendedProcBbsTransformIndex.append(transI);

                        tmpExtendedProcBbsOrigProc.append(procI);
                    }
                }
            }
        }
        else if (nTrans == 1)
        {
            label& i = permutationIndices[0];

            for (i = -1; i <= 1; i++)
            {
                if (i == 0 && procI == Pstream::myProcNo())
                {
                    // Skip this processor's extended boundBox when it
                    // has no transformation
                    continue;
                }

                label transI = globalTransforms.encodeTransformIndex
                (
                    permutationIndices
                );

                const vector& transform = globalTransforms.transform(transI);

                treeBoundBox extendedReferredProcBb
                (
                    allExtendedProcBbs[procI].min() + transform,
                    allExtendedProcBbs[procI].max() + transform
                );

                if (procBb.overlaps(extendedReferredProcBb))
                {
                    tmpExtendedProcBbsInRange.append
                    (
                        extendedReferredProcBb
                    );

                    tmpExtendedProcBbsTransformIndex.append(transI);

                    tmpExtendedProcBbsOrigProc.append(procI);
                }
            }
        }
    }

    extendedProcBbsInRange = tmpExtendedProcBbsInRange.xfer();
    extendedProcBbsTransformIndex = tmpExtendedProcBbsTransformIndex.xfer();
    extendedProcBbsOrigProc = tmpExtendedProcBbsOrigProc.xfer();
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::buildMap
(
    const List<label>& toProc
)
{
    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // 1. Count
    labelList nSend(Pstream::nProcs(), 0);

    forAll(toProc, i)
    {
        label procI = toProc[i];

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
    forAll(toProc, i)
    {
        label procI = toProc[i];

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

    mapPtr_.reset
    (
        new mapDistribute
        (

            constructSize,
            sendMap.xfer(),
            constructMap.xfer()
        )
    );
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::prepareParticlesToRefer
(
    const List<DynamicList<ParticleType*> >& cellOccupancy
)
{
    referredParticles_.setSize(cellIndexAndTransformToDistribute_.size());

    // Clear all existing referred particles

    forAll(referredParticles_, i)
    {
        referredParticles_[i].clear();
    }

    forAll(cellIndexAndTransformToDistribute_, i)
    {
        const labelPair giat = cellIndexAndTransformToDistribute_[i];

        label cellIndex = globalTransforms_.index(giat);

        List<ParticleType*> realParticles = cellOccupancy[cellIndex];

        IDLList<ParticleType>& particlesToRefer = referredParticles_[i];

        forAll (realParticles, rM)
        {
            const ParticleType& particle = *realParticles[rM];

            particlesToRefer.append(particle.clone().ptr());

            prepareParticleToBeReferred(particlesToRefer.last(), giat);
        }
    }
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::prepareParticleToBeReferred
(
    ParticleType* particle,
    labelPair giat
)
{
    const vector& transform = globalTransforms_.transform
    (
        globalTransforms_.transformIndex(giat)
    );

    particle->position() -= transform;
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::writeReferredParticleCloud()
{
    bool writeCloud = true;

    if (writeCloud)
    {
        cloud_.clear();

        forAll(referredParticles_, refCellI)
        {
            const IDLList<ParticleType>& refCell = referredParticles_[refCellI];

            forAllConstIter(typename IDLList<ParticleType>, refCell, iter)
            {
                cloud_.addParticle(iter().clone().ptr());
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::InteractionLists<ParticleType>::InteractionLists
(
    const polyMesh& mesh,
    scalar maxDistance
)
:
    mesh_(mesh),
    cloud_(mesh_, "referredParticleCloud", IDLList<ParticleType>()),
    mapPtr_(),
    globalTransforms_(mesh_),
    maxDistance_(maxDistance),
    dil_(),
    directWallFaces_(),
    ril_(),
    rilInverse_(),
    cellIndexAndTransformToDistribute_(),
    referredParticles_()
{
    Info<< "Building InteractionLists" << endl;

    Random rndGen(419715);

    const vector interactionVec = maxDistance_*vector::one;

    treeBoundBox procBb(treeBoundBox(mesh_.points()));

    treeBoundBox extendedProcBb
    (
        procBb.min() - interactionVec,
        procBb.max() + interactionVec
    );

    treeBoundBoxList allExtendedProcBbs(Pstream::nProcs());

    allExtendedProcBbs[Pstream::myProcNo()] = extendedProcBb;

    Pstream::gatherList(allExtendedProcBbs);

    Pstream::scatterList(allExtendedProcBbs);

    List<treeBoundBox> extendedProcBbsInRange;
    List<label> extendedProcBbsTransformIndex;
    List<label> extendedProcBbsOrigProc;

    findExtendedProcBbsInRange
    (
        procBb,
        allExtendedProcBbs,
        globalTransforms_,
        extendedProcBbsInRange,
        extendedProcBbsTransformIndex,
        extendedProcBbsOrigProc
    );

    treeBoundBoxList cellBbs(mesh_.nCells());

    forAll(cellBbs, cellI)
    {
        cellBbs[cellI] = treeBoundBox
        (
            mesh_.cells()[cellI].points
            (
                mesh_.faces(),
                mesh_.points()
            )
        );
    }

    // Recording which cells are in range of an extended boundBox, as
    // only these cells will need to be tested to determine which
    // referred cells that they interact with.
    PackedBoolList cellInRangeOfCoupledPatch(mesh_.nCells(), false);

    // IAndT: index and transform
    DynamicList<labelPair> globalIAndTToExchange;

    DynamicList<treeBoundBox> bbsToExchange;

    DynamicList<label> procToDistributeTo;

    forAll(extendedProcBbsInRange, ePBIRI)
    {
        const treeBoundBox& otherExtendedProcBb =
        extendedProcBbsInRange[ePBIRI];

        label transformIndex = extendedProcBbsTransformIndex[ePBIRI];

        label origProc = extendedProcBbsOrigProc[ePBIRI];

        forAll(cellBbs, cellI)
        {
            const treeBoundBox& cellBb = cellBbs[cellI];

            if (cellBb.overlaps(otherExtendedProcBb))
            {
                // This cell is in range of the Bb of the other
                // processor Bb, and so needs to be referred to it

                cellInRangeOfCoupledPatch[cellI] = true;

                globalIAndTToExchange.append
                (
                    globalTransforms_.encode(cellI, transformIndex)
                );

                bbsToExchange.append(cellBb);

                procToDistributeTo.append(origProc);
            }
        }
    }

    buildMap(procToDistributeTo);

    // Needed for reverseDistribute
    label preDistributionSize = procToDistributeTo.size();

    map().distribute(bbsToExchange);

    map().distribute(globalIAndTToExchange);

    // Determine labelList specifying only cells that are in range of
    // a coupled boundary to build an octree limited to these cells.
    DynamicList<label> coupledPatchRangeCells;

    forAll(cellInRangeOfCoupledPatch, cellI)
    {
        if (cellInRangeOfCoupledPatch[cellI])
        {
            coupledPatchRangeCells.append(cellI);
        }
    }

    treeBoundBox procBbRndExt
    (
        treeBoundBox(mesh_.points()).extend(rndGen, 1e-4)
    );

    indexedOctree<treeDataCell> coupledPatchRangeTree
    (
        treeDataCell(true, mesh_, coupledPatchRangeCells),
        procBbRndExt,
        8,              // maxLevel,
        10,             // leafSize,
        100.0
    );

    ril_.setSize(bbsToExchange.size());

    // This needs to be a boolList, not PackedBoolList if
    // reverseDistribute is called.
    boolList bbRequiredByAnyCell(bbsToExchange.size(), false);

    Info<< "    Building referred interaction lists" << endl;

    forAll(bbsToExchange, bbI)
    {
        const labelPair& giat = globalIAndTToExchange[bbI];

        const vector& transform = globalTransforms_.transform
        (
            globalTransforms_.transformIndex(giat)
        );

        treeBoundBox extendedBb
        (
            bbsToExchange[bbI].min() - interactionVec - transform,
            bbsToExchange[bbI].max() + interactionVec - transform
        );

        // Find all elements intersecting box.
        labelList interactingElems
        (
            coupledPatchRangeTree.findBox(extendedBb)
        );

        if (!interactingElems.empty())
        {
            bbRequiredByAnyCell[bbI] = true;
        }

        ril_[bbI].setSize(interactingElems.size(), -1);

        forAll(interactingElems, i)
        {
            label elemI = interactingElems[i];

            // Here, a more detailed geometric test could be applied,
            // i.e. a more accurate bounding volume like a OBB or
            // convex hull or an exact geometrical test.

            label c = coupledPatchRangeTree.shapes().cellLabels()[elemI];

            ril_[bbI][i] = c;
        }
    }

    // Perform subset of ril_, to remove any referred cells that do
    // not interact.  They will not be sent from originating
    // processors.  This assumes that the disappearance of the cell
    // from the sending list of the source processor, simply removes
    // the referred cell from the ril_, all of the subsequent indices
    // shuffle down one, but the order and structure is preserved,
    // i.e. it, is as if the cell had never been referred in the first
    // place.

    inplaceSubset(bbRequiredByAnyCell, ril_);

    // Send information about which cells are actually required back
    // to originating processors.

    // At this point, bbsToExchange does not need to be maintained
    // or distributed as it is not longer needed.

    bbsToExchange.setSize(0);

    map().reverseDistribute
    (
        preDistributionSize,
        bbRequiredByAnyCell
    );

    map().reverseDistribute
    (
        preDistributionSize,
        globalIAndTToExchange
    );

    // Perform ordering of cells to send, this invalidates the
    // previous value of preDistributionSize.

    preDistributionSize = -1;

    // Move all of the used cells to refer to the start of the list
    // and truncate it

    inplaceSubset(bbRequiredByAnyCell, globalIAndTToExchange);

    inplaceSubset(bbRequiredByAnyCell, procToDistributeTo);

    preDistributionSize = procToDistributeTo.size();

    // Rebuild mapDistribute with only required referred cells
    buildMap(procToDistributeTo);

    // Store cellIndexAndTransformToDistribute
    cellIndexAndTransformToDistribute_.transfer(globalIAndTToExchange);

    // Determine inverse addressing for referred cells

    rilInverse_.setSize(mesh_.nCells());

    // Temporary Dynamic lists for accumulation
    List<DynamicList<label> > rilInverseTemp(rilInverse_.size());

    // Loop over all referred cells
    forAll(ril_, refCellI)
    {
        const labelList& realCells = ril_[refCellI];

        // Loop over all real cells in that the referred cell is to
        // supply interactions to and record the index of this
        // referred cell in the real cells entry in rilInverse

        forAll(realCells, realCellI)
        {
            rilInverseTemp[realCells[realCellI]].append(refCellI);
        }
    }

    forAll(rilInverse_, cellI)
    {
        rilInverse_[cellI].transfer(rilInverseTemp[cellI]);
    }

    // Direct interaction list

    Info<< "    Building direct interaction lists" << endl;

    indexedOctree<treeDataCell> allCellsTree
    (
        treeDataCell(true, mesh_),
        procBbRndExt,
        8,              // maxLevel,
        10,             // leafSize,
        100.0
    );

    dil_.setSize(mesh_.nCells());

    forAll(cellBbs, cellI)
    {
        const treeBoundBox& cellBb = cellBbs[cellI];

        treeBoundBox extendedBb
        (
            cellBb.min() - interactionVec,
            cellBb.max() + interactionVec
        );

        // Find all elements intersecting box.
        labelList interactingElems
        (
            allCellsTree.findBox(extendedBb)
        );

        // Reserve space to avoid multiple resizing
        DynamicList<label> cellDIL(interactingElems.size());

        forAll(interactingElems, i)
        {
            label elemI = interactingElems[i];

            label c = allCellsTree.shapes().cellLabels()[elemI];

            // Here, a more detailed geometric test could be applied,
            // i.e. a more accurate bounding volume like a OBB or
            // convex hull, or an exact geometrical test.

            // The higher index cell is added to the lower index
            // cell's DIL.  A cell is not added to its own DIL.
            if (c > cellI)
            {
                cellDIL.append(c);
            }
        }

        dil_[cellI].transfer(cellDIL);
    }

    // Direct wall faces

    // DynamicLists for data gathering
    DynamicList<label> thisCellOnlyWallFaces;
    DynamicList<label> otherCellOnlyWallFaces;
    List<DynamicList<label> > wallFacesTemp(mesh_.nCells());

    const labelList& patchID = mesh_.boundaryMesh().patchID();

    label nInternalFaces = mesh_.nInternalFaces();

    forAll(wallFacesTemp, thisCellI)
    {
        // Find all of the wall faces for the current cell
        const labelList& thisCellFaces = mesh_.cells()[thisCellI];

        DynamicList<label>& thisCellWallFaces = wallFacesTemp[thisCellI];

        thisCellOnlyWallFaces.clear();

        forAll(thisCellFaces, tCFI)
        {
            label faceI = thisCellFaces[tCFI];

            if (!mesh_.isInternalFace(faceI))
            {
                label patchI = patchID[faceI - nInternalFaces];

                const polyPatch& patch = mesh_.boundaryMesh()[patchI];

                if (isA<wallPolyPatch>(patch))
                {
                    thisCellOnlyWallFaces.append(faceI);
                }
            }
        }

        // Add all the found wall faces to this cell's list, and
        // retain the wall faces for this cell only to add to other
        // cells.
        thisCellWallFaces.append(thisCellOnlyWallFaces);

        // Loop over all of the cells in the dil for this cell, adding
        // the wallFaces for this cell to the other cell's wallFace
        // list, and all of the wallFaces for the other cell to this
        // cell's list

        const labelList& cellDil = dil_[thisCellI];

        forAll(cellDil, i)
        {
            label otherCellI = cellDil[i];

            const labelList& otherCellFaces = mesh_.cells()[otherCellI];

            DynamicList<label>& otherCellWallFaces = wallFacesTemp[otherCellI];

            otherCellOnlyWallFaces.clear();

            forAll(otherCellFaces, oCFI)
            {
                label faceI = otherCellFaces[oCFI];

                if (!mesh_.isInternalFace(faceI))
                {
                    label patchI = patchID[faceI - nInternalFaces];

                    const polyPatch& patch = mesh_.boundaryMesh()[patchI];

                    if (isA<wallPolyPatch>(patch))
                    {
                        otherCellOnlyWallFaces.append(faceI);
                    }
                }
            }

            thisCellWallFaces.append(otherCellOnlyWallFaces);

            otherCellWallFaces.append(thisCellOnlyWallFaces);
        }
    }

    directWallFaces_.setSize(mesh_.nCells());

    forAll(directWallFaces_, i)
    {
        directWallFaces_[i].transfer(wallFacesTemp[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::InteractionLists<ParticleType>::~InteractionLists()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::InteractionLists<ParticleType>::sendReferredParticles
(
    const List<DynamicList<ParticleType*> >& cellOccupancy,
    PstreamBuffers& pBufs
)
{
    prepareParticlesToRefer(cellOccupancy);

    // Stream data into buffer
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& subMap = map().subMap()[domain];

        if (subMap.size())
        {
            // Put data into send buffer
            UOPstream toDomain(domain, pBufs);

            UIndirectList<IDLList<ParticleType> > subMappedParticles
            (
                referredParticles_,
                subMap
            );

            forAll(subMappedParticles, i)
            {
                toDomain << subMappedParticles[i];
            }
        }
    }

    // Start sending and receiving but do not block.
    pBufs.finishedSends(false);
};


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::receiveReferredParticles
(
    PstreamBuffers& pBufs
)
{
    Pstream::waitRequests();

    referredParticles_.setSize(map().constructSize());

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& constructMap = map().constructMap()[domain];

        if (constructMap.size())
        {
            UIPstream str(domain, pBufs);

            forAll (constructMap, i)
            {
                referredParticles_[constructMap[i]] = IDLList<ParticleType>
                (
                    str,
                    typename ParticleType::iNew(cloud_)
                );
            }
        }
    }

    writeReferredParticleCloud();
}


// ************************************************************************* //
