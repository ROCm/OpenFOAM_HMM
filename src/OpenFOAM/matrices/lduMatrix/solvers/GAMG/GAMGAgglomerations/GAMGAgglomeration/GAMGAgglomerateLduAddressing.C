/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "GAMGAgglomeration.H"
#include "GAMGInterface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GAMGAgglomeration::agglomerateLduAddressing
(
    const label fineLevelIndex
)
{
    const lduMesh& fineMesh = meshLevel(fineLevelIndex);
    const lduAddressing& fineMeshAddr = fineMesh.lduAddr();

    const labelUList& upperAddr = fineMeshAddr.upperAddr();
    const labelUList& lowerAddr = fineMeshAddr.lowerAddr();

    label nFineFaces = upperAddr.size();

    // Get restriction map for current level
    const labelField& restrictMap = restrictAddressing(fineLevelIndex);

    if (min(restrictMap) == -1)
    {
        FatalErrorIn("GAMGAgglomeration::agglomerateLduAddressing")
            << "min(restrictMap) == -1" << exit(FatalError);
    }

    if (restrictMap.size() != fineMeshAddr.size())
    {
        FatalErrorIn
        (
            "GAMGAgglomeration::agglomerateLduAddressing"
            "(const label fineLevelIndex)"
        )   << "restrict map does not correspond to fine level. " << endl
            << " Sizes: restrictMap: " << restrictMap.size()
            << " nEqns: " << fineMeshAddr.size()
            << abort(FatalError);
    }


    // Get the number of coarse cells
    const label nCoarseCells = nCells_[fineLevelIndex];

    // Storage for coarse cell neighbours and coefficients

    // Guess initial maximum number of neighbours in coarse cell
    label maxNnbrs = 10;

    // Number of faces for each coarse-cell
    labelList cCellnFaces(nCoarseCells, 0);

    // Setup initial packed storage for coarse-cell faces
    labelList cCellFaces(maxNnbrs*nCoarseCells);

    // Create face-restriction addressing
    faceRestrictAddressing_.set(fineLevelIndex, new labelList(nFineFaces));
    labelList& faceRestrictAddr = faceRestrictAddressing_[fineLevelIndex];

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineFaces);

    // Counter for coarse faces
    label& nCoarseFaces = nFaces_[fineLevelIndex];
    nCoarseFaces = 0;

    // Loop through all fine faces
    forAll(upperAddr, fineFacei)
    {
        label rmUpperAddr = restrictMap[upperAddr[fineFacei]];
        label rmLowerAddr = restrictMap[lowerAddr[fineFacei]];

        if (rmUpperAddr == rmLowerAddr)
        {
            // For each fine face inside of a coarse cell keep the address
            // of the cell corresponding to the face in the faceRestrictAddr
            // as a negative index
            faceRestrictAddr[fineFacei] = -(rmUpperAddr + 1);
        }
        else
        {
            // this face is a part of a coarse face

            label cOwn = rmUpperAddr;
            label cNei = rmLowerAddr;

            // get coarse owner and neighbour
            if (rmUpperAddr > rmLowerAddr)
            {
                cOwn = rmLowerAddr;
                cNei = rmUpperAddr;
            }

            // check the neighbour to see if this face has already been found
            label* ccFaces = &cCellFaces[maxNnbrs*cOwn];

            bool nbrFound = false;
            label& ccnFaces = cCellnFaces[cOwn];

            for (int i=0; i<ccnFaces; i++)
            {
                if (initCoarseNeighb[ccFaces[i]] == cNei)
                {
                    nbrFound = true;
                    faceRestrictAddr[fineFacei] = ccFaces[i];
                    break;
                }
            }

            if (!nbrFound)
            {
                if (ccnFaces >= maxNnbrs)
                {
                    label oldMaxNnbrs = maxNnbrs;
                    maxNnbrs *= 2;

                    cCellFaces.setSize(maxNnbrs*nCoarseCells);

                    forAllReverse(cCellnFaces, i)
                    {
                        label* oldCcNbrs = &cCellFaces[oldMaxNnbrs*i];
                        label* newCcNbrs = &cCellFaces[maxNnbrs*i];

                        for (int j=0; j<cCellnFaces[i]; j++)
                        {
                            newCcNbrs[j] = oldCcNbrs[j];
                        }
                    }

                    ccFaces = &cCellFaces[maxNnbrs*cOwn];
                }

                ccFaces[ccnFaces] = nCoarseFaces;
                initCoarseNeighb[nCoarseFaces] = cNei;
                faceRestrictAddr[fineFacei] = nCoarseFaces;
                ccnFaces++;

                // new coarse face created
                nCoarseFaces++;
            }
        }
    } // end for all fine faces


    // Renumber into upper-triangular order

    // All coarse owner-neighbour storage
    labelList coarseOwner(nCoarseFaces);
    labelList coarseNeighbour(nCoarseFaces);
    labelList coarseFaceMap(nCoarseFaces);

    label coarseFacei = 0;

    forAll(cCellnFaces, cci)
    {
        label* cFaces = &cCellFaces[maxNnbrs*cci];
        label ccnFaces = cCellnFaces[cci];

        for (int i=0; i<ccnFaces; i++)
        {
            coarseOwner[coarseFacei] = cci;
            coarseNeighbour[coarseFacei] = initCoarseNeighb[cFaces[i]];
            coarseFaceMap[cFaces[i]] = coarseFacei;
            coarseFacei++;
        }
    }

    forAll(faceRestrictAddr, fineFacei)
    {
        if (faceRestrictAddr[fineFacei] >= 0)
        {
            faceRestrictAddr[fineFacei] =
                coarseFaceMap[faceRestrictAddr[fineFacei]];
        }
    }


    // Create face-flip status
    faceFlipMap_.set(fineLevelIndex, new boolList(nFineFaces, false));
    boolList& faceFlipMap = faceFlipMap_[fineLevelIndex];


    label nFlipped = 0;
    label nDissapear = 0;

    forAll(faceRestrictAddr, fineFacei)
    {
        label coarseFacei = faceRestrictAddr[fineFacei];

        if (coarseFacei >= 0)
        {
            // Maps to coarse face
            label cOwn = coarseOwner[coarseFacei];
            label cNei = coarseNeighbour[coarseFacei];

            label rmUpperAddr = restrictMap[upperAddr[fineFacei]];
            label rmLowerAddr = restrictMap[lowerAddr[fineFacei]];

            if (cOwn == rmUpperAddr && cNei == rmLowerAddr)
            {
                faceFlipMap[fineFacei] = true;
                nFlipped++;
            }
            else if (cOwn == rmLowerAddr && cNei == rmUpperAddr)
            {
                //faceFlipMap[fineFacei] = false;
            }
            else
            {
                FatalErrorIn("GAMGAgglomeration::agglomerateLduAddressing(..)")
                    << "problem."
                    << " fineFacei:" << fineFacei
                    << " rmUpperAddr:" << rmUpperAddr
                    << " rmLowerAddr:" << rmLowerAddr
                    << " coarseFacei:" << coarseFacei
                    << " cOwn:" << cOwn
                    << " cNei:" << cNei
                    << exit(FatalError);
            }
        }
        else
        {
            nDissapear++;
        }
    }

    Pout<< "Total faces:" << faceRestrictAddr.size()
        << " nDissapear:" << nDissapear
        << " unflipped:" << faceRestrictAddr.size()-nFlipped
        << " flipped:" << nFlipped
        << endl;



    // Clear the temporary storage for the coarse cell data
    cCellnFaces.setSize(0);
    cCellFaces.setSize(0);
    initCoarseNeighb.setSize(0);
    coarseFaceMap.setSize(0);


    // Create coarse-level interfaces

    // Get reference to fine-level interfaces
    const lduInterfacePtrsList& fineInterfaces =
        interfaceLevels_[fineLevelIndex];

    // Create coarse-level interfaces
    primitiveInterfaces_.set
    (
        fineLevelIndex + 1,
        new PtrList<const lduInterface>(fineInterfaces.size())
    );
    PtrList<const lduInterface>& coarsePrimInterfaces =
        primitiveInterfaces_[fineLevelIndex + 1];

    interfaceLevels_.set
    (
        fineLevelIndex + 1,
        new lduInterfacePtrsList(fineInterfaces.size())
    );
    lduInterfacePtrsList& coarseInterfaces =
        interfaceLevels_[fineLevelIndex + 1];

    labelListList coarseInterfaceAddr(fineInterfaces.size());

    nPatchFaces_.set(fineLevelIndex, new labelList(fineInterfaces.size(), 0));
    labelList& nPatchFaces = nPatchFaces_[fineLevelIndex];

    patchFaceRestrictAddressing_.set
    (
        fineLevelIndex,
        new labelListList(fineInterfaces.size())
    );
    labelListList& patchFineToCoarse =
        patchFaceRestrictAddressing_[fineLevelIndex];


    // Initialise transfer of restrict addressing on the interface
    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            fineInterfaces[inti].initInternalFieldTransfer
            (
                Pstream::nonBlocking,
                restrictMap
            );
        }
    }

    if (Pstream::parRun())
    {
        Pstream::waitRequests();
    }


    // Add the coarse level
    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            coarsePrimInterfaces.set
            (
                inti,
                GAMGInterface::New
                (
                    inti,
                    coarseInterfaces,
                    fineInterfaces[inti],
                    fineInterfaces[inti].interfaceInternalField(restrictMap),
                    fineInterfaces[inti].internalFieldTransfer
                    (
                        Pstream::nonBlocking,
                        restrictMap
                    ),
                    fineLevelIndex,
                    fineMesh.comm()
                ).ptr()
            );

            coarseInterfaces.set(inti, &coarsePrimInterfaces[inti]);
            coarseInterfaceAddr[inti] = coarseInterfaces[inti].faceCells();
            nPatchFaces[inti] = coarseInterfaceAddr[inti].size();
            patchFineToCoarse[inti] = refCast<const GAMGInterface>
            (
                coarseInterfaces[inti]
            ).faceRestrictAddressing();
        }
    }


    if (debug)
    {
        Pout<< "GAMGAgglomeration :"
            << " agglomerated level " << fineLevelIndex
            << " from nCells:" << fineMeshAddr.size()
            << " nFaces:" << upperAddr.size()
            << " to nCells:" << nCoarseCells
            << " nFaces:" << nCoarseFaces
            << endl;
    }


    // Set the coarse ldu addressing onto the list
    meshLevels_.set
    (
        fineLevelIndex,
        new lduPrimitiveMesh
        (
            nCoarseCells,
            coarseOwner,
            coarseNeighbour,
            coarseInterfaceAddr,
            coarseInterfaces,
            fineMeshAddr.patchSchedule(),
            fineMesh.comm(),
            true
        )
    );
}


void Foam::GAMGAgglomeration::gatherMeshes
(
    const label meshComm,
    const labelList& procAgglomMap,
    const labelList& procIDs,
    const label allMeshComm,

    const lduMesh& myMesh,

    autoPtr<lduPrimitiveMesh>& allMesh,
    labelList& procCellOffsets,
    labelListList& procFaceMap,
    labelListList& procBoundaryMap,
    labelListListList& procBoundaryFaceMap
)
{
    label oldWarn = UPstream::warnComm;
    UPstream::warnComm = meshComm;

    //Pout<< "GAMGAgglomeration : gathering myMesh (level="
    //    << levelIndex
    //    << ") using communicator " << meshComm << endl;

    PtrList<lduMesh> otherMeshes;
    lduPrimitiveMesh::gather(meshComm, myMesh, procIDs, otherMeshes);

    if (Pstream::myProcNo(meshComm) == procIDs[0])
    {
        // Combine all addressing
        // ~~~~~~~~~~~~~~~~~~~~~~

        //Pout<< "Own Mesh " << myMesh.lduAddr().size() << endl;
        //forAll(otherMeshes, i)
        //{
        //    Pout<< "    otherMesh " << i << " "
        //        << otherMeshes[i].lduAddr().size()
        //        << endl;
        //}

        labelList procFaceOffsets;

        allMesh.reset
        (
            new lduPrimitiveMesh
            (
                allMeshComm,
                procAgglomMap,

                procIDs,
                myMesh,
                otherMeshes,

                procCellOffsets,
                procFaceOffsets,
                procFaceMap,
                procBoundaryMap,
                procBoundaryFaceMap
            )
        );

        //Pout<< "** Agglomerated Mesh " << allMesh().info();
    }

    UPstream::warnComm = oldWarn;
}


void Foam::GAMGAgglomeration::procAgglomerateLduAddressing
(
    const label meshComm,
    const labelList& procAgglomMap,
    const labelList& procIDs,
    const label allMeshComm,

    const label levelIndex
)
{
    Pout<< "** GAMGAgglomeration:procAgglomerateLduAddressing **" << endl;
    Pout<< "level:" << levelIndex
        << " havefinemesh:" << meshLevels_.set(levelIndex-1) << endl;

    const lduMesh& myMesh = meshLevels_[levelIndex-1];

    Pout<< "my mesh                     :" << myMesh.info();
    Pout<< "nCoarseCells                :" << nCells_[levelIndex] << endl;
    Pout<< "restrictAddressing_         :"
        << restrictAddressing_[levelIndex].size()
        << " max:" << max(restrictAddressing_[levelIndex]) << endl;
    Pout<< "faceRestrictAddressing_     :"
        << faceRestrictAddressing_[levelIndex].size()
        << " max:" << max(faceRestrictAddressing_[levelIndex]) << endl;
    Pout<< "patchFaceRestrictAddressing_:" << endl;
    forAll(patchFaceRestrictAddressing_[levelIndex], patchI)
    {
        const labelList& map =
            patchFaceRestrictAddressing_[levelIndex][patchI];

        if (map.size())
        {
            Pout<< "    patch:" << patchI
                << " size:" << map.size()
                << " max:" << max(map)
                << endl;
        }
    }
    Pout<< "interfaceLevels_:" << endl;
    forAll(interfaceLevels_[levelIndex], patchI)
    {
        if (interfaceLevels_[levelIndex].set(patchI))
        {
            Pout<< "    patch:" << patchI
                << " interface:" << interfaceLevels_[levelIndex][patchI].type()
                << " size:"
                << interfaceLevels_[levelIndex][patchI].faceCells().size()
                << endl;
        }
    }




    label oldWarn = UPstream::warnComm;
    UPstream::warnComm = meshComm;


    procAgglomMap_.set(levelIndex, new labelList(procAgglomMap));
    agglomProcIDs_.set(levelIndex, new labelList(procIDs));
    procCommunicator_[levelIndex] = allMeshComm;

    // These could only be set on the master procs but it is
    // quite convenient to also have them on the slaves
    procCellOffsets_.set(levelIndex, new labelList(0));
    procFaceMap_.set(levelIndex, new labelListList(0));
    procBoundaryMap_.set(levelIndex, new labelListList(0));
    procBoundaryFaceMap_.set(levelIndex, new labelListListList(0));

    autoPtr<lduPrimitiveMesh> allMesh;
    gatherMeshes
    (
        meshComm,
        procAgglomMap,
        procIDs,
        allMeshComm,

        myMesh,

        allMesh,
        procCellOffsets_[levelIndex],
        procFaceMap_[levelIndex],
        procBoundaryMap_[levelIndex],
        procBoundaryFaceMap_[levelIndex]
    );


    if (Pstream::myProcNo(meshComm) == procIDs[0])
    {
        // Combine all addressing
        // ~~~~~~~~~~~~~~~~~~~~~~

        meshLevels_.set(levelIndex-1, allMesh);

        const lduPrimitiveMesh& levelMesh = meshLevels_[levelIndex-1];
        const lduInterfacePtrsList& rawInterfaces = levelMesh.rawInterfaces();
        lduInterfacePtrsList& levelInterfaces = interfaceLevels_[levelIndex];
        levelInterfaces.setSize(rawInterfaces.size());
        forAll(rawInterfaces, inti)
        {
            if (rawInterfaces.set(inti))
            {
                levelInterfaces.set(inti, &rawInterfaces[inti]);
            }
            else
            {
                levelInterfaces.set(inti, NULL);
            }
        }
    }

    // Combine restrict addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~


    procAgglomerateRestrictAddressing
    (
        meshComm,
        procIDs,
        levelIndex
    );

    if (Pstream::myProcNo(meshComm) != procIDs[0])
    {
        clearLevel(levelIndex);
    }
    else
    {
        //const lduPrimitiveMesh& levelMesh = meshLevels_[levelIndex-1];
        //Pout<< "** DONE GAMGAgglomeration:procAgglomerateLduAddressing **"
        //    << endl;
        //Pout<< "my mesh:" << levelMesh.info();
        //Pout<< "nCoarseCells:" << nCells_[levelIndex] << endl;
        //Pout<< "restrictAddressing_:"
        //    << restrictAddressing_[levelIndex].size()
        //    << " max:" << max(restrictAddressing_[levelIndex]) << endl;
        //Pout<< "** DONE GAMGAgglomeration:procAgglomerateLduAddressing **"
        //    << endl;
    }

    UPstream::warnComm = oldWarn;
}


void Foam::GAMGAgglomeration::procAgglomerateRestrictAddressing
(
    const label comm,
    const labelList& procIDs,
    const label levelIndex
)
{
    //Pout<< "** GAMGAgglomeration:procAgglomerateRestrictAddressing **"
    //    << endl;
    //Pout<< "level:" << levelIndex << endl;
    //Pout<< "procIDs:" << procIDs << endl;
    //Pout<< "myProcNo:" << UPstream::myProcNo(comm) << endl;
    //Pout<< "procCellOffsets_:" << procCellOffsets_[levelIndex] << endl;

    //const lduMesh& levelMesh = meshLevel(levelIndex);
    //Pout<< "fine:" << restrictAddressing_[levelIndex].size() << endl;

    // Collect number of cells
    labelList nFineCells;
    gatherList
    (
        comm,
        procIDs,
        restrictAddressing_[levelIndex].size(),
        nFineCells
    );
    //const labelList& offsets = procCellOffsets_[levelIndex];
    //Pout<< "offsets:" << offsets << endl;
    //Pout<< "nFineCells:" << nFineCells << endl;

    labelList offsets(nFineCells.size()+1);
    {
        offsets[0] = 0;
        forAll(nFineCells, i)
        {
            offsets[i+1] = offsets[i] + nFineCells[i];
        }
    }
    //Pout<< "offsets:" << offsets << endl;

    // Combine and renumber nCoarseCells
    //Pout<< "Starting nCells_ ..." << endl;
    labelList nCoarseCells;
    gatherList
    (
        comm,
        procIDs,
        nCells_[levelIndex],
        nCoarseCells
    );

    // (cell)restrictAddressing
    //Pout<< "Starting restrictAddressing_ ..." << endl;
    const globalIndex cellOffsetter(offsets);

    labelList procRestrictAddressing;
    cellOffsetter.gather
    (
        comm,
        procIDs,
        restrictAddressing_[levelIndex],
        procRestrictAddressing
    );


    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        labelList coarseCellOffsets(procIDs.size()+1);
        {
            coarseCellOffsets[0] = 0;
            forAll(procIDs, i)
            {
                coarseCellOffsets[i+1] = coarseCellOffsets[i]+nCoarseCells[i];
            }
        }
        //label nOldCoarseCells = nCells_[levelIndex];

        nCells_[levelIndex] = coarseCellOffsets.last();

        //Pout<< "Finished nCoarseCells_ ..."
        //    << " was:" << nOldCoarseCells
        //    << " now:" << nCells_[levelIndex] << endl;


        // Renumber consecutively
        for (label procI = 1; procI < procIDs.size(); procI++)
        {
            SubList<label> procSlot
            (
                procRestrictAddressing,
                offsets[procI+1]-offsets[procI],
                offsets[procI]
            );
            forAll(procSlot, i)
            {
                procSlot[i] += coarseCellOffsets[procI];
            }
        }

        //Pout<< "coarseCellOffsets:" << coarseCellOffsets << endl;
        restrictAddressing_[levelIndex].transfer(procRestrictAddressing);
        //Pout<< "restrictAddressing_:"
        //    << restrictAddressing_[levelIndex].size()
        //    << " min:" << min(restrictAddressing_[levelIndex])
        //    << " max:" << max(restrictAddressing_[levelIndex])
        //    << endl;
        //Pout<< "Finished restrictAddressing_ ..." << endl;
    }
    //Pout<< "** DONE GAMGAgglomeration:procAgglomerateRestrictAddressing **"
    //    << endl;
}


//void Foam::GAMGAgglomeration::gatherList
//(
//    const label comm,
//    const labelList& procIDs,
//
//    const label myVal,
//    labelList& vals,
//    const int tag
//)
//{
//    vals.setSize(procIDs.size());
//
//    if (Pstream::myProcNo(comm) == procIDs[0])
//    {
//        vals[0] = myVal;
//
//        for (label i = 1; i < procIDs.size(); i++)
//        {
//            label& slaveVal = vals[i];
//            IPstream::read
//            (
//                Pstream::scheduled,
//                procIDs[i],
//                reinterpret_cast<char*>(&slaveVal),
//                sizeof(slaveVal),
//                tag,
//                comm
//            );
//        }
//    }
//    else
//    {
//        OPstream::write
//        (
//            Pstream::scheduled,
//            procIDs[0],
//            reinterpret_cast<const char*>(&myVal),
//            sizeof(myVal),
//            tag,
//            comm
//        );
//    }
//}


void Foam::GAMGAgglomeration::calculateRegionMaster
(
    const label comm,
    const labelList& procAgglomMap,
    labelList& masterProcs,
    List<int>& agglomProcIDs
)
{
    // Determine the master processors
    Map<label> agglomToMaster(procAgglomMap.size());

    forAll(procAgglomMap, procI)
    {
        label coarseI = procAgglomMap[procI];

        Map<label>::iterator fnd = agglomToMaster.find(coarseI);
        if (fnd == agglomToMaster.end())
        {
            agglomToMaster.insert(coarseI, procI);
        }
        else
        {
            fnd() = min(fnd(), procI);
        }
    }

    masterProcs.setSize(agglomToMaster.size());
    forAllConstIter(Map<label>, agglomToMaster, iter)
    {
        masterProcs[iter.key()] = iter();
    }


    // Collect all the processors in my agglomeration
    label myProcID = Pstream::myProcNo(comm);
    label myAgglom = procAgglomMap[myProcID];

    // Get all processors agglomerating to the same coarse
    // processor
    agglomProcIDs = findIndices(procAgglomMap, myAgglom);
    // Make sure the master is the first element.
    label index = findIndex
    (
        agglomProcIDs,
        agglomToMaster[myAgglom]
    );
    Swap(agglomProcIDs[0], agglomProcIDs[index]);
}



// ************************************************************************* //
