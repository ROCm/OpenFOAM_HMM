/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "globalIndex.H"
#include "lduPrimitiveMesh.H"
#include "processorGAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkUpperTriangular
(
    const label size,
    const labelUList& l,
    const labelUList& u
)
{
    forAll(l, faceI)
    {
        if (u[faceI] < l[faceI])
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Reversed face. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI] << abort(FatalError);
        }
        if (l[faceI] < 0 || u[faceI] < 0 || u[faceI] >= size)
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Illegal cell label. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI] << abort(FatalError);
        }
    }

    for (label faceI=1; faceI < l.size(); faceI++)
    {
        if (l[faceI-1] > l[faceI])
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Lower not in incremental cell order."
                << " Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << " previous l:" << l[faceI-1] << abort(FatalError);
        }
        else if (l[faceI-1] == l[faceI])
        {
            // Same cell.
            if (u[faceI-1] > u[faceI])
            {
                FatalErrorIn
                (
                    "checkUpperTriangular"
                    "(const label, const labelUList&, const labelUList&)"
                )   << "Upper not in incremental cell order."
                    << " Problem at face " << faceI
                    << " l:" << l[faceI] << " u:" << u[faceI]
                    << " previous u:" << u[faceI-1] << abort(FatalError);
            }
        }
    }
}


void print(const string& msg, const lduMesh& mesh)
{
    const lduAddressing& addressing = mesh.lduAddr();
    const lduInterfacePtrsList interfaces = mesh.interfaces();

    Pout<< "Mesh:" << msg.c_str() << nl
        << "    cells:" << addressing.size() << nl
        << "    faces:" << addressing.lowerAddr().size() << nl
        << "    patches:" << interfaces.size() << nl;


    const labelUList& l = addressing.lowerAddr();
    const labelUList& startAddr = addressing.losortStartAddr();
    const labelUList& addr = addressing.losortAddr();

    forAll(addressing, cellI)
    {
        Pout<< "    cell:" << cellI << nl;

        label start = startAddr[cellI];
        label end = startAddr[cellI+1];

        for (label index = start; index < end; index++)
        {
            Pout<< "        nbr:" << l[addr[index]] << nl;
        }
    }

    Pout<< "Patches:" << nl;
    forAll(interfaces, i)
    {
        if (interfaces.set(i))
        {
            if (isA<processorLduInterface>(interfaces[i]))
            {
                const processorLduInterface& pldui =
                    refCast<const processorLduInterface>(interfaces[i]);
                Pout<< "    " << i
                    << " me:" << pldui.myProcNo()
                    << " nbr:" << pldui.neighbProcNo()
                    << " comm:" << pldui.comm()
                    << " tag:" << pldui.tag()
                    << nl;
            }

            {
                Pout<< "    " << i << " addressing:" << nl;
                const labelUList& faceCells = interfaces[i].faceCells();
                forAll(faceCells, i)
                {
                    Pout<< "\t\t" << i << '\t' << faceCells[i] << nl;
                }
            }
        }
    }
}



template<class ProcPatch>
lduSchedule nonBlockingSchedule
(
    const lduInterfacePtrsList& interfaces
)
{
    lduSchedule schedule(2*interfaces.size());
    label slotI = 0;

    forAll(interfaces, i)
    {
        if (interfaces.set(i) && !isA<ProcPatch>(interfaces[i]))
        {
            schedule[slotI].patch = i;
            schedule[slotI].init = true;
            slotI++;
            schedule[slotI].patch = i;
            schedule[slotI].init = false;
            slotI++;
        }
    }

    forAll(interfaces, i)
    {
        if (interfaces.set(i) && isA<ProcPatch>(interfaces[i]))
        {
            schedule[slotI].patch = i;
            schedule[slotI].init = true;
            slotI++;
        }
    }

    forAll(interfaces, i)
    {
        if (interfaces.set(i) && isA<ProcPatch>(interfaces[i]))
        {
            schedule[slotI].patch = i;
            schedule[slotI].init = false;
            slotI++;
        }
    }

    return schedule;
}



// Say combining procs 1,13,9,4.
// - cells get added in order 1,4,9,13
// - internal faces get added in order of processor and then order of
//   neighbouring processor:
//     - internal faces of 1
//     - faces between 1 and 4 (keep orientation)
//     - faces between 1 and 9
//     - faces between 1 and 13
//     - interfaces between 1 and other processors
//     and then
//     - internal faces of 4
//     - faces between 4 and 1 (skip; already added)
//     - faces between 4 and 9
//     - faces between 4 and 13
//     - interfaces between 4 and other processors
//     etc.
// - this still does not guarantee that upper-triangular order is preserved:
// A cell has two processor faces. These now become internal faces and if
// the cell numbering on the other side is not in the same order as the internal
// faces we have upper-triangular.
autoPtr<lduPrimitiveMesh> combineMeshes
(
    const label newComm,
    const labelList& procIDs,
    const PtrList<lduMesh>& procMeshes,

    labelList& cellOffsets,     // per mesh the starting cell
    labelList& faceOffsets,     // per mesh the starting face
    labelList& interfaceOffsets // per mesh,per interface the starting face
)
{
    // Sanity check.
    for (label i = 1; i < procIDs.size(); i++)
    {
        if (procIDs[i] <= procIDs[i-1])
        {
            FatalErrorIn
            (
                "combineMeshes(const labelList&, const PtrList<lduMesh>&)"
            )   << "Processor " << procIDs[i] << " at index " << i
                << " should be higher numbered than its predecessor "
                << procIDs[i-1]
                << exit(FatalError);
        }
    }


    label currentComm = procMeshes[0].comm();


    // Cells get added in order.
    cellOffsets.setSize(procMeshes.size()+1);
    cellOffsets[0] = 0;
    forAll(procMeshes, procMeshI)
    {
        const lduMesh& mesh = procMeshes[procMeshI];
        cellOffsets[procMeshI+1] =
            cellOffsets[procMeshI]
          + mesh.lduAddr().size();
    }

    // Faces initially get added in order, sorted later
    labelList internalFaceOffsets(procMeshes.size()+1);
    internalFaceOffsets[0] = 0;
    forAll(procMeshes, procMeshI)
    {
        const lduMesh& mesh = procMeshes[procMeshI];
        internalFaceOffsets[procMeshI+1] =
            internalFaceOffsets[procMeshI]
          + mesh.lduAddr().lowerAddr().size();
    }

    // Count how faces get added. Interfaces inbetween get merged.

    // Merged interfaces: map from two processors back to
    // - procMeshes
    // - interface in procMesh
    EdgeMap<labelPairList> mergedMap
    (
        2*procMeshes[0].interfaces().size()
    );


    // Unmerged interfaces: map from two processors back to
    // - procMeshes
    // - interface in procMesh
    EdgeMap<labelPairList> unmergedMap
    (
        2*procMeshes[0].interfaces().size()
    );
    // Per interface the size
    //DynamicList<label> interfaceSize(2*procMeshes[0].interfaces().size());

    label nBoundaryFaces = 0;
    label nInterfaces = 0;
    labelList nCoupledFaces(procMeshes.size(), 0);

    forAll(procMeshes, procMeshI)
    {
        const lduInterfacePtrsList interfaces =
            procMeshes[procMeshI].interfaces();

        // Get sorted order of processors
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                if (isA<processorGAMGInterface>(interfaces[intI]))
                {
                    const processorGAMGInterface& pldui =
                        refCast<const processorGAMGInterface>
                        (
                            interfaces[intI]
                        );
                    if (pldui.myProcNo() != procIDs[procMeshI])
                    {
                        FatalErrorIn("combineMeshes()")
                            << "proc:" << procIDs[procMeshI]
                            << " myProcNo:" << pldui.myProcNo()
                            << abort(FatalError);
                    }


                    const edge procEdge
                    (
                        min(pldui.myProcNo(), pldui.neighbProcNo()),
                        max(pldui.myProcNo(), pldui.neighbProcNo())
                    );


                    label index = findIndex(procIDs, pldui.neighbProcNo());

                    if (index == -1)
                    {
                        // Still external interface
                        Pout<< "external interface: myProcNo:"
                            << pldui.myProcNo()
                            << " nbr:" << pldui.neighbProcNo()
                            << " size:" << pldui.faceCells().size()
                            << endl;

                        nBoundaryFaces += pldui.faceCells().size();

                        EdgeMap<labelPairList>::iterator iter =
                            unmergedMap.find(procEdge);

                        if (iter != unmergedMap.end())
                        {
                            iter().append(labelPair(procMeshI, intI));
                        }
                        else
                        {
                            unmergedMap.insert
                            (
                                procEdge,
                                labelPairList(1, labelPair(procMeshI, intI))
                            );
                        }

                        nInterfaces++;
                    }
                    else //if (pldui.myProcNo() < pldui.neighbProcNo())
                    {
                        // Merged interface
                        Pout<< "merged interface: myProcNo:" << pldui.myProcNo()
                            << " nbr:" << pldui.neighbProcNo()
                            << " size:" << pldui.faceCells().size()
                            << endl;
                        if (pldui.myProcNo() < pldui.neighbProcNo())
                        {
                            nCoupledFaces[procMeshI] +=
                                pldui.faceCells().size();
                        }

                        EdgeMap<labelPairList>::iterator iter =
                            mergedMap.find(procEdge);

                        if (iter != mergedMap.end())
                        {
                            iter().append(labelPair(procMeshI, intI));
                        }
                        else
                        {
                            mergedMap.insert
                            (
                                procEdge,
                                labelPairList(1, labelPair(procMeshI, intI))
                            );
                        }
                    }
                }
                else
                {
                    // Still external (non proc) interface
                    FatalErrorIn("combineMeshes()")
                        << "At mesh from processor " << procIDs[procMeshI]
                        << " have interface " << intI
                        << " of unhandled type " << interfaces[intI].type()
                        << exit(FatalError);

                    nBoundaryFaces += interfaces[intI].faceCells().size();
                    nInterfaces++;
                }
            }
        }
    }



    {
        Pout<< "Remaining interfaces:" << endl;
        forAllConstIter(EdgeMap<labelPairList>, unmergedMap, iter)
        {
            Pout<< "    procEdge:" << iter.key() << endl;
            const labelPairList& elems = iter();
            forAll(elems, i)
            {
                label procMeshI = elems[i][0];
                label interfaceI = elems[i][1];
                const lduInterfacePtrsList interfaces =
                    procMeshes[procMeshI].interfaces();
                const processorGAMGInterface& pldui =
                    refCast<const processorGAMGInterface>
                    (
                        interfaces[interfaceI]
                    );

                Pout<< "        proc:" << procIDs[procMeshI]
                    << " interfaceI:" << interfaceI
                    << " between:" << pldui.myProcNo()
                    << " and:" << pldui.neighbProcNo()
                    << endl;
            }
            Pout<< endl;
        }
    }
    {
        Pout<< "Merged interfaces:" << endl;
        forAllConstIter(EdgeMap<labelPairList>, mergedMap, iter)
        {
            Pout<< "    procEdge:" << iter.key() << endl;
            const labelPairList& elems = iter();
            forAll(elems, i)
            {
                label procMeshI = elems[i][0];
                label interfaceI = elems[i][1];
                const lduInterfacePtrsList interfaces =
                    procMeshes[procMeshI].interfaces();
                const processorGAMGInterface& pldui =
                    refCast<const processorGAMGInterface>
                    (
                        interfaces[interfaceI]
                    );

                Pout<< "        proc:" << procIDs[procMeshI]
                    << " interfaceI:" << interfaceI
                    << " between:" << pldui.myProcNo()
                    << " and:" << pldui.neighbProcNo()
                    << endl;
            }
            Pout<< endl;
        }
    }


    // Adapt faceOffsets for internal interfaces
    faceOffsets.setSize(procMeshes.size()+1);
    faceOffsets[0] = 0;
    forAll(procMeshes, procMeshI)
    {
        faceOffsets[procMeshI+1] =
            faceOffsets[procMeshI]
          + procMeshes[procMeshI].lduAddr().lowerAddr().size() //internal
          + nCoupledFaces[procMeshI];   // resolved coupled faces
    }


    // Combine upper and lower
    labelList allLower(faceOffsets.last(), -1);
    labelList allUpper(allLower.size(), -1);


    // Old internal faces and resolved coupled interfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(procMeshes, procMeshI)
    {
        const labelUList& l = procMeshes[procMeshI].lduAddr().lowerAddr();
        const labelUList& u = procMeshes[procMeshI].lduAddr().upperAddr();

        // Add internal faces
        label allFaceI = faceOffsets[procMeshI];

        forAll(l, faceI)
        {
            allLower[allFaceI] = cellOffsets[procMeshI]+l[faceI];
            allUpper[allFaceI] = cellOffsets[procMeshI]+u[faceI];
            allFaceI++;
        }

        // Add converted interfaces
        const lduInterfacePtrsList interfaces =
            procMeshes[procMeshI].interfaces();

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                if (isA<processorGAMGInterface>(interfaces[intI]))
                {
                    const processorGAMGInterface& pldui =
                        refCast<const processorGAMGInterface>
                        (
                            interfaces[intI]
                        );

                    // Look up corresponding interfaces
                    label myP = pldui.myProcNo();
                    label nbrP = pldui.neighbProcNo();

                    if (myP < nbrP)
                    {
                        EdgeMap<labelPairList>::const_iterator fnd =
                            mergedMap.find(edge(myP, nbrP));

                        if (fnd != mergedMap.end())
                        {
                            const labelPairList& elems = fnd();

                            // Find nbrP in elems
                            label nbrProcMeshI = -1;
                            label nbrIntI = -1;
                            if (procIDs[elems[0][0]] == nbrP)
                            {
                                nbrProcMeshI = elems[0][0];
                                nbrIntI = elems[0][1];
                            }
                            else
                            {
                                nbrProcMeshI = elems[1][0];
                                nbrIntI = elems[1][1];
                            }

                            if
                            (
                                elems.size() != 2
                             || procIDs[nbrProcMeshI] != nbrP
                            )
                            {
                                FatalErrorIn("combineMeshes()")
                                    << "elems:" << elems << abort(FatalError);
                            }


                            const lduInterfacePtrsList nbrInterfaces =
                                procMeshes[nbrProcMeshI].interfaces();

                            const processorGAMGInterface& nbrPldui =
                                refCast<const processorGAMGInterface>
                                (
                                    nbrInterfaces[nbrIntI]
                                );
                            const labelUList& faceCells = pldui.faceCells();
                            const labelUList& nbrFaceCells =
                                nbrPldui.faceCells();

                            forAll(faceCells, pfI)
                            {
                                allLower[allFaceI] =
                                    cellOffsets[procMeshI]+faceCells[pfI];
                                allUpper[allFaceI] =
                                    cellOffsets[nbrProcMeshI]+nbrFaceCells[pfI];
                                allFaceI++;
                            }
                        }
                    }
                }
            }
        }
    }


    // Kept interfaces
    // ~~~~~~~~~~~~~~~

    lduInterfacePtrsList coarseInterfaces(nInterfaces);
    label allInterfaceI = 0;

    forAllConstIter(EdgeMap<labelPairList>, unmergedMap, iter)
    {
        Pout<< "procEdge:" << iter.key() << endl;
        const labelPairList& elems = iter();

        // Count
        label n = 0;

        forAll(elems, i)
        {
            label procMeshI = elems[i][0];
            label interfaceI = elems[i][1];
            const lduInterfacePtrsList interfaces =
                procMeshes[procMeshI].interfaces();
            n += interfaces[interfaceI].faceCells().size();
        }

        labelField allFaceCells(n);
        labelField allFaceRestrictAddressing(n);
        n = 0;

        forAll(elems, i)
        {
            label procMeshI = elems[i][0];
            label interfaceI = elems[i][1];
            const lduInterfacePtrsList interfaces =
                procMeshes[procMeshI].interfaces();

            const labelUList& l = interfaces[interfaceI].faceCells();
            forAll(l, faceI)
            {
                allFaceCells[n] = cellOffsets[procMeshI]+l[faceI];
                allFaceRestrictAddressing[n] = n;
                n++;
            }
        }


        // Find out local and remote processor in new communicator
        label minProcI = UPstream::procNo
        (
            newComm,
            currentComm,
            iter.key()[0]
        );
        label maxProcI = UPstream::procNo
        (
            newComm,
            currentComm,
            iter.key()[1]
        );

        coarseInterfaces.set
        (
            allInterfaceI,
            new processorGAMGInterface
            (
                allInterfaceI,
                coarseInterfaces,
                allFaceCells,
                allFaceRestrictAddressing,
                newComm,
                minProcI,
                maxProcI,
                tensorField(),          // forwardT
                Pstream::msgType()      // tag
            )
        );
        allInterfaceI++;
    }


    // Extract faceCells from coarseInterfaces
    labelListList coarseInterfaceAddr(coarseInterfaces.size());
    forAll(coarseInterfaces, coarseIntI)
    {
        if (coarseInterfaces.set(coarseIntI))
        {
            coarseInterfaceAddr[coarseIntI] =
                coarseInterfaces[coarseIntI].faceCells();
        }
    }


    Pout<< "allLower:" << allLower << endl;
    Pout<< "allUpper:" << allUpper << endl;
    checkUpperTriangular(cellOffsets.last(), allLower, allUpper);


    autoPtr<lduPrimitiveMesh> meshPtr
    (
        new lduPrimitiveMesh
        (
            cellOffsets.last(),         //nCells
            allLower,                   //lower
            allUpper,                   //upper
            coarseInterfaceAddr,        //faceCells
            coarseInterfaces,           //interfaces
            nonBlockingSchedule<processorGAMGInterface>(coarseInterfaces),
            newComm,                    //communicator
            true                        //reuse
        )
    );
    return meshPtr;
}



void sendReceive
(
    const label comm,
    const label tag,
    const globalIndex& offsets,
    const scalarField& field,

    scalarField& allField
)
{
    label nProcs = Pstream::nProcs(comm);

    if (Pstream::master(comm))
    {
        allField.setSize(offsets.size());

        // Assign master slot
        SubList<scalar>
        (
            allField,
            offsets.localSize(0),
            offsets.offset(0)
        ).assign(field);

        // Assign slave slots
        for (label procI = 1; procI < nProcs; procI++)
        {
            SubList<scalar> procSlot
            (
                allField,
                offsets.localSize(procI),
                offsets.offset(procI)
            );

            Pout<< "Receiving allField from " << procI
                << " at offset:" << offsets.offset(procI)
                << " size:" << offsets.size()
                << endl;

            IPstream::read
            (
                Pstream::nonBlocking,
                procI,
                reinterpret_cast<char*>(procSlot.begin()),
                procSlot.byteSize(),
                tag,
                comm
            );
        }
    }
    else
    {
        OPstream::write
        (
            Pstream::nonBlocking,
            0,                          // master
            reinterpret_cast<const char*>(field.begin()),
            field.byteSize(),
            tag,
            comm
        );
    }
}


void sendReceive
(
    const label comm,
    const label tag,
    const globalIndex& offsets,
    const FieldField<Field, scalar>& field,

    FieldField<Field, scalar>& allField
)
{
    PstreamBuffers pBufs(Pstream::nonBlocking, Pstream::msgType(), comm);

    if (!Pstream::master(comm))
    {
        UOPstream toMaster(Pstream::masterNo(), pBufs);

        Pout<< "To 0 sending " << field.size()
            << " fields." << endl;

        forAll(field, intI)
        {
            toMaster << field[intI];
        }
    }
    pBufs.finishedSends();
    if (Pstream::master(comm))
    {
        allField.setSize(offsets.size());
        forAll(allField, i)
        {
            allField.set(i, new scalarField(0));
        }

        // Insert master values
        forAll(field, intI)
        {
            allField[intI] = field[intI];
        }


        // Receive and insert slave values
        label nProcs = Pstream::nProcs(comm);

        for (label procI = 1; procI < nProcs; procI++)
        {
            UIPstream fromSlave(procI, pBufs);

            label nSlaveInts = offsets.localSize(procI);

            Pout<< "From " << procI << " receiving "
                << nSlaveInts << " fields." << endl;

            for (label intI = 0; intI < nSlaveInts; intI++)
            {
                label slotI = offsets.toGlobal(procI, intI);

                Pout<< "    int:" << intI << " goes into slot " << slotI
                    << endl;

                fromSlave >> allField[slotI];
            }
        }
    }
}



void collect
(
    const label comm,
    const globalIndex& cellOffsets,
    const globalIndex& faceOffsets,

    const scalarField& diagonal,
    const scalarField& upper,
    const scalarField& lower,

    scalarField& allDiagonal,
    scalarField& allUpper,
    scalarField& allLower
)
{
    label nOutstanding = Pstream::nRequests();
    int allDiagonalTag = Pstream::allocateTag("allDiagonal:" __FILE__);
    int allUpperTag = Pstream::allocateTag("allUpper:" __FILE__);
    int allLowerTag = Pstream::allocateTag("allLower:" __FILE__);


    sendReceive
    (
        comm,
        allDiagonalTag,
        cellOffsets,
        diagonal,
        allDiagonal
    );

    sendReceive
    (
        comm,
        allUpperTag,
        faceOffsets,
        upper,
        allUpper
    );

    sendReceive
    (
        comm,
        allLowerTag,
        faceOffsets,
        lower,
        allLower
    );

    Pstream::waitRequests(nOutstanding);

    Pstream::freeTag("allDiagonal:" __FILE__, allDiagonalTag);
    Pstream::freeTag("allUpper:" __FILE__, allUpperTag);
    Pstream::freeTag("allLower:" __FILE__, allLowerTag);
}


void setCommunicator(fvMesh& mesh, const label newComm)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // The current state is consistent with the mesh so check where the new
    // communicator is and adjust accordingly.

    forAll(pbm, patchI)
    {
        if (isA<processorPolyPatch>(pbm[patchI]))
        {
            processorPolyPatch& ppp = const_cast<processorPolyPatch&>
            (
                refCast
                <
                    const processorPolyPatch
                >(pbm[patchI])
            );

            label thisRank = UPstream::procNo
            (
                newComm,
                ppp.comm(),
                ppp.myProcNo()
            );
            label nbrRank = UPstream::procNo
            (
                newComm,
                ppp.comm(),
                ppp.neighbProcNo()
            );

            //ppp.comm() = newComm;
            ppp.myProcNo() = thisRank;
            ppp.neighbProcNo() = nbrRank;
        }
    }
    mesh.polyMesh::comm() = newComm;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    simpleControl simple(mesh);

    //const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;


    // Get a subset of processors
    labelList subProcs(3);
    subProcs[0] = 0;
    subProcs[1] = 1;
    subProcs[2] = 2;


    const UPstream::communicator newComm
    (
        UPstream::worldComm,
        subProcs,
        true
    );


    Pout<< "procIDs world  :" << UPstream::procID(UPstream::worldComm) << endl;
    Pout<< "procIDs newComm:" << UPstream::procID(newComm) << endl;


// On ALL processors: get the interfaces:
lduInterfacePtrsList interfaces(mesh.interfaces());
const lduAddressing& addressing = mesh.lduAddr();

if (Pstream::myProcNo(newComm) != -1)
{
    print("InitialMesh", mesh);

    labelList procIDs(2);
    procIDs[0] = 0;
    procIDs[1] = 1;

    DynamicList<label> compactMap(interfaces.size());
    forAll(interfaces, intI)
    {
        if (interfaces.set(intI))
        {
            compactMap.append(intI);
        }
    }

    // Collect meshes from procs 0,1 (in newComm) on 1.
    if (Pstream::myProcNo(newComm) == procIDs[0])
    {
        PtrList<lduMesh> procMeshes(procIDs.size());

        // Master mesh (copy for now. TBD)
        {
            // Convert to processorGAMGInterfaces in new communicator
            lduInterfacePtrsList compactInterfaces(compactMap.size());
            labelListList compactPatchAddr(compactMap.size());

            forAll(compactMap, compactI)
            {
                label intI = compactMap[compactI];

                const processorLduInterface& pldui =
                    refCast<const processorLduInterface>(interfaces[intI]);

                compactPatchAddr[compactI] = interfaces[intI].faceCells();
                labelList faceRestrictAddressing
                (
                    identity(interfaces[intI].faceCells().size())
                );


                Pout<< "for interface:" << intI
                    << " comm:" << pldui.comm()
                    << " myProcNo:" << pldui.myProcNo()
                    << " in newComm:" << UPstream::procNo
                        (
                            newComm,
                            pldui.comm(),
                            pldui.myProcNo()
                        )
                    << endl;
                Pout<< "for interface:" << intI
                    << " comm:" << pldui.comm()
                    << " neighbProcNo:" << pldui.neighbProcNo()
                    << " in newComm:" << UPstream::procNo
                        (
                            newComm,
                            pldui.comm(),
                            pldui.neighbProcNo()
                        )
                    << endl;


                compactInterfaces.set
                (
                    compactI,
                    new processorGAMGInterface
                    (
                        compactI,
                        compactInterfaces,
                        compactPatchAddr[compactI],
                        faceRestrictAddressing,

                        newComm,                        //pldui.comm(),
                        UPstream::procNo
                        (
                            newComm,
                            pldui.comm(),
                            pldui.myProcNo()
                        ),
                        UPstream::procNo
                        (
                            newComm,
                            pldui.comm(),
                            pldui.neighbProcNo()
                        ),
                        pldui.forwardT(),
                        pldui.tag()
                    )
                );
            }



            procMeshes.set
            (
                0,
                new lduPrimitiveMesh
                (
                    addressing.size(),
                    addressing.lowerAddr(),
                    addressing.upperAddr(),
                    compactPatchAddr,
                    compactInterfaces,
                    nonBlockingSchedule<processorGAMGInterface>
                    (
                        compactInterfaces
                    ),  //patch schedule,
                    newComm
                )
            );
        }

        // Slave meshes
        for (label i = 1; i < procIDs.size(); i++)
        {
            //Pout<< "on master :"
            //    << " receiving from slave " << procIDs[i]
            //    << endl;

            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,          // bufSize
                Pstream::msgType(),
                newComm
            );

            label nCells = readLabel(fromSlave);
            labelList lowerAddr(fromSlave);
            labelList upperAddr(fromSlave);
            labelListList patchAddr(fromSlave);
            labelList comm(fromSlave);
            labelList myProcNo(fromSlave);
            labelList neighbProcNo(fromSlave);
            List<tensorField> forwardT(fromSlave);
            labelList tag(fromSlave);

            // Convert to processorGAMGInterfaces
            lduInterfacePtrsList compactInterfaces(patchAddr.size());
            forAll(patchAddr, compactI)
            {
                labelList faceRestrictAddressing
                (
                    identity(patchAddr[compactI].size())
                );

                compactInterfaces.set
                (
                    compactI,
                    new processorGAMGInterface
                    (
                        compactI,
                        compactInterfaces,
                        patchAddr[compactI],
                        faceRestrictAddressing,

                        comm[compactI],
                        UPstream::procNo
                        (
                            newComm,
                            comm[compactI],
                            myProcNo[compactI]
                        ),
                        UPstream::procNo
                        (
                            newComm,
                            comm[compactI],
                            neighbProcNo[compactI]
                        ),

                        forwardT[compactI],
                        tag[compactI]
                    )
                );
            }


            procMeshes.set
            (
                i,
                new lduPrimitiveMesh
                (
                    nCells,
                    lowerAddr,
                    upperAddr,
                    patchAddr,
                    compactInterfaces,
                    nonBlockingSchedule<processorGAMGInterface>
                    (
                        compactInterfaces
                    ),
                    newComm
                )
            );
        }


        // Print a bit
        forAll(procMeshes, i)
        {
            const lduMesh& pMesh = procMeshes[i];
            print("procMesh" + Foam::name(i), pMesh);

            const lduAddressing& addr = pMesh.lduAddr();
            checkUpperTriangular
            (
                addr.size(),
                addr.lowerAddr(),
                addr.upperAddr()
            );
        }


        // Combine
        labelList cellOffsets;
        labelList faceOffsets;
        labelList interfaceOffsets;
        autoPtr<lduPrimitiveMesh> allMeshPtr = combineMeshes
        (
            newComm,
            procIDs,
            procMeshes,

            cellOffsets,     // per mesh the starting cell
            faceOffsets,     // per mesh the starting face
            interfaceOffsets // per mesh,per interface the starting face
        );

        print("ALLMESH", allMeshPtr());
    }
    else if (findIndex(procIDs, Pstream::myProcNo(newComm)) != -1)
    {
        // Send to master

        //- Extract info from processorGAMGInterface
        labelListList patchAddr(compactMap.size());
        //labelListList faceRestrictAddr(compactMap.size());
        labelList comm(compactMap.size());
        labelList myProcNo(compactMap.size());
        labelList neighbProcNo(compactMap.size());
        List<tensorField> forwardT(compactMap.size());
        labelList tag(compactMap.size());

        forAll(compactMap, compactI)
        {
            label intI = compactMap[compactI];

            const processorLduInterface& pldui =
                refCast<const processorLduInterface>(interfaces[intI]);

            patchAddr[compactI] = interfaces[intI].faceCells();
            //faceRestrictAddr[compactI] = pldui.faceRestrictAddresssing();
            comm[compactI] = pldui.comm();
            myProcNo[compactI] = pldui.myProcNo();
            neighbProcNo[compactI] =  pldui.neighbProcNo();
            forwardT[compactI] =  pldui.forwardT();
            tag[compactI] =  pldui.tag();
        }

        //Pout<< "on slave sending to " << Pstream::masterNo() << endl;
        //Pout<< "on slave nCells:" << addressing.size() << endl;
        //Pout<< "on slave lowerAdd:" << addressing.lowerAddr().size() << endl;
        //Pout<< "on slave upperAddr:" << addressing.upperAddr().size() << endl;
        //Pout<< "on slave patchAddr:" << patchAddr.size() << endl;
        //Pout<< "on slave comm:" << comm << endl;
        //Pout<< "on slave myProcNo:" << myProcNo << endl;
        //Pout<< "on slave neighbProcNo:" << neighbProcNo << endl;
        //Pout<< "on slave forwardT:" << forwardT << endl;
        //Pout<< "on slave tag:" << tag << endl;

        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            Pstream::msgType(),
            newComm
        );
        toMaster
            << addressing.size()
            << addressing.lowerAddr()
            << addressing.upperAddr()
            << patchAddr
            //<< faceRestrictAddr
            << comm
            << myProcNo
            << neighbProcNo
            << forwardT
            << tag;
    }
}
return 0;



//if (Pstream::myProcNo(newComm) != -1)
//{
//    scalarField diagonal(11+Pstream::myProcNo(), 11+Pstream::myProcNo());
//    scalarField upper(22+Pstream::myProcNo(), 22+Pstream::myProcNo());
//    scalarField lower(upper.size(), 22+Pstream::myProcNo());
//
//    // Collect the local sizes
//    const globalIndex cellOffsets
//    (
//        diagonal.size(),
//        Pstream::msgType(),
//        newComm,
//        true
//    );
//
//    const globalIndex faceOffsets
//    (
//        upper.size(),
//        Pstream::msgType(),
//        newComm,
//        true
//    );
//
//
//    scalarField allDiagonal;
//    scalarField allUpper;
//    scalarField allLower;
//
//    collect
//    (
//        newComm,
//        cellOffsets,
//        faceOffsets,
//
//        diagonal,
//        upper,
//        lower,
//
//        allDiagonal,
//        allUpper,
//        allLower
//    );
//
//    Pout<< "diagonal:" << diagonal << endl;
//    Pout<< "allDiagonal:" << allDiagonal << endl;
//
//    Pout<< "upper:" << upper << endl;
//    Pout<< "allUpper:" << allUpper << endl;
//
//    Pout<< "lower:" << lower << endl;
//    Pout<< "allLower:" << allLower << endl;
//
//}



//if (Pstream::myProcNo(newComm) != -1)
//{
//    // Interface coefficients. The destination number of interfaces can differ
//    // We have pre-sorted the overall match so we can just slot them in.
//    FieldField<Field, scalar> interfaceBouCoeffs;
//    //FieldField<Field, scalar> interfaceIntCoeffs;
//    {
//        label nInt = 3+Pstream::myProcNo();
//        interfaceBouCoeffs.setSize(nInt);
//        forAll(interfaceBouCoeffs, intI)
//        {
//            interfaceBouCoeffs.set(intI, new scalarField(intI+1));
//            scalarField& fld = interfaceBouCoeffs[intI];
//            fld = Pstream::myProcNo();
//        }
//    }
//
//    const globalIndex interfaceBouOffsets
//    (
//        interfaceBouCoeffs.size(),
//        Pstream::msgType(),
//        newComm,
//        true
//    );
//
//    FieldField<Field, scalar> allInterfaceBouCoeffs;
//    sendReceive
//    (
//        newComm,
//        Pstream::msgType(),
//        interfaceBouOffsets,
//        interfaceBouCoeffs,
//
//        allInterfaceBouCoeffs
//    );
//
//    Pout<< "interfaceBouCoeffs   :" << interfaceBouCoeffs << endl;
//    Pout<< "allInterfaceBouCoeffs:" << allInterfaceBouCoeffs << endl;
//
//}

    {
        Pout<< "World:" << UPstream::worldComm
            << " procID:" << 2
            << " subComm:" << newComm
            << " rank1:" << UPstream::procNo(newComm, UPstream::worldComm, 1)
            << " rank2:" << UPstream::procNo(newComm, UPstream::worldComm, 2)
            << endl;
    }
    //
    //
    //{
    //    setCommunicator(mesh, newComm);
    //    Pout<< "Mesh comm:" << mesh.comm() << endl;
    //    forAll(pbm, patchI)
    //    {
    //        if (isA<processorPolyPatch>(pbm[patchI]))
    //        {
    //            const processorPolyPatch& ppp =
    //            refCast
    //            <
    //                const processorPolyPatch
    //            >(pbm[patchI]);
    //
    //
    //            Pout<< "    " << ppp.name()
    //                << " myRank:" << ppp.myProcNo()
    //                << " nbrRank:" << ppp.neighbProcNo()
    //                << endl;
    //        }
    //    }
    //    setCommunicator(mesh, UPstream::worldComm);
    //}


    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix Teqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
            );


            {
                label oldWarn = UPstream::warnComm;
                UPstream::warnComm = newComm;

                label oldComm = mesh.comm();
                setCommunicator(mesh, newComm);
                Pout<< "** oldcomm:" << oldComm
                    << "  newComm:" << mesh.comm() << endl;

                if (Pstream::myProcNo(mesh.comm()) != -1)
                {
                    solve(Teqn);
                }

                setCommunicator(mesh, oldComm);
                Pout<< "** reset mesh to:" << mesh.comm() << endl;

                UPstream::warnComm = oldWarn;
            }
        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
