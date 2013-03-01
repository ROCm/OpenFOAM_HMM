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

\*---------------------------------------------------------------------------*/

#include "lduPrimitiveMesh.H"
#include "processorLduInterface.H"
#include "EdgeMap.H"
#include "labelPair.H"
#include "processorGAMGInterface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduPrimitiveMesh::checkUpperTriangular
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
            //FatalErrorIn
            WarningIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Reversed face. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                //<< abort(FatalError);
                << endl;
        }
        if (l[faceI] < 0 || u[faceI] < 0 || u[faceI] >= size)
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Illegal cell label. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << abort(FatalError);
        }
    }

    for (label faceI=1; faceI < l.size(); faceI++)
    {
        if (l[faceI-1] > l[faceI])
        {
            //FatalErrorIn
            WarningIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Lower not in incremental cell order."
                << " Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << " previous l:" << l[faceI-1]
                //<< abort(FatalError);
                << endl;
        }
        else if (l[faceI-1] == l[faceI])
        {
            // Same cell.
            if (u[faceI-1] > u[faceI])
            {
                //FatalErrorIn
                WarningIn
                (
                    "checkUpperTriangular"
                    "(const label, const labelUList&, const labelUList&)"
                )   << "Upper not in incremental cell order."
                    << " Problem at face " << faceI
                    << " l:" << l[faceI] << " u:" << u[faceI]
                    << " previous u:" << u[faceI-1]
                    //<< abort(FatalError);
                    << endl;
            }
        }
    }
}


Foam::label Foam::lduPrimitiveMesh::size(const PtrList<lduMesh>& meshes)
{
    label size = 0;

    forAll(meshes, i)
    {
        size += meshes[i].lduAddr().size();
    }
    return size;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label nCells,
    const labelUList& l,
    const labelUList& u,
    const labelListList& pa,
    const lduInterfacePtrsList interfaces,
    const lduSchedule& ps,
    const label comm
)
:
    lduAddressing(nCells),
    lowerAddr_(l),
    upperAddr_(u),
    patchAddr_(pa),
    interfaces_(interfaces),
    patchSchedule_(ps),
    comm_(comm)
{
    Pout<< "lduPrimitiveMesh :"
        << " nCells:" << nCells
        << " l:" << lowerAddr_.size()
        << " u:" << upperAddr_.size()
        << " pa:" << patchAddr_.size()
        << " interfaces:" << interfaces_.size()
        << " comm:" << comm_
        << endl;
    forAll(interfaces_, i)
    {
        if (interfaces_.set(i))
        {
            if (isA<processorLduInterface>(interfaces_[i]))
            {
                const processorLduInterface& pi = refCast
                <
                    const processorLduInterface
                >(interfaces_[i]);

                Pout<< "    patch:" << i
                    << " size:" << patchAddr_[i].size()
                    << " myProcNo:" << pi.myProcNo()
                    << " neighbProcNo:" << pi.neighbProcNo()
                    << " comm:" << pi.comm()
                    << endl;
            }
        }
    }
}


Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label nCells,
    labelList& l,
    labelList& u,
    labelListList& pa,
    lduInterfacePtrsList interfaces,
    const lduSchedule& ps,
    const label comm,
    bool reUse
)
:
    lduAddressing(nCells),
    lowerAddr_(l, reUse),
    upperAddr_(u, reUse),
    patchAddr_(pa, reUse),
    interfaces_(interfaces, reUse),
    patchSchedule_(ps),
    comm_(comm)
{
    Pout<< "lduPrimitiveMesh :"
        << " nCells:" << nCells
        << " l:" << lowerAddr_.size()
        << " u:" << upperAddr_.size()
        << " pa:" << patchAddr_.size()
        << " interfaces:" << interfaces_.size()
        << " comm:" << comm_
        << endl;
    forAll(interfaces_, i)
    {
        if (interfaces_.set(i))
        {
            if (isA<processorLduInterface>(interfaces_[i]))
            {
                const processorLduInterface& pi = refCast
                <
                    const processorLduInterface
                >(interfaces_[i]);

                Pout<< "    patch:" << i
                    << " size:" << patchAddr_[i].size()
                    << " myProcNo:" << pi.myProcNo()
                    << " neighbProcNo:" << pi.neighbProcNo()
                    << " comm:" << pi.comm()
                    << endl;
            }
        }
    }
}


Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label comm,
    const labelList& procIDs,
    const lduMesh& myMesh,
    const PtrList<lduMesh>& otherMeshes,

    labelList& cellOffsets,
    labelListList& faceMap,
    labelListList& boundaryMap,
    labelListListList& boundaryFaceMap
)
:
    lduAddressing(myMesh.lduAddr().size() + size(otherMeshes)),
    lowerAddr_(0),
    upperAddr_(0),
    patchAddr_(0),
    interfaces_(0),
    patchSchedule_(0),
    comm_(comm)
{
    // Sanity check.
    for (label i = 1; i < procIDs.size(); i++)
    {
        if (procIDs[i] <= procIDs[i-1])
        {
            FatalErrorIn
            (
                "lduPrimitiveMesh::lduPrimitiveMesh(..)"
            )   << "Processor " << procIDs[i] << " at index " << i
                << " should be higher numbered than its predecessor "
                << procIDs[i-1]
                << exit(FatalError);
        }
    }

    const label currentComm = myMesh.comm();

    forAll(otherMeshes, i)
    {
        if (otherMeshes[i].comm() != currentComm)
        {
            FatalErrorIn
            (
                "lduPrimitiveMesh::lduPrimitiveMesh(..)"
            )   << "Communicator " << otherMeshes[i].comm()
                << " at index " << i
                << " differs from that of predecessor "
                << currentComm
                << exit(FatalError);
        }
    }

    const label nMeshes = otherMeshes.size()+1;

    // Cells get added in order.
    cellOffsets.setSize(nMeshes+1);
    cellOffsets[0] = 0;
    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduMesh& procMesh = mesh(myMesh, otherMeshes, procMeshI);

        cellOffsets[procMeshI+1] =
            cellOffsets[procMeshI]
          + procMesh.lduAddr().size();
    }

    // Faces initially get added in order, sorted later
    labelList internalFaceOffsets(nMeshes+1);
    internalFaceOffsets[0] = 0;
    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduMesh& procMesh = mesh(myMesh, otherMeshes, procMeshI);

        internalFaceOffsets[procMeshI+1] =
            internalFaceOffsets[procMeshI]
          + procMesh.lduAddr().lowerAddr().size();
    }

    // Count how faces get added. Interfaces inbetween get merged.

    // Merged interfaces: map from two processors back to
    // - procMeshes
    // - interface in procMesh
    EdgeMap<labelPairList> mergedMap
    (
        2*myMesh.interfaces().size()
    );

    // Unmerged interfaces: map from two processors back to
    // - procMeshes
    // - interface in procMesh
    EdgeMap<labelPairList> unmergedMap(mergedMap.size());

    boundaryMap.setSize(nMeshes);
    boundaryFaceMap.setSize(nMeshes);


    label nBoundaryFaces = 0;
    label nInterfaces = 0;
    labelList nCoupledFaces(nMeshes, 0);

    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduInterfacePtrsList interfaces =
            mesh(myMesh, otherMeshes, procMeshI).interfaces();

        // Inialise all boundaries as merged
        boundaryMap[procMeshI].setSize(interfaces.size(), -1);
        boundaryFaceMap[procMeshI].setSize(interfaces.size());

        // Get sorted order of processors
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                if (isA<processorLduInterface>(interfaces[intI]))
                {
                    const processorLduInterface& pldui =
                        refCast<const processorLduInterface>
                        (
                            interfaces[intI]
                        );
                    if (pldui.myProcNo() != procIDs[procMeshI])
                    {
                        FatalErrorIn
                        (
                            "lduPrimitiveMesh::lduPrimitiveMesh(..)"
                        )   << "proc:" << procIDs[procMeshI]
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
                            << " size:" << interfaces[intI].faceCells().size()
                            << endl;

                        nBoundaryFaces += interfaces[intI].faceCells().size();

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
                    else
                    {
                        // Merged interface
                        Pout<< "merged interface: myProcNo:" << pldui.myProcNo()
                            << " nbr:" << pldui.neighbProcNo()
                            << " size:" << interfaces[intI].faceCells().size()
                            << endl;
                        if (pldui.myProcNo() < pldui.neighbProcNo())
                        {
                            nCoupledFaces[procMeshI] +=
                                interfaces[intI].faceCells().size();
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
                    FatalErrorIn("lduPrimitiveMesh::lduPrimitiveMesh(..)")
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
                    mesh(myMesh, otherMeshes, procMeshI).interfaces();
                const processorLduInterface& pldui =
                    refCast<const processorLduInterface>
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
                    mesh(myMesh, otherMeshes, procMeshI).interfaces();
                const processorLduInterface& pldui =
                    refCast<const processorLduInterface>
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
    labelList faceOffsets(nMeshes+1);
    faceOffsets[0] = 0;
    faceMap.setSize(nMeshes);
    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduMesh& procMesh = mesh(myMesh, otherMeshes, procMeshI);
        label nInternal = procMesh.lduAddr().lowerAddr().size();

        faceOffsets[procMeshI+1] =
            faceOffsets[procMeshI]
          + nInternal
          + nCoupledFaces[procMeshI];

        labelList& map = faceMap[procMeshI];
        map.setSize(nInternal);
        forAll(map, i)
        {
            map[i] = faceOffsets[procMeshI] + i;
        }
    }


    // Combine upper and lower
    lowerAddr_.setSize(faceOffsets.last(), -1);
    upperAddr_.setSize(lowerAddr_.size(), -1);


    // Old internal faces and resolved coupled interfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduMesh& procMesh = mesh(myMesh, otherMeshes, procMeshI);

        const labelUList& l = procMesh.lduAddr().lowerAddr();
        const labelUList& u = procMesh.lduAddr().upperAddr();

        // Add internal faces
        label allFaceI = faceOffsets[procMeshI];

        forAll(l, faceI)
        {
            lowerAddr_[allFaceI] = cellOffsets[procMeshI]+l[faceI];
            upperAddr_[allFaceI] = cellOffsets[procMeshI]+u[faceI];
            allFaceI++;
        }

        // Add merged interfaces
        const lduInterfacePtrsList interfaces = procMesh.interfaces();

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                if (isA<processorLduInterface>(interfaces[intI]))
                {
                    const processorLduInterface& pldui =
                        refCast<const processorLduInterface>
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
                                FatalErrorIn
                                (
                                    "lduPrimitiveMesh::lduPrimitiveMesh(..)"
                                )   << "elems:" << elems << abort(FatalError);
                            }


                            const lduInterfacePtrsList nbrInterfaces =
                                mesh
                                (
                                    myMesh,
                                    otherMeshes,
                                    nbrProcMeshI
                                ).interfaces();


                            const labelUList& faceCells =
                                interfaces[intI].faceCells();
                            const labelUList& nbrFaceCells =
                                nbrInterfaces[nbrIntI].faceCells();

                            labelList& bfMap =
                                boundaryFaceMap[procMeshI][intI];
                            labelList& nbrBfMap =
                                boundaryFaceMap[nbrProcMeshI][nbrIntI];

                            bfMap.setSize(faceCells.size(), -1);
                            nbrBfMap.setSize(faceCells.size(), -1);

                            forAll(faceCells, pfI)
                            {
                                lowerAddr_[allFaceI] =
                                    cellOffsets[procMeshI]+faceCells[pfI];
                                bfMap[pfI] = allFaceI;
                                upperAddr_[allFaceI] =
                                    cellOffsets[nbrProcMeshI]+nbrFaceCells[pfI];
                                nbrBfMap[pfI] = (-allFaceI-1);
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

    interfaces_.setSize(nInterfaces);
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
                mesh
                (
                    myMesh,
                    otherMeshes,
                    procMeshI
                ).interfaces();

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
                mesh
                (
                    myMesh,
                    otherMeshes,
                    procMeshI
                ).interfaces();

            boundaryMap[procMeshI][interfaceI] = allInterfaceI;
            labelList& bfMap = boundaryFaceMap[procMeshI][interfaceI];

            const labelUList& l = interfaces[interfaceI].faceCells();
            forAll(l, faceI)
            {
                allFaceCells[n] = cellOffsets[procMeshI]+l[faceI];
                allFaceRestrictAddressing[n] = n;
                bfMap[faceI] = faceI;
                n++;
            }
        }


        // Find out local and remote processor in new communicator

        label myProcNo = -1;
        label neighbProcNo = -1;

        if (findIndex(procIDs, iter.key()[0]) != -1)
        {
            myProcNo = UPstream::procNo
            (
                comm_,
                currentComm,
                iter.key()[0]
            );
            neighbProcNo = UPstream::procNo
            (
                comm_,
                currentComm,
                iter.key()[1]
            );
        }
        else
        {
            myProcNo = UPstream::procNo
            (
                comm_,
                currentComm,
                iter.key()[1]
            );
            neighbProcNo = UPstream::procNo
            (
                comm_,
                currentComm,
                iter.key()[0]
            );
        }


        interfaces_.set
        (
            allInterfaceI,
            new processorGAMGInterface
            (
                allInterfaceI,
                interfaces_,
                allFaceCells,
                allFaceRestrictAddressing,
                comm_,
                myProcNo,
                neighbProcNo,
                tensorField(),          // forwardT
                Pstream::msgType()      // tag
            )
        );
        allInterfaceI++;
    }


    // Extract faceCells from interfaces_
    labelListList patchAddr_(interfaces_.size());
    forAll(interfaces_, coarseIntI)
    {
        if (interfaces_.set(coarseIntI))
        {
            patchAddr_[coarseIntI] =
                interfaces_[coarseIntI].faceCells();
        }
    }

    patchSchedule_ = nonBlockingSchedule<processorGAMGInterface>(interfaces_);

    Pout<< "lowerAddr_:" << lowerAddr_ << endl;
    Pout<< "upperAddr_:" << upperAddr_ << endl;
    checkUpperTriangular(cellOffsets.last(), lowerAddr_, upperAddr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduMesh& Foam::lduPrimitiveMesh::mesh
(
    const lduMesh& myMesh,
    const PtrList<lduMesh>& otherMeshes,
    const label meshI
)
{
    return (meshI == 0 ? myMesh : otherMeshes[meshI-1]);
}


void Foam::lduPrimitiveMesh::gather
(
    const lduMesh& mesh,
    const labelList& procIDs,
    PtrList<lduMesh>& otherMeshes
)
{
    // Force calculation of schedule (since does parallel comms)
    (void)mesh.lduAddr().patchSchedule();


    const label meshComm = mesh.comm();

    if (Pstream::myProcNo(meshComm) == procIDs[0])
    {
        otherMeshes.setSize(procIDs.size()-1);

        // Slave meshes
        for (label i = 1; i < procIDs.size(); i++)
        {
            //Pout<< "on master :"
            //    << " receiving from slave " << procIDs[i] << endl;

            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,          // bufSize
                Pstream::msgType(),
                meshComm
            );

            label nCells = readLabel(fromSlave);
            labelList lowerAddr(fromSlave);
            labelList upperAddr(fromSlave);
            boolList validInterface(fromSlave);

            // Construct GAMGInterfaces
            lduInterfacePtrsList newInterfaces(validInterface.size());
            labelListList patchAddr(validInterface.size());
            forAll(validInterface, intI)
            {
                if (validInterface[intI])
                {
                    word coupleType(fromSlave);

                    newInterfaces.set
                    (
                        intI,
                        GAMGInterface::New
                        (
                            coupleType,
                            intI,
                            newInterfaces,
                            fromSlave
                        ).ptr()
                    );
                    patchAddr[intI] = newInterfaces[intI].faceCells();
                }
            }

            otherMeshes.set
            (
                i-1,
                new lduPrimitiveMesh
                (
                    nCells,
                    lowerAddr,
                    upperAddr,
                    patchAddr,
                    newInterfaces,
                    nonBlockingSchedule<processorGAMGInterface>
                    (
                        newInterfaces
                    ),
                    meshComm,
                    true
                )
            );
        }
    }
    else if (findIndex(procIDs, Pstream::myProcNo(meshComm)) != -1)
    {
        // Send to master

        const lduAddressing& addressing = mesh.lduAddr();
        lduInterfacePtrsList interfaces(mesh.interfaces());
        boolList validInterface(interfaces.size());
        forAll(interfaces, intI)
        {
            validInterface[intI] = interfaces.set(intI);
        }

        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            Pstream::msgType(),
            meshComm
        );

        //Pout<< "sent nCells:" << addressing.size()
        //    << "sent lowerAddr:" << addressing.lowerAddr().size()
        //    << "sent upperAddr:" << addressing.upperAddr()
        //    << "sent validInterface:" << validInterface
        //    << endl;

        toMaster
            << addressing.size()
            << addressing.lowerAddr()
            << addressing.upperAddr()
            << validInterface;

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const GAMGInterface& interface = refCast<const GAMGInterface>
                (
                    interfaces[intI]
                );

                toMaster << interface.type();
                interface.write(toMaster);
            }
        }
    }
}


// ************************************************************************* //
