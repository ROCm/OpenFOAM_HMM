/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "zoneDistribute.H"
#include "dummyTransform.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"

#include "globalPoints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneDistribute, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::zoneDistribute::coupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<processorPolyPatch>(pp))
        {
            nCoupled += pp.size();
        }
    }
    labelList coupledFaces(nCoupled);
    nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<processorPolyPatch>(pp))
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                coupledFaces[nCoupled++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            coupledFaces
        ),
        mesh_.points()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneDistribute::zoneDistribute(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::TopologicalMeshObject, zoneDistribute>(mesh),
    coupledBoundaryPoints_(coupledFacesPatch()().meshPoints()),
    send_(UPstream::nProcs()),
    stencil_(zoneCPCStencil::New(mesh)),
    globalNumbering_(stencil_.globalNumbering()),
    sendTo_(),   // Initial zero-sized
    recvFrom_()  // Initial zero-sized
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

Foam::zoneDistribute& Foam::zoneDistribute::New(const fvMesh& mesh)
{
    auto* ptr = mesh.thisDb().getObjectPtr<zoneDistribute>("zoneDistribute");

    if (!ptr)
    {
        ptr = new zoneDistribute(mesh);
        regIOobject::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneDistribute::updateStencil(const boolList& zone)
{
    zoneCPCStencil::New(mesh_).updateStencil(zone);
}


void Foam::zoneDistribute::setUpCommforZone
(
    const boolList& zone,
    bool updateStencil
)
{
    zoneCPCStencil& stencil = zoneCPCStencil::New(mesh_);

    if (updateStencil)
    {
        stencil.updateStencil(zone);
    }

    if (UPstream::parRun())
    {
        if (sendTo_.empty())
        {
            // First time
            sendTo_.resize(UPstream::nProcs());
            recvFrom_.resize(UPstream::nProcs());
            sendTo_ = false;
            recvFrom_ = false;
        }

        const labelHashSet& comms = stencil.needsComm();

        List<labelHashSet> needed(UPstream::nProcs());

        for (const label celli : comms)
        {
            if (zone[celli])
            {
                for (const label gblIdx : stencil_[celli])
                {
                    if (!globalNumbering_.isLocal(gblIdx))
                    {
                        const label procID =
                            globalNumbering_.whichProcID(gblIdx);
                        needed[procID].insert(gblIdx);
                    }
                }
            }
        }


        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);
        labelList recvSizes;

        // Stream data into buffer
        for (const int proci : UPstream::allProcs())
        {
            if (proci != UPstream::myProcNo() && !needed[proci].empty())
            {
                // Put data into send buffer
                UOPstream toProc(proci, pBufs);

                toProc << needed[proci].sortedToc();
            }
        }


        // Need update, or use existing partial communication info?
        bool fullUpdate = false;
        for (const int proci : UPstream::allProcs())
        {
            if
            (
                proci != UPstream::myProcNo()
             && (!sendTo_[proci] && !needed[proci].empty())
            )
            {
                // Changed state sendTo_ from false -> true
                sendTo_[proci] = true;
                fullUpdate = true;
            }
        }


        if (returnReduce(fullUpdate, orOp<bool>()))
        {
            pBufs.finishedSends(recvSizes);

            // Update which ones receive
            for (const int proci : UPstream::allProcs())
            {
                recvFrom_[proci] = (recvSizes[proci] > 0);
            }
        }
        else
        {
            // No change in senders...
            // - can communicate with a subset of processors
            DynamicList<label> sendProcs;
            DynamicList<label> recvProcs;

            for (const int proci : UPstream::allProcs())
            {
                if (sendTo_[proci])
                {
                    sendProcs.append(proci);
                }
                if (recvFrom_[proci])
                {
                    recvProcs.append(proci);
                }
            }

            // Wait until everything is written
            pBufs.finishedSends(sendProcs, recvProcs, recvSizes);
        }

        for (const int proci : UPstream::allProcs())
        {
            send_[proci].clear();

            if
            (
                proci != UPstream::myProcNo()
             && recvFrom_[proci]   // Or: (recvSizes[proci] > 0)
            )
            {
                UIPstream fromProc(proci, pBufs);
                fromProc >> send_[proci];
            }
        }
    }
}


// ************************************************************************* //
