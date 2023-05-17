/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

Note
    The algorithm NBX (Nonblocking consensus exchange) is described by

    "Scalable Communication Protocols for Dynamic Sparse Data Exchange",
    Hoeffler, Siebert, Lumsdaine
    May 2010 ACM SIGPLAN Notices 45(5):159-168
    https://doi.org/10.1145/1837853.1693476

    http://unixer.de/publications/img/hoefler-dsde-protocols.pdf

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "contiguous.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

//- Exchange \em contiguous data using non-blocking consensus exchange
//- with optional tracking of the receive sizes.
//
//  No internal guards or resizing - data containers are all properly
//  sized before calling.
//
//  \param[in]  sendBufs  The send buffers list (size: numProcs)
//  \param[out] recvBufs  The recv buffers list (size: numProcs)
//  \param[out] recvSizes The recv sizes (size: 0 or numProcs).
//     This parameter can be an empty list, in which case the receive sizes
//     are not returned.
//  \param tag   The message tag
//  \param comm  The communicator
//  \param wait  Wait for non-blocking receives to complete
//  \param recvCommType If blocking or (default) non-blocking

template<class Container, class Type>
void exchangeConsensus
(
    const UList<Container>& sendBufs,
    UList<Container>& recvBufs,
    labelUList& recvSizes,
    const int tag,
    const label comm,
    const bool wait = true,
    const UPstream::commsTypes recvCommType = UPstream::commsTypes::nonBlocking
)
{
    static_assert(is_contiguous<Type>::value, "Contiguous data only!");

    const bool initialBarrier = (UPstream::tuning_NBX_ > 0);

    const label startOfRequests = UPstream::nRequests();
    const label myProci = UPstream::myProcNo(comm);
    const label numProc = UPstream::nProcs(comm);

    // Initial: clear all receive information
    for (auto& buf : recvBufs)
    {
        buf.clear();
    }
    recvSizes = Zero;

    if (!UPstream::is_rank(comm))
    {
        return;  // Process not in communicator
    }

    // #ifdef FULLDEBUG
    if (sendBufs.size() > numProc)
    {
        FatalErrorInFunction
            << "Send buffers:" << sendBufs.size() << " > numProcs:" << numProc
            << Foam::abort(FatalError);
    }
    if (recvBufs.size() < numProc)
    {
        FatalErrorInFunction
            << "Recv buffers:" << recvBufs.size() << " < numProcs:" << numProc
            << Foam::abort(FatalError);
    }
    // #endif

    if (!UPstream::is_parallel(comm))
    {
        // Do myself
        recvBufs[myProci] = sendBufs[myProci];
        if (myProci < recvSizes.size())
        {
            recvSizes[myProci] = recvBufs.size();
        }
        return;
    }

    // An initial barrier may help to avoid synchronisation problems
    // caused elsewhere
    if (initialBarrier)
    {
        UPstream::barrier(comm);
    }

    // Algorithm NBX: Nonblocking consensus with List containers

    DynamicList<UPstream::Request> sendRequests(sendBufs.size());

    // Start nonblocking synchronous send to processor dest
    for (label proci = 0; proci < numProc; ++proci)
    {
        const auto& sendData = sendBufs[proci];

        if (sendData.empty())
        {
            // Do not send/recv empty data
        }
        else if (proci == myProci)
        {
            // Do myself
            recvBufs[proci] = sendData;
            if (proci < recvSizes.size())
            {
                recvSizes[proci] = sendData.size();
            }
        }
        else
        {
            // Has data to send.
            // The MPI send requests are tracked on a local list
            UOPstream::write
            (
                sendRequests.emplace_back(),
                proci,
                sendData.cdata_bytes(),
                sendData.size_bytes(),
                tag,
                comm,
                UPstream::sendModes::sync
            );
        }
    }


    // Probe and receive

    UPstream::Request barrierRequest;

    for (bool barrier_active = false, done = false; !done; /*nil*/)
    {
        std::pair<int, int> probed =
            UPstream::probeMessage
            (
                UPstream::commsTypes::nonBlocking,
                -1,  // ANY_SOURCE
                tag,
                comm
            );

        if (probed.second > 0)
        {
            // Message found and had size.
            // - receive into dest buffer location

            const label proci = probed.first;
            const label count = (probed.second / sizeof(Type));

            auto& recvData = recvBufs[proci];
            recvData.resize_nocopy(count);

            if (proci < recvSizes.size())
            {
                recvSizes[proci] = count;
            }

            // Any non-blocking MPI recv requests are tracked on internal stack
            UIPstream::read
            (
                recvCommType,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes(),
                tag,
                comm
            );
        }

        if (barrier_active)
        {
            // Test barrier for completion
            // - all received, or nothing to receive
            if (UPstream::finishedRequest(barrierRequest))
            {
                done = true;
            }
        }
        else
        {
            // Check if all sends have arrived
            if (UPstream::finishedRequests(sendRequests))
            {
                UPstream::barrier(comm, &barrierRequest);
                barrier_active = true;
            }
        }
    }

    // Wait for non-blocking receives to finish
    if (wait && recvCommType == UPstream::commsTypes::nonBlocking)
    {
        UPstream::waitRequests(startOfRequests);
    }
}


//- Exchange \em contiguous data using non-blocking consensus exchange.
//
//  No internal guards - the sending Map corresponds to a segment of
//  0-numProcs.
//
//  \param[in]  sendBufs  The send buffers map (addr: 0-numProcs)
//  \param[out] recvBufs  The recv buffers map
//  \param tag   The message tag
//  \param comm  The communicator
//  \param wait  Wait for non-blocking receives to complete
//  \param recvCommType If blocking or (default) non-blocking

template<class Container, class Type>
void exchangeConsensus
(
    const Map<Container>& sendBufs,
    Map<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait = true,
    const UPstream::commsTypes recvCommType = UPstream::commsTypes::nonBlocking
)
{
    static_assert(is_contiguous<Type>::value, "Contiguous data only!");

    const label startOfRequests = UPstream::nRequests();
    const label myProci = UPstream::myProcNo(comm);

    // Initial: clear out receive 'slots'
    // Preferrable to clear out the map entries instead of the map itself
    // since this can potentially preserve allocated space
    // (eg DynamicList entries) between calls

    forAllIters(recvBufs, iter)
    {
        iter.val().clear();
    }

    if (!UPstream::is_parallel(comm))
    {
        // Do myself
        const auto iter = sendBufs.find(myProci);
        if (iter.good())
        {
            const auto& sendData = iter.val();

            if (!sendData.empty())
            {
                // Do myself: insert_or_assign
                recvBufs(iter.key()) = sendData;
            }
        }
        return;
    }


    // Algorithm NBX: Nonblocking consensus with Map (HashTable) containers

    DynamicList<UPstream::Request> sendRequests(sendBufs.size());

    // Start nonblocking synchronous send to process dest
    forAllConstIters(sendBufs, iter)
    {
        const label proci = iter.key();
        const auto& sendData = iter.val();

        #ifdef FULLDEBUG
        if (proci >= UPstream::nProcs(comm))
        {
            FatalErrorInFunction
                << "Send buffer:" << proci << " >= numProcs:"
                << UPstream::nProcs(comm)
                << Foam::abort(FatalError);
        }
        #endif

        if (sendData.empty())
        {
            // Do not send/recv empty data
        }
        else if (proci == myProci)
        {
            // Do myself: insert_or_assign
            recvBufs(proci) = sendData;
        }
        else
        {
            // Has data to send.
            // The MPI send requests are tracked on a local list
            UOPstream::write
            (
                sendRequests.emplace_back(),
                proci,
                sendData.cdata_bytes(),
                sendData.size_bytes(),
                tag,
                comm,
                UPstream::sendModes::sync
            );
        }
    }


    // Probe and receive

    UPstream::Request barrierRequest;

    for (bool barrier_active = false, done = false; !done; /*nil*/)
    {
        std::pair<int, int> probed =
            UPstream::probeMessage
            (
                UPstream::commsTypes::nonBlocking,
                -1,  // ANY_SOURCE
                tag,
                comm
            );

        if (probed.second > 0)
        {
            // Message found and had size.
            // - receive into dest buffer location

            const label proci = probed.first;
            const label count = (probed.second / sizeof(Type));

            auto& recvData = recvBufs(proci);
            recvData.resize_nocopy(count);

            // Any non-blocking MPI recv requests are tracked on internal stack
            UIPstream::read
            (
                recvCommType,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes(),
                tag,
                comm
            );
        }

        if (barrier_active)
        {
            // Test barrier for completion
            if (UPstream::finishedRequest(barrierRequest))
            {
                done = true;
            }
        }
        else
        {
            // Check if all sends have arrived
            if (UPstream::finishedRequests(sendRequests))
            {
                UPstream::barrier(comm, &barrierRequest);
                barrier_active = true;
            }
        }
    }

    // Wait for non-blocking receives to finish
    if (wait && recvCommType == UPstream::commsTypes::nonBlocking)
    {
        UPstream::waitRequests(startOfRequests);
    }
}

} // namespace PstreamDetail
} // namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, class Type>
void Foam::Pstream::exchangeConsensus
(
    const UList<Container>& sendBufs,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    static_assert(is_contiguous<Type>::value, "Contiguous data only!");

    if (sendBufs.size() != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
            << "Send buffers size:" << sendBufs.size()
            << " != numProcs:" << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    // Resize receive buffers. Individual clearing is done internally
    recvBufs.resize_nocopy(sendBufs.size());
    labelList dummyRecvSizes;

    PstreamDetail::exchangeConsensus
    (
        sendBufs,
        recvBufs,
        dummyRecvSizes,
        tag,
        comm,
        wait
    );
}


template<class Container, class Type>
void Foam::Pstream::exchangeConsensus
(
    const Map<Container>& sendBufs,
    Map<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    static_assert(is_contiguous<Type>::value, "Contiguous data only!");

    PstreamDetail::exchangeConsensus
    (
        sendBufs,
        recvBufs,
        tag,
        comm,
        wait
    );
}


// ************************************************************************* //
