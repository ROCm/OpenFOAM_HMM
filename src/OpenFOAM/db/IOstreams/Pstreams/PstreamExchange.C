/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "Pstream.H"
#include "contiguous.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

//- Setup sends and receives, each specified as [rank, span] tuple
//  The serial list of tuples can be populated from other lists, from maps
//  of data or subsets of lists/maps etc.
template<class Type>
void exchangeBuf
(
    const UList<std::pair<int, stdFoam::span<const Type>>>& sends,
    const UList<std::pair<int, stdFoam::span<Type>>>& recvs,

    const int tag,
    const label comm,
    const bool wait
)
{
    const label startOfRequests = UPstream::nRequests();
    const int myProci = UPstream::myProcNo(comm);

    // Set up receives
    // ~~~~~~~~~~~~~~~

    for (auto& slot : recvs)
    {
        // [rank, span]
        const auto proci = slot.first;
        auto& payload = slot.second;

        if (proci != myProci && !payload.empty())
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                payload.data_bytes(),
                payload.size_bytes(),
                tag,
                comm
            );
        }
    }

    // Set up sends
    // ~~~~~~~~~~~~

    for (const auto& slot : sends)
    {
        // [rank, span]
        const auto proci = slot.first;
        const auto& payload = slot.second;

        if (proci != myProci && !payload.empty())
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    payload.cdata_bytes(),
                    payload.size_bytes(),
                    tag,
                    comm
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message to:"
                    << proci << " nBytes:"
                    << label(payload.size_bytes())
                    << Foam::abort(FatalError);
            }
        }
    }

    // Wait for all to finish
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (wait)
    {
        UPstream::waitRequests(startOfRequests);
    }
}


//- Chunked exchange of \em contiguous data.
//- Setup sends and receives, each specified as [rank, span] tuple.
//  The serial list of tuples can be populated from other lists, from
//  maps of data or subsets of lists/maps etc.
template<class Type>
void exchangeChunkedBuf
(
    const UList<std::pair<int, stdFoam::span<const Type>>>& sends,
    const UList<std::pair<int, stdFoam::span<Type>>>& recvs,

    const int tag,
    const label comm,
    const bool wait
)
{
    typedef std::pair<int, stdFoam::span<const Type>> sendTuple;
    typedef std::pair<int, stdFoam::span<Type>> recvTuple;

    // Caller already checked for parRun and maxChunkSize > 0

    {
        // Determine the number of chunks to send. Note that we
        // only have to look at the sending data since we are
        // guaranteed that some processor's sending size is some other
        // processor's receive size. Also we can ignore any local comms.
        //
        // We need to send chunks so the number of iterations:
        //  maxChunkSize                        iterations
        //  ------------                        ----------
        //  0                                   0
        //  1..maxChunkSize                     1
        //  maxChunkSize+1..2*maxChunkSize      2
        //  ...

        const label maxChunkSize =
        (
            max
            (
                static_cast<label>(1),
                static_cast<label>(UPstream::maxCommsSize/sizeof(Type))
            )
        );

        const int myProci = UPstream::myProcNo(comm);

        label nChunks(0);
        {
            // Get max send count (elements)
            auto maxCount = static_cast<stdFoam::span<char>::size_type>(0);

            for (const auto& slot : sends)
            {
                // [rank, span]
                const auto proci = slot.first;
                const auto count = slot.second.size();

                if (proci != myProci && count > maxCount)
                {
                    // Note: using max() can be ambiguous
                    maxCount = count;
                }
            }

            // Convert from send count (elements) to number of chunks.
            // Can normally calculate with (count-1), but add some safety
            if (maxCount)
            {
                nChunks = 1 + label(maxCount/maxChunkSize);
            }

            // MPI reduce (message tag is irrelevant)
            reduce(nChunks, maxOp<label>(), UPstream::msgType(), comm);
        }


        // Dispatch the exchanges chunk-wise
        List<sendTuple> sendChunks(sends);
        List<recvTuple> recvChunks(recvs);

        // Dispatch
        for (label iter = 0; iter < nChunks; ++iter)
        {
            // The begin/end for the data window
            const auto beg = static_cast<std::size_t>(iter*maxChunkSize);
            const auto end = static_cast<std::size_t>((iter+1)*maxChunkSize);

            forAll(sendChunks, sloti)
            {
                const auto& baseline = sends[sloti].second;
                auto& payload = sendChunks[sloti].second;

                // Window the data
                if (beg < baseline.size())
                {
                    payload =
                    (
                        (end < baseline.size())
                      ? baseline.subspan(beg, end - beg)
                      : baseline.subspan(beg)
                    );
                }
                else
                {
                    payload = baseline.first(0);  // zero-sized
                }
            }

            forAll(recvChunks, sloti)
            {
                const auto& baseline = recvs[sloti].second;
                auto& payload = recvChunks[sloti].second;

                // Window the data
                if (beg < baseline.size())
                {
                    payload =
                    (
                        (end < baseline.size())
                      ? baseline.subspan(beg, end - beg)
                      : baseline.subspan(beg)
                    );
                }
                else
                {
                    payload = baseline.first(0);  // zero-sized
                }
            }


            // Exchange data chunks
            PstreamDetail::exchangeBuf<Type>
            (
                sendChunks,
                recvChunks,
                tag,
                comm,
                wait
            );

            // Debugging output - can report on master only...
            #if 0 // ifdef Foam_PstreamExchange_debug_chunks
            do
            {
                labelList sendStarts(sends.size());
                labelList sendCounts(sends.size());

                forAll(sendChunks, sloti)
                {
                    const auto& baseline = sends[sloti].second;
                    const auto& payload = sendChunks[sloti].second;

                    sendStarts[sloti] = (payload.data() - baseline.data());
                    sendCounts[sloti] = (payload.size());
                }

                Info<< "iter " << iter
                    << ": beg=" << flatOutput(sendStarts)
                    << " len=" << flatOutput(sendCounts) << endl;
            } while (false);
            #endif
        }
    }
}


//- Exchange \em contiguous data using point-to-point communication.
//- Sends sendBufs, receives into recvBufs.
//  Data provided and received as container all of which have been
//  properly sized before calling
//
// No internal guards or resizing.
template<class Container, class Type>
void exchangeContainer
(
    const UList<Container>& sendBufs,
    UList<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait               //!< Wait for requests to complete
)
{
    const label startOfRequests = UPstream::nRequests();
    const label myProci = UPstream::myProcNo(comm);

    // Set up receives
    // ~~~~~~~~~~~~~~~

    forAll(recvBufs, proci)
    {
        auto& recvData = recvBufs[proci];

        if (proci != myProci && !recvData.empty())
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes(),
                tag,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    forAll(sendBufs, proci)
    {
        const auto& sendData = sendBufs[proci];

        if (proci != myProci && !sendData.empty())
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendData.cdata_bytes(),
                    sendData.size_bytes(),
                    tag,
                    comm
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message. "
                    << "to:" << proci << " nBytes:"
                    << label(sendData.size_bytes())
                    << Foam::abort(FatalError);
            }
        }
    }

    // Wait for all to finish
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (wait)
    {
        UPstream::waitRequests(startOfRequests);
    }
}


//- Exchange \em contiguous data using point-to-point communication.
//- Sends sendBufs, receives into recvBufs.
//  Data provided and received as container all of which have been
//  properly sized before calling
//
// No internal guards or resizing.
template<class Container, class Type>
void exchangeContainer
(
    const Map<Container>& sendBufs,
    Map<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait               //!< Wait for requests to complete
)
{
    const label startOfRequests = UPstream::nRequests();
    const label myProci = UPstream::myProcNo(comm);

    // Set up receives
    // ~~~~~~~~~~~~~~~

    forAllIters(recvBufs, iter)
    {
        const label proci = iter.key();
        auto& recvData = iter.val();

        if (proci != myProci && !recvData.empty())
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes(),
                tag,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    forAllConstIters(sendBufs, iter)
    {
        const label proci = iter.key();
        const auto& sendData = iter.val();

        if (proci != myProci && !sendData.empty())
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendData.cdata_bytes(),
                    sendData.size_bytes(),
                    tag,
                    comm
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message to:"
                    << proci << " nBytes:"
                    << label(sendData.size_bytes())
                    << Foam::abort(FatalError);
            }
        }
    }

    // Wait for all to finish
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (wait)
    {
        UPstream::waitRequests(startOfRequests);
    }
}

} // namespace PstreamDetail
} // namespace Foam

#include "PstreamExchangeConsensus.C"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, class Type>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    const labelUList& recvSizes,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    static_assert(is_contiguous<Type>::value, "Contiguous data only!");

    if (!UPstream::is_rank(comm))
    {
        return;  // Process not in communicator
    }

    const label myProci = UPstream::myProcNo(comm);
    const label numProcs = UPstream::nProcs(comm);

    if (sendBufs.size() != numProcs)
    {
        FatalErrorInFunction
            << "Size of list " << sendBufs.size()
            << " does not equal the number of processors " << numProcs
            << Foam::abort(FatalError);
    }

    recvBufs.resize_nocopy(numProcs);

    if (UPstream::is_parallel(comm))
    {
        // Presize all receive buffers
        forAll(recvSizes, proci)
        {
            const label count = recvSizes[proci];

            if (proci != myProci && count > 0)
            {
                recvBufs[proci].resize_nocopy(count);
            }
            else
            {
                recvBufs[proci].clear();
            }
        }

        typedef std::pair<int, stdFoam::span<const Type>> sendTuple;
        typedef std::pair<int, stdFoam::span<Type>> recvTuple;

        if (UPstream::maxCommsSize <= 0)
        {
            // Do the exchanging in one go
            PstreamDetail::exchangeContainer<Container, Type>
            (
                sendBufs,
                recvBufs,
                tag,
                comm,
                wait
            );
        }
        else
        {
            // Dispatch using chunk-wise exchanges
            // Populate send sequence
            DynamicList<sendTuple> sends(sendBufs.size());
            forAll(sendBufs, proci)
            {
                const auto& sendData = sendBufs[proci];

                if (proci != myProci && !sendData.empty())
                {
                    sends.push_back
                    (
                        sendTuple
                        (
                            proci,
                            { sendData.cdata(), std::size_t(sendData.size()) }
                        )
                    );
                }
            }

            // Populate recv sequence
            DynamicList<recvTuple> recvs(recvBufs.size());
            forAll(recvBufs, proci)
            {
                auto& recvData = recvBufs[proci];

                if (proci != myProci && !recvData.empty())
                {
                    recvs.push_back
                    (
                        recvTuple
                        (
                            proci,
                            { recvData.data(), std::size_t(recvData.size()) }
                        )
                    );
                }
            }

            // Exchange buffers in chunks
            PstreamDetail::exchangeChunkedBuf<Type>
            (
                sends,
                recvs,
                tag,
                comm,
                wait
            );
        }
    }

    // Do myself. Already checked if in communicator
    recvBufs[myProci] = sendBufs[myProci];
}


template<class Container, class Type>
void Foam::Pstream::exchange
(
    const Map<Container>& sendBufs,
    const Map<label>& recvSizes,
    Map<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    static_assert(is_contiguous<Type>::value, "Contiguous data only!");

    const int myProci = UPstream::myProcNo(comm);

    // Initial: clear out receive 'slots'
    // Preferrable to clear out the map entries instead of the map itself
    // since this can potentially preserve allocated space
    // (eg DynamicList entries) between calls

    forAllIters(recvBufs, iter)
    {
        iter.val().clear();
    }

    if (UPstream::is_parallel(comm))
    {
        // Presize all receive buffers
        forAllIters(recvSizes, iter)
        {
            const label proci = iter.key();
            const label count = iter.val();

            if (proci != myProci && count > 0)
            {
                recvBufs(proci).resize_nocopy(count);
            }
        }

        // Define the exchange sequences as a flattened list.
        // We add an additional step of ordering the send/recv list
        // by message size, which can help with transfer speeds.

        typedef std::pair<int, stdFoam::span<const Type>> sendTuple;
        typedef std::pair<int, stdFoam::span<Type>> recvTuple;

        // Populate send sequences
        DynamicList<sendTuple> sends(sendBufs.size());
        forAllConstIters(sendBufs, iter)
        {
            const auto proci = iter.key();
            const auto& sendData = iter.val();

            if (proci != myProci && !sendData.empty())
            {
                sends.push_back
                (
                    sendTuple
                    (
                        proci,
                        { sendData.cdata(), std::size_t(sendData.size()) }
                    )
                );
            }
        }

        // Shorter messages first
        std::sort
        (
            sends.begin(),
            sends.end(),
            [=](const sendTuple& a, const sendTuple& b)
            {
                return (a.second.size() < b.second.size());
            }
        );

        // Populate recv sequences
        DynamicList<recvTuple> recvs(recvBufs.size());
        forAllIters(recvBufs, iter)
        {
            const auto proci = iter.key();
            auto& recvData = recvBufs[proci];

            if (proci != myProci && !recvData.empty())
            {
                recvs.push_back
                (
                    recvTuple
                    (
                        proci,
                        { recvData.data(), std::size_t(recvData.size()) }
                    )
                );
            }
        }

        // Shorter messages first
        std::sort
        (
            recvs.begin(),
            recvs.end(),
            [=](const recvTuple& a, const recvTuple& b)
            {
                return (a.second.size() < b.second.size());
            }
        );


        if (UPstream::maxCommsSize <= 0)
        {
            // Do the exchanging in a single go
            PstreamDetail::exchangeBuf<Type>
            (
                sends,
                recvs,
                tag,
                comm,
                wait
            );
        }
        else
        {
            // Exchange buffers in chunks
            PstreamDetail::exchangeChunkedBuf<Type>
            (
                sends,
                recvs,
                tag,
                comm,
                wait
            );
        }
    }

    // Do myself (if actually in the communicator)
    if (UPstream::is_rank(comm))
    {
        const auto iter = sendBufs.find(myProci);

        bool needsCopy = iter.good();

        if (needsCopy)
        {
            const auto& sendData = iter.val();

            needsCopy = !sendData.empty();
            if (needsCopy)
            {
                // insert_or_assign
                recvBufs(myProci) = sendData;
            }
        }

        if (!needsCopy)
        {
            recvBufs.erase(myProci);
        }
    }
}


template<class Container>
void Foam::Pstream::exchangeSizes
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const Container& sendBufs,
    labelList& recvSizes,
    const label tag,
    const label comm
)
{
    if (!UPstream::is_rank(comm))
    {
        recvSizes.clear();
        return;  // Process not in communicator
    }

    const label myProci = UPstream::myProcNo(comm);
    const label numProcs = UPstream::nProcs(comm);

    if (sendBufs.size() != numProcs)
    {
        FatalErrorInFunction
            << "Size of container " << sendBufs.size()
            << " does not equal the number of processors " << numProcs
            << Foam::abort(FatalError);
    }

    labelList sendSizes(numProcs);
    for (label proci = 0; proci < numProcs; ++proci)
    {
        sendSizes[proci] = sendBufs[proci].size();
    }

    recvSizes.resize_nocopy(numProcs);
    recvSizes = 0;  // Ensure non-received entries are properly zeroed

    // Preserve self-send, even if not described by neighbourhood
    recvSizes[myProci] = sendSizes[myProci];

    const label startOfRequests = UPstream::nRequests();

    for (const label proci : recvProcs)
    {
        if (proci != myProci)
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                reinterpret_cast<char*>(&recvSizes[proci]),
                sizeof(label),
                tag,
                comm
            );
        }
    }

    for (const label proci : sendProcs)
    {
        if (proci != myProci)
        {
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                reinterpret_cast<char*>(&sendSizes[proci]),
                sizeof(label),
                tag,
                comm
            );
        }
    }

    UPstream::waitRequests(startOfRequests);
}


template<class Container>
void Foam::Pstream::exchangeSizes
(
    const labelUList& neighProcs,
    const Container& sendBufs,
    labelList& recvSizes,
    const label tag,
    const label comm
)
{
    if (!UPstream::is_rank(comm))
    {
        recvSizes.clear();
        return;  // Process not in communicator
    }

    Pstream::exchangeSizes<Container>
    (
        neighProcs,  // send
        neighProcs,  // recv
        sendBufs,
        recvSizes,
        tag,
        comm
    );
}


// Sparse sending
template<class Container>
void Foam::Pstream::exchangeSizes
(
    const Map<Container>& sendBufs,
    Map<label>& recvSizes,
    const label tag,
    const label comm
)
{
    Map<label> sendSizes(2*sendBufs.size());
    recvSizes.clear();  // Done in allToAllConsensus too, but be explicit here

    if (!UPstream::is_rank(comm))
    {
        return;  // Process not in communicator
    }

    forAllConstIters(sendBufs, iter)
    {
        const label proci = iter.key();
        const label count = iter.val().size();

        if (count)
        {
            sendSizes.emplace(proci, count);
        }
    }

    UPstream::allToAllConsensus
    (
        sendSizes,
        recvSizes,
        (tag + 314159),  // some unique tag?
        comm
    );
}


template<class Container>
void Foam::Pstream::exchangeSizes
(
    const Container& sendBufs,
    labelList& recvSizes,
    const label comm
)
{
    if (!UPstream::is_rank(comm))
    {
        recvSizes.clear();
        return;  // Process not in communicator
    }

    const label numProcs = UPstream::nProcs(comm);

    if (sendBufs.size() != numProcs)
    {
        FatalErrorInFunction
            << "Size of container " << sendBufs.size()
            << " does not equal the number of processors " << numProcs
            << Foam::abort(FatalError);
    }

    labelList sendSizes(numProcs);
    forAll(sendBufs, proci)
    {
        sendSizes[proci] = sendBufs[proci].size();
    }
    recvSizes.resize_nocopy(sendSizes.size());

    if
    (
        UPstream::nProcsNonblockingExchange > 1
     && UPstream::nProcsNonblockingExchange <= numProcs
    )
    {
        // Use algorithm NBX: Nonblocking Consensus Exchange

        UPstream::allToAllConsensus
        (
            sendSizes,
            recvSizes,
            (UPstream::msgType() + 314159),  // some unique tag?
            comm
        );
        return;
    }

    UPstream::allToAll(sendSizes, recvSizes, comm);
}


template<class Container, class Type>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    // Algorithm PEX: Personalized Exchange
    // - Step 1: each process writes the data sizes to each peer and
    //   redistributes the vector (eg, MPI_Alltoall or non-blocking consensus)
    // - Step 2: size receive buffers and setup receives for all
    //   non-zero sendcounts. Post all sends and wait.

    labelList recvSizes;
    exchangeSizes(sendBufs, recvSizes, comm);

    exchange<Container, Type>(sendBufs, recvSizes, recvBufs, tag, comm, wait);
}


template<class Container, class Type>
void Foam::Pstream::exchange
(
    const Map<Container>& sendBufs,
    Map<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    // Algorithm PEX: Personalized Exchange
    // but using nonblocking consensus exchange for the sizes

    Map<label> recvSizes;
    exchangeSizes(sendBufs, recvSizes, tag, comm);

    exchange<Container, Type>(sendBufs, recvSizes, recvBufs, tag, comm, wait);
}


// ************************************************************************* //
