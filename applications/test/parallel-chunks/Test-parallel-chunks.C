/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Test-parallel-chunks

Description
    Test for sending contiguous data in chunk-wise.
    Largely mirrors Pstream::exchange or vice versa

\*---------------------------------------------------------------------------*/

#define Foam_PstreamExchange_debug_chunks

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"

using namespace Foam;


// Looks like Pstream::exchangeBuf
template<class T>
void do_exchangeBuf
(
    const label sendSize,
    const char* sendData,
    const label recvSize,
    char* recvData,
    const int tag,
    const label comm,
    const bool wait
)
{
    const label startOfRequests = UPstream::nRequests();

    // Set up receives
    // ~~~~~~~~~~~~~~~

    // forAll(recvSizes, proci)
    {
        // if (proci != Pstream::myProcNo(comm) && recvSizes[proci] > 0)
        if (!Pstream::master(comm) && recvSize > 0)
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                UPstream::myProcNo(comm),   // proci,
                recvData,
                recvSize*sizeof(T),
                tag,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    // forAll(sendBufs, proci)
    for (const int proci : Pstream::subProcs(comm))
    {
        if (sendSize > 0)
        // if (proci != Pstream::myProcNo(comm) && sendSizes[proci] > 0)
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendData,
                    sendSize*sizeof(T),
                    tag,
                    comm
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message. "
                    << "to:" << proci << " nBytes:"
                    << label(sendSize*sizeof(T))
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


// Looks like Pstream::exchangeContainer
template<class Container, class T>
void do_exchangeContainer
(
    const Container& sendData,
    const label recvSize,
    Container& recvData,
    const int tag,
    const label comm,
    const bool wait
)
{
    const label startOfRequests = UPstream::nRequests();

    // Set up receives
    // ~~~~~~~~~~~~~~~

    // for (const int proci : Pstream::allProcs(comm))
    {
        if (!Pstream::master(comm) && recvSize > 0)
        // if (proci != Pstream::myProcNo(comm) && recvSize > 0)
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                UPstream::myProcNo(comm),  // proci,
                recvData.data_bytes(),
                recvSize*sizeof(T),
                tag,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    if (Pstream::master(comm) && sendData.size() > 0)
    {
        for (const int proci : Pstream::subProcs(comm))
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


template<class Container, class T>
void broadcast_chunks
(
    Container& sendData,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm,
    const bool wait = true
)
{
    // OR  static_assert(is_contiguous<T>::value, "Contiguous data only!")
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Contiguous data only." << sizeof(T) << Foam::abort(FatalError);
    }

    if (UPstream::maxCommsSize <= 0)
    {
        // Do in one go
        Info<< "send " << sendData.size() << " elements in one go" << endl;
        Pstream::broadcast(sendData, comm);
        return;
    }

    label sendSize(sendData.size());
    Pstream::broadcast(sendSize, comm);

    label recvSize(sendSize);

    sendData.resize_nocopy(recvSize);  // A no-op on master

    // Determine the number of chunks to send. Note that we
    // only have to look at the sending data since we are
    // guaranteed that some processor's sending size is some other
    // processor's receive size. Also we can ignore any local comms.

    // We need to send chunks so the number of iterations:
    //  maxChunkSize                        iterations
    //  ------------                        ----------
    //  0                                   0
    //  1..maxChunkSize                     1
    //  maxChunkSize+1..2*maxChunkSize      2
    //  ...

    const label maxChunkSize
    (
        max
        (
            static_cast<label>(1),
            static_cast<label>(UPstream::maxCommsSize/sizeof(T))
        )
    );

    label nChunks(0);
    {
        // Get max send count (elements)
        // forAll(sendBufs, proci)
        // {
        //     if (proci != Pstream::myProcNo(comm))
        //     {
        //         nChunks = max(nChunks, sendBufs[proci].size());
        //     }
        // }
        nChunks = sendSize;

        // Convert from send count (elements) to number of chunks.
        // Can normally calculate with (count-1), but add some safety
        if (nChunks)
        {
            nChunks = 1 + (nChunks/maxChunkSize);
        }
        reduce(nChunks, maxOp<label>(), tag, comm);

        Info
            << "send " << sendSize << " elements ("
            << (sendSize*sizeof(T)) << " bytes) in " << nChunks
            << " chunks of " << maxChunkSize << " elements ("
            << (maxChunkSize*sizeof(T)) << " bytes) for maxCommsSize:"
            << Pstream::maxCommsSize
            << endl;
    }

    // stress-test with shortened sendSize
    // will produce useless loops, but no calls
    // sendSize /= 2;

    label nSend(0);
    label startSend(0);
    char* charPtrSend;

    for (label iter = 0; iter < nChunks; ++iter)
    {
        nSend = min
        (
            maxChunkSize,
            sendSize-startSend
        );

        charPtrSend =
        (
            nSend > 0
          ? reinterpret_cast<char*>(&(sendData[startSend]))
          : nullptr
        );

        Info<< "iter " << iter
            << ": beg=" << startSend << " len=" << nSend
            << " (" << (nSend*sizeof(T)) << " bytes)" << endl;

        UPstream::broadcast(charPtrSend, nSend*sizeof(T), comm);

        // forAll(nSend, proci)
        {
            startSend += nSend;
        }
    }

    Info<< "final: " << startSend << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addOption("comms-size", "int", "override Pstream::maxCommsSize");

    #include "setRootCase.H"

    if (!Pstream::parRun())
    {
        Info<< "\nWarning: not parallel - skipping further tests\n" << endl;
        return 0;
    }

    labelList input1;
    if (Pstream::master())
    {
        input1 = identity(500);
    }
    broadcast_chunks<labelList, label>(input1);

    Pstream::maxCommsSize = 33;

    args.readIfPresent("comms-size", Pstream::maxCommsSize);

    broadcast_chunks<labelList, label>(input1);

    // Mostly the same with PstreamBuffers
    if (false)
    {
        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

        labelList sendData;
        if (Pstream::master())
        {
            sendData = identity(500);

            for (const int proci : Pstream::subProcs())
            {
                UOPstream os(proci, pBufs);
                os << sendData;
            }
        }

        Info<< "call finishedSends()" << endl;
        pBufs.finishedScatters();

        if (!Pstream::master())
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> sendData;
        }
    }

    // Manually
    Info<< "perform list exchange" << endl;
    {
        labelListList sendBufs(UPstream::nProcs());
        labelListList recvBufs(UPstream::nProcs());
        labelList recvSizes;

        if (Pstream::master())
        {
            for (const int proci : Pstream::allProcs())
            {
                if (proci != Pstream::myProcNo())
                {
                    sendBufs[proci] = identity(500);
                }
            }
        }

        Pstream::exchangeSizes(sendBufs, recvSizes);

        Pstream::exchange<labelList, label>
        (
            sendBufs,
            recvSizes,
            recvBufs
        );
    }


    Info<< "perform Map exchange" << endl;
    {
        Map<labelList> sendBufs;
        Map<labelList> recvBufs;
        Map<label> recvSizes;

        if (Pstream::master())
        {
            for (const int proci : Pstream::allProcs())
            {
                if (proci != Pstream::myProcNo())
                {
                    sendBufs(proci) = identity(500);
                }
            }
        }

        Pstream::exchangeSizes(sendBufs, recvSizes);

        Pstream::exchange<labelList, label>
        (
            sendBufs,
            recvSizes,
            recvBufs
        );
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
