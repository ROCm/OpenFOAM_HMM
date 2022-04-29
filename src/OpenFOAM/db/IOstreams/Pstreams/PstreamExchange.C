/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

Description
    Exchange data.

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "contiguous.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, class T>
void Foam::Pstream::exchangeContainer
(
    const UList<Container>& sendBufs,
    const labelUList& recvSizes,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    const label startOfRequests = UPstream::nRequests();

    // Set up receives
    // ~~~~~~~~~~~~~~~

    forAll(recvSizes, proci)
    {
        if (proci != Pstream::myProcNo(comm) && recvSizes[proci] > 0)
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                recvBufs[proci].data_bytes(),
                recvSizes[proci]*sizeof(T),
                tag,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    forAll(sendBufs, proci)
    {
        if (proci != Pstream::myProcNo(comm) && sendBufs[proci].size() > 0)
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendBufs[proci].cdata_bytes(),
                    sendBufs[proci].size_bytes(),
                    tag,
                    comm
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message. "
                    << "to:" << proci << " nBytes:"
                    << label(sendBufs[proci].size_bytes())
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


template<class T>
void Foam::Pstream::exchangeBuf
(
    const labelUList& sendSizes,
    const UList<const char*>& sendBufs,
    const labelUList& recvSizes,
    List<char*>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    const label startOfRequests = UPstream::nRequests();

    // Set up receives
    // ~~~~~~~~~~~~~~~

    forAll(recvSizes, proci)
    {
        if (proci != Pstream::myProcNo(comm) && recvSizes[proci] > 0)
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                recvBufs[proci],
                recvSizes[proci]*sizeof(T),
                tag,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    forAll(sendBufs, proci)
    {
        if (proci != Pstream::myProcNo(comm) && sendSizes[proci] > 0)
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendBufs[proci],
                    sendSizes[proci]*sizeof(T),
                    tag,
                    comm
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message. "
                    << "to:" << proci << " nBytes:"
                    << label(sendSizes[proci]*sizeof(T))
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
    // OR  static_assert(is_contiguous<T>::value, "Contiguous data only!")
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Contiguous data only." << sizeof(T) << Foam::abort(FatalError);
    }

    if (sendBufs.size() != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
            << "Size of list " << sendBufs.size()
            << " does not equal the number of processors "
            << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    recvBufs.resize_nocopy(sendBufs.size());

    if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        // Presize all receive buffers
        forAll(recvSizes, proci)
        {
            const label nRecv = recvSizes[proci];

            if (proci != Pstream::myProcNo(comm) && nRecv > 0)
            {
                recvBufs[proci].resize_nocopy(nRecv);
            }
        }

        if (UPstream::maxCommsSize <= 0)
        {
            // Do the exchanging in one go
            exchangeContainer<Container, T>
            (
                sendBufs,
                recvSizes,
                recvBufs,
                tag,
                comm,
                wait
            );
        }
        else
        {
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
                forAll(sendBufs, proci)
                {
                    if (proci != Pstream::myProcNo(comm))
                    {
                        nChunks = max(nChunks, sendBufs[proci].size());
                    }
                }

                // Convert from send count (elements) to number of chunks.
                // Can normally calculate with (count-1), but add some safety
                if (nChunks)
                {
                    nChunks = 1 + (nChunks/maxChunkSize);
                }
                reduce(nChunks, maxOp<label>(), tag, comm);
            }

            labelList nRecv(sendBufs.size());
            labelList nSend(sendBufs.size());
            labelList startRecv(sendBufs.size(), Zero);
            labelList startSend(sendBufs.size(), Zero);

            List<const char*> charPtrSend(sendBufs.size());
            List<char*> charPtrRecv(sendBufs.size());

            for (label iter = 0; iter < nChunks; ++iter)
            {
                forAll(sendBufs, proci)
                {
                    nSend[proci] = min
                    (
                        maxChunkSize,
                        sendBufs[proci].size()-startSend[proci]
                    );
                    nRecv[proci] = min
                    (
                        maxChunkSize,
                        recvBufs[proci].size()-startRecv[proci]
                    );

                    charPtrSend[proci] =
                    (
                        nSend[proci] > 0
                      ? reinterpret_cast<const char*>
                        (
                            &(sendBufs[proci][startSend[proci]])
                        )
                      : nullptr
                    );
                    charPtrRecv[proci] =
                    (
                        nRecv[proci] > 0
                      ? reinterpret_cast<char*>
                        (
                            &(recvBufs[proci][startRecv[proci]])
                        )
                      : nullptr
                    );
                }

                /// Info<< "iter " << iter
                ///     << ": beg=" << flatOutput(startSend)
                ///     << " len=" << flatOutput(nSend) << endl;

                exchangeBuf<T>
                (
                    nSend,
                    charPtrSend,
                    nRecv,
                    charPtrRecv,
                    tag,
                    comm,
                    wait
                );

                forAll(nSend, proci)
                {
                    startSend[proci] += nSend[proci];
                    startRecv[proci] += nRecv[proci];
                }
            }
        }
    }

    // Do myself
    recvBufs[Pstream::myProcNo(comm)] = sendBufs[Pstream::myProcNo(comm)];
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
    if (sendBufs.size() != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
            << "Size of container " << sendBufs.size()
            << " does not equal the number of processors "
            << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    labelList sendSizes(sendProcs.size());
    forAll(sendProcs, i)
    {
        sendSizes[i] = sendBufs[sendProcs[i]].size();
    }

    recvSizes.resize_nocopy(sendBufs.size());
    recvSizes = 0;  // Ensure non-received entries are properly zeroed

    const label startOfRequests = UPstream::nRequests();

    for (const label proci : recvProcs)
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

    forAll(sendProcs, i)
    {
        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            sendProcs[i],
            reinterpret_cast<char*>(&sendSizes[i]),
            sizeof(label),
            tag,
            comm
        );
    }

    UPstream::waitRequests(startOfRequests);
}


/// FUTURE?
///
/// template<class Container>
/// void Foam::Pstream::exchangeSizes
/// (
///     const labelUList& neighProcs,
///     const Container& sendBufs,
///     labelList& recvSizes,
///     const label tag,
///     const label comm
/// )
/// {
///     exchangeSizes<Container>(neighProcs, neighProcs, sendBufs, tag, comm);
/// }


template<class Container>
void Foam::Pstream::exchangeSizes
(
    const Container& sendBufs,
    labelList& recvSizes,
    const label comm
)
{
    if (sendBufs.size() != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
            << "Size of container " << sendBufs.size()
            << " does not equal the number of processors "
            << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    labelList sendSizes(sendBufs.size());
    forAll(sendBufs, proci)
    {
        sendSizes[proci] = sendBufs[proci].size();
    }
    recvSizes.resize_nocopy(sendSizes.size());
    UPstream::allToAll(sendSizes, recvSizes, comm);
}


template<class Container, class T>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    labelList recvSizes;
    exchangeSizes(sendBufs, recvSizes, comm);

    exchange<Container, T>(sendBufs, recvSizes, recvBufs, tag, comm, wait);
}


// ************************************************************************* //
