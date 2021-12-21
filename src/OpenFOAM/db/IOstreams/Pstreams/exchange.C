/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    const bool block
)
{
    const label startOfRequests = Pstream::nRequests();

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
                    << label(sendBufs[proci].size()*sizeof(T))
                    << Foam::abort(FatalError);
            }
        }
    }


    // Wait for all to finish
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (block)
    {
        Pstream::waitRequests(startOfRequests);
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
    const bool block
)
{
    const label startOfRequests = Pstream::nRequests();

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

    if (block)
    {
        Pstream::waitRequests(startOfRequests);
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
    const bool block
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

    recvBufs.setSize(sendBufs.size());

    if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        // Presize all receive buffers
        forAll(recvSizes, proci)
        {
            const label nRecv = recvSizes[proci];

            if (proci != Pstream::myProcNo(comm) && nRecv > 0)
            {
                recvBufs[proci].setSize(nRecv);
            }
        }

        if (Pstream::maxCommsSize <= 0)
        {
            // Do the exchanging in one go
            exchangeContainer<Container, T>
            (
                sendBufs,
                recvSizes,
                recvBufs,
                tag,
                comm,
                block
            );
        }
        else
        {
            // Determine the number of chunks to send. Note that we
            // only have to look at the sending data since we are
            // guaranteed that some processor's sending size is some other
            // processor's receive size. Also we can ignore any local comms.

            label maxNSend = 0;
            forAll(sendBufs, proci)
            {
                if (proci != Pstream::myProcNo(comm))
                {
                    maxNSend = max(maxNSend, sendBufs[proci].size());
                }
            }

            const label maxNBytes = sizeof(T)*maxNSend;

            // We need to send maxNBytes bytes so the number of iterations:
            //  maxNBytes                           iterations
            //  ---------                           ----------
            //  0                                   0
            //  1..maxCommsSize                     1
            //  maxCommsSize+1..2*maxCommsSize      2
            //      etc.

            label nIter;
            if (maxNBytes == 0)
            {
                nIter = 0;
            }
            else
            {
                nIter = (maxNBytes-1)/Pstream::maxCommsSize+1;
            }
            reduce(nIter, maxOp<label>(), tag, comm);


            List<const char*> charSendBufs(sendBufs.size());
            List<char*> charRecvBufs(sendBufs.size());

            labelList nRecv(sendBufs.size());
            labelList startRecv(sendBufs.size(), Zero);
            labelList nSend(sendBufs.size());
            labelList startSend(sendBufs.size(), Zero);

            for (label iter = 0; iter < nIter; iter++)
            {
                forAll(sendBufs, proci)
                {
                    nSend[proci] = min
                    (
                        Pstream::maxCommsSize,
                        sendBufs[proci].size()-startSend[proci]
                    );
                    charSendBufs[proci] =
                    (
                        nSend[proci] > 0
                      ? reinterpret_cast<const char*>
                        (
                            &(sendBufs[proci][startSend[proci]])
                        )
                      : nullptr
                    );

                    nRecv[proci] = min
                    (
                        Pstream::maxCommsSize,
                        recvBufs[proci].size()-startRecv[proci]
                    );

                    charRecvBufs[proci] =
                    (
                        nRecv[proci] > 0
                      ? reinterpret_cast<char*>
                        (
                            &(recvBufs[proci][startRecv[proci]])
                        )
                      : nullptr
                    );
                }

                exchangeBuf<T>
                (
                    nSend,
                    charSendBufs,
                    nRecv,
                    charRecvBufs,
                    tag,
                    comm,
                    block
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
    recvSizes.setSize(sendSizes.size());
    allToAll(sendSizes, recvSizes, comm);
}


template<class Container, class T>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool block
)
{
    labelList recvSizes;
    exchangeSizes(sendBufs, recvSizes, comm);

    exchange<Container, T>(sendBufs, recvSizes, recvBufs, tag, comm, block);
}


// ************************************************************************* //
