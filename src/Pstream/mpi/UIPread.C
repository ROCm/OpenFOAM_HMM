/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Read from UIPstream

\*---------------------------------------------------------------------------*/

#include "mpi.h"

#include "UIPstream.H"
#include "PstreamGlobals.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UIPstream::UIPstream
(
    const commsTypes commsType,
    const int fromProcNo,
    DynamicList<char>& externalBuf,
    const int tag,
    streamFormat format,
    versionNumber version
)
:
    UPstream(commsType),
    Istream(format, version),
    fromProcNo_(fromProcNo),
    externalBuf_(externalBuf),
    externalBufPosition_(0),
    tag_(tag),
    messageSize_(0)
{
    setOpened();
    setGood();

    if (commsType == UPstream::nonBlocking)
    {
        // Message is already received into externalBuf
    }
    else
    {
        MPI_Status status;

        label wantedSize = externalBuf_.capacity();

        // If the buffer size is not specified, probe the incomming message
        // and set it
        if (!wantedSize)
        {
            MPI_Probe(procID(fromProcNo_), tag_, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_BYTE, &messageSize_);

            externalBuf_.setCapacity(messageSize_);
            wantedSize = messageSize_;
        }

        messageSize_ = UIPstream::read
        (
            commsType,
            fromProcNo_,
            externalBuf_.begin(),
            wantedSize,
            tag_
        );

        // Set addressed size. Leave actual allocated memory intact.
        externalBuf_.setSize(messageSize_);

        if (!messageSize_)
        {
            FatalErrorIn
            (
                "UIPstream::UIPstream(const commsTypes, const int, "
                "DynamicList<char>&, streamFormat, versionNumber)"
            )   << "read failed"
                << Foam::abort(FatalError);
        }
    }
}


Foam::UIPstream::UIPstream(const int fromProcNo, PstreamBuffers& buffers)
:
    UPstream(buffers.commsType_),
    Istream(buffers.format_, buffers.version_),
    fromProcNo_(fromProcNo),
    externalBuf_(buffers.recvBuf_[fromProcNo]),
    externalBufPosition_(0),
    tag_(buffers.tag_),
    messageSize_(0)
{
    if (commsType() != UPstream::scheduled && !buffers.finishedSendsCalled_)
    {
        FatalErrorIn("UIPstream::UIPstream(const int, PstreamBuffers&)")
            << "PstreamBuffers::finishedSends() never called." << endl
            << "Please call PstreamBuffers::finishedSends() after doing"
            << " all your sends (using UOPstream) and before doing any"
            << " receives (using UIPstream)" << Foam::exit(FatalError);
    }

    setOpened();
    setGood();

    if (commsType() == UPstream::nonBlocking)
    {
        // Message is already received into externalBuf
    }
    else
    {
        MPI_Status status;

        label wantedSize = externalBuf_.capacity();

        // If the buffer size is not specified, probe the incomming message
        // and set it
        if (!wantedSize)
        {
            MPI_Probe(procID(fromProcNo_), tag_, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_BYTE, &messageSize_);

            externalBuf_.setCapacity(messageSize_);
            wantedSize = messageSize_;
        }

        messageSize_ = UIPstream::read
        (
            commsType(),
            fromProcNo_,
            externalBuf_.begin(),
            wantedSize,
            tag_
        );

        // Set addressed size. Leave actual allocated memory intact.
        externalBuf_.setSize(messageSize_);

        if (!messageSize_)
        {
            FatalErrorIn
            (
                "UIPstream::UIPstream(const int, PstreamBuffers&)"
            )   << "read failed"
                << Foam::abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UIPstream::read
(
    const commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize,
    const int tag
)
{
    if (commsType == blocking || commsType == scheduled)
    {
        MPI_Status status;

        if
        (
            MPI_Recv
            (
                buf,
                bufSize,
                MPI_PACKED,
                procID(fromProcNo),
                tag,
                MPI_COMM_WORLD,
                &status
            )
        )
        {
            FatalErrorIn
            (
                "UIPstream::read"
                "(const int fromProcNo, char* buf, std::streamsize bufSize)"
            )   << "MPI_Recv cannot receive incomming message"
                << Foam::abort(FatalError);

            return 0;
        }


        // Check size of message read

        label messageSize;
        MPI_Get_count(&status, MPI_BYTE, &messageSize);

        if (messageSize > bufSize)
        {
            FatalErrorIn
            (
                "UIPstream::read"
                "(const int fromProcNo, char* buf, std::streamsize bufSize)"
            )   << "buffer (" << label(bufSize)
                << ") not large enough for incomming message ("
                << messageSize << ')'
                << Foam::abort(FatalError);
        }

        return messageSize;
    }
    else if (commsType == nonBlocking)
    {
        MPI_Request request;

        if
        (
            MPI_Irecv
            (
                buf,
                bufSize,
                MPI_PACKED,
                procID(fromProcNo),
                tag,
                MPI_COMM_WORLD,
                &request
            )
        )
        {
            FatalErrorIn
            (
                "UIPstream::read"
                "(const int fromProcNo, char* buf, std::streamsize bufSize)"
            )   << "MPI_Recv cannot start non-blocking receive"
                << Foam::abort(FatalError);

            return 0;
        }

        PstreamGlobals::outstandingRequests_.append(request);

        // Assume the message is completely received.
        return bufSize;
    }
    else
    {
        FatalErrorIn
        (
            "UIPstream::read"
            "(const int fromProcNo, char* buf, std::streamsize bufSize)"
        )   << "Unsupported communications type " << commsType
            << Foam::abort(FatalError);

        return 0;
    }
}


// ************************************************************************* //
