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
    Read token and binary block from IPstream

\*---------------------------------------------------------------------------*/

#include "mpi.h"

#include "IPstream.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Outstanding non-blocking operations.
//! @cond fileScope
//Foam::DynamicList<MPI_Request> IPstream_outstandingRequests_;
//! @endcond fileScope

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::IPstream::IPstream
(
    const commsTypes commsType,
    const int fromProcNo,
    const label bufSize,
    streamFormat format,
    versionNumber version
)
:
    Pstream(commsType, bufSize),
    Istream(format, version),
    fromProcNo_(fromProcNo),
    messageSize_(0)
{
    setOpened();
    setGood();

    MPI_Status status;

    // Cannot use buf_.size() since appends a few bytes extra
    label realBufSize = bufSize;

    // If the buffer size is not specified, probe the incomming message
    // and set it
    if (!bufSize)
    {
        if (commsType == nonBlocking)
        {
            FatalErrorIn
            (
                "IPstream::IPstream(const commsTypes, const int, "
                "const label, streamFormat, versionNumber)"
            )   << "Can use nonBlocking mode only with pre-allocated buffers"
                << Foam::abort(FatalError);
        }

        MPI_Probe(procID(fromProcNo_), msgType(), MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_BYTE, &messageSize_);

        buf_.setSize(messageSize_);
        realBufSize = buf_.size();
    }

    messageSize_ = read(commsType, fromProcNo_, buf_.begin(), realBufSize);

    if (!messageSize_)
    {
        FatalErrorIn
        (
            "IPstream::IPstream(const commsTypes, const int, "
            "const label, streamFormat, versionNumber)"
        )   << "read failed"
            << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::IPstream::read
(
    const commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize
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
                msgType(),
                MPI_COMM_WORLD,
                &status
            )
        )
        {
            FatalErrorIn
            (
                "IPstream::read"
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
                "IPstream::read"
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
                msgType(),
                MPI_COMM_WORLD,
                &request
            )
        )
        {
            FatalErrorIn
            (
                "IPstream::read"
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
            "IPstream::read"
            "(const int fromProcNo, char* buf, std::streamsize bufSize)"
        )   << "Unsupported communications type " << commsType
            << Foam::abort(FatalError);

        return 0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
