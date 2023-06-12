/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "UIPstream.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"
#include "IOstreams.H"

// FUTURE? probe and receive message
// - as of 2023-06 appears to be broken with INTELMPI + PMI-2 (slurm)
//   and perhaps other places so currently avoid

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// General blocking/non-blocking MPI receive, optionally with probed
// message information.
static Foam::label UPstream_mpi_receive
(
    const Foam::UPstream::commsTypes commsType,
    char* buf,
    const std::streamsize bufSize,
    const int fromProcNo,
    const int tag,
    const Foam::label communicator,
    Foam::UPstream::Request* req
)
{
    using namespace Foam;

    PstreamGlobals::reset_request(req);

    if (UPstream::debug)
    {
        Pout<< "UIPstream::read : starting read from:" << fromProcNo
            << " tag:" << tag << " comm:" << communicator
            << " wanted size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << Foam::endl;
    }
    if (UPstream::warnComm >= 0 && communicator != UPstream::warnComm)
    {
        Pout<< "UIPstream::read : starting read from:" << fromProcNo
            << " tag:" << tag << " comm:" << communicator
            << " wanted size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Pout);
    }

    profilingPstream::beginTiming();

    if
    (
        commsType == UPstream::commsTypes::blocking
     || commsType == UPstream::commsTypes::scheduled
    )
    {
        int returnCode = 0;
        MPI_Status status;

        {
            returnCode = MPI_Recv
            (
                buf,
                bufSize,
                MPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &status
            );
        }

        if (returnCode != MPI_SUCCESS)
        {
            FatalErrorInFunction
                << "MPI_Recv cannot receive incoming message"
                << Foam::abort(FatalError);
            return 0;
        }

        profilingPstream::addGatherTime();

        // Check size of message read

        int messageSize;
        MPI_Get_count(&status, MPI_BYTE, &messageSize);

        if (UPstream::debug)
        {
            Pout<< "UIPstream::read : finished read from:" << fromProcNo
                << " tag:" << tag << " read size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << Foam::endl;
        }

        if (messageSize > bufSize)
        {
            FatalErrorInFunction
                << "buffer (" << label(bufSize)
                << ") not large enough for incoming message ("
                << messageSize << ')'
                << Foam::abort(FatalError);
        }

        return messageSize;
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        int returnCode = 0;
        MPI_Request request;

        {
            returnCode = MPI_Irecv
            (
                buf,
                bufSize,
                MPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &request
            );
        }

        if (returnCode != MPI_SUCCESS)
        {
            FatalErrorInFunction
                << "MPI_Irecv cannot start non-blocking receive"
                << Foam::abort(FatalError);

            return 0;
        }

        if (UPstream::debug)
        {
            Pout<< "UIPstream::read : started read from:" << fromProcNo
                << " tag:" << tag << " read size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << " request:" <<
                (req ? label(-1) : PstreamGlobals::outstandingRequests_.size())
                << Foam::endl;
        }

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();

        // Assume the message will be completely received.
        return bufSize;
    }

    FatalErrorInFunction
        << "Unsupported communications type " << int(commsType)
        << Foam::abort(FatalError);

    return 0;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UIPstream::bufferIPCrecv()
{
    // Called by constructor
    if (UPstream::debug)
    {
        Pout<< "UIPstream IPC read buffer :"
            << " from:" << fromProcNo_
            << " tag:" << tag_ << " comm:" << comm_
            << " wanted size:" << recvBuf_.capacity()
            << Foam::endl;
    }

    // No buffer size allocated/specified - probe size of incoming message
    if (!recvBuf_.capacity())
    {
        profilingPstream::beginTiming();

        MPI_Status status;

        MPI_Probe
        (
            fromProcNo_,
            tag_,
            PstreamGlobals::MPICommunicators_[comm_],
           &status
        );

        MPI_Get_count(&status, MPI_BYTE, &messageSize_);

        profilingPstream::addProbeTime();

        recvBuf_.resize(messageSize_);

        if (UPstream::debug)
        {
            Pout<< "UIPstream::UIPstream : probed size:"
                << messageSize_ << Foam::endl;
        }
    }

    messageSize_ = UPstream_mpi_receive
    (
        commsType(),
        recvBuf_.data(),
        recvBuf_.capacity(),
        fromProcNo_,
        tag_,
        comm_,
        nullptr   // UPstream::Request
    );

    // Set addressed size. Leave actual allocated memory intact.
    recvBuf_.resize(messageSize_);

    if (!messageSize_)
    {
        setEof();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UIPstream::read
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize,
    const int tag,
    const label communicator,
    UPstream::Request* req
)
{
    return UPstream_mpi_receive
    (
        commsType,
        buf,
        bufSize,
        fromProcNo,
        tag,
        communicator,
        req
    );
}


// ************************************************************************* //
