/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UIPBstream::bufferIPCrecv()
{
    // Uses double broadcast. Symmetric with UOPBstream::bufferIPCsend()
    // 1. for the data size
    // 2. for the data itself

    // Expected message size, similar to MPI_Probe
    // Same type must be expected in UOPBstream::bufferIPCsend()
    label bufSize(0);

    // Broadcast #1 - data size
    if
    (
        !UPstream::broadcast
        (
            reinterpret_cast<char*>(&bufSize),
            sizeof(label),
            comm_,
            fromProcNo_  //< is actually rootProcNo
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Bcast failure receiving buffer size" << nl
            << Foam::abort(FatalError);
    }

    if (UPstream::debug)
    {
        Pout<< "UOPBstream IPC read buffer :"
            << " root:" << fromProcNo_
            << " comm:" << comm_
            << " probed size:" << bufSize
            << " wanted size:" << recvBuf_.capacity()
            << Foam::endl;
    }

    // No buffer size allocated/specified
    if (!recvBuf_.capacity())
    {
        recvBuf_.resize(bufSize);
    }

    // This is the only real information we can trust
    messageSize_ = bufSize;

    // Broadcast #2 - data content
    // - skip if there is no data to receive

    if (messageSize_)
    {
        if
        (
            !UPstream::broadcast
            (
                recvBuf_.data(),
                messageSize_,  // same as bufSize
                comm_,
                fromProcNo_  //< is actually rootProcNo
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Bcast failure receiving buffer data:" << bufSize << nl
                << Foam::abort(FatalError);
        }
    }

    // Set addressed size. Leave actual allocated memory intact.
    recvBuf_.resize(messageSize_);

    if (!messageSize_)
    {
        setEof();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UIPBstream::read
(
    const int rootProcNo,
    char* buf,
    const std::streamsize bufSize,
    const label comm
)
{
    if
    (
        !UPstream::broadcast(buf, bufSize, comm, rootProcNo)
    )
    {
        FatalErrorInFunction
            << "MPI_Bcast failure receiving data:" << label(bufSize) << nl
            << Foam::abort(FatalError);
        return 0;
    }

    return bufSize;
}


// ************************************************************************* //
