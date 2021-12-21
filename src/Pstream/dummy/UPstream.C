/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{}


bool Foam::UPstream::initNull()
{
    WarningInFunction
        << "The dummy Pstream library cannot be used in parallel mode"
        << endl;

    return false;
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    FatalErrorInFunction
        << "The dummy Pstream library cannot be used in parallel mode"
        << endl
        << Foam::exit(FatalError);

    return false;
}


void Foam::UPstream::shutdown(int errNo)
{}


void Foam::UPstream::exit(int errNo)
{
    // No MPI - just exit
    std::exit(errNo);
}


void Foam::UPstream::abort()
{
    // No MPI - just abort
    std::abort();
}


void Foam::reduce(scalar&, const sumOp<scalar>&, const int, const label)
{}


void Foam::reduce(scalar&, const minOp<scalar>&, const int, const label)
{}


void Foam::reduce(vector2D&, const sumOp<vector2D>&, const int, const label)
{}


void Foam::sumReduce
(
    scalar&,
    label&,
    const int,
    const label
)
{}


void Foam::reduce(scalar&, const sumOp<scalar>&, const int, const label, label&)
{}


void Foam::reduce
(
    scalar[],
    const int,
    const sumOp<scalar>&,
    const int,
    const label,
    label&
)
{}


#if defined(WM_SPDP)
void Foam::reduce
(
    solveScalar& Value,
    const sumOp<solveScalar>& bop,
    const int tag,
    const label comm
)
{}
void Foam::reduce
(
    solveScalar& Value,
    const minOp<solveScalar>& bop,
    const int tag,
    const label comm
)
{}
void Foam::reduce
(
    Vector2D<solveScalar>& Value,
    const sumOp<Vector2D<solveScalar>>& bop,
    const int tag,
    const label comm
)
{}
void Foam::sumReduce
(
    solveScalar& Value,
    label& Count,
    const int tag,
    const label comm
)
{}
void Foam::reduce
(
    solveScalar& Value,
    const sumOp<solveScalar>& bop,
    const int tag,
    const label comm,
    label& request
)
{}
void Foam::reduce
(
    solveScalar[],
    const int,
    const sumOp<solveScalar>&,
    const int,
    const label,
    label&
)
{}
#endif



void Foam::UPstream::allToAll
(
    const labelUList& sendData,
    labelUList& recvData,
    const label communicator
)
{
    recvData.deepCopy(sendData);
}


void Foam::UPstream::mpiGather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    std::memmove(recvData, sendData, sendSize);
}


void Foam::UPstream::mpiScatter
(
    const char* sendData,
    int sendSize,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    std::memmove(recvData, sendData, sendSize);
}


void Foam::UPstream::gather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,
    const label communicator
)
{
    std::memmove(recvData, sendData, sendSize);
}


void Foam::UPstream::scatter
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    std::memmove(recvData, sendData, recvSize);
}


void Foam::UPstream::allocatePstreamCommunicator
(
    const label,
    const label
)
{}


void Foam::UPstream::freePstreamCommunicator(const label)
{}


Foam::label Foam::UPstream::nRequests()
{
    return 0;
}


void Foam::UPstream::resetRequests(const label i)
{}


void Foam::UPstream::waitRequests(const label start)
{}


void Foam::UPstream::waitRequest(const label i)
{}


bool Foam::UPstream::finishedRequest(const label i)
{
    NotImplemented;
    return false;
}


// ************************************************************************* //
