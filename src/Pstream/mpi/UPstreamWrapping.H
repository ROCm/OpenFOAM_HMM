/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

Namespace
    Foam::PstreamDetail

Description
    Some implementation details for Pstream and/or MPI.

InNamespace
    Foam::PstreamDetail

Description
    Functions to wrap MPI_Bcast, MPI_Allreduce, MPI_Iallreduce etc.

SourceFiles
    UPstreamWrappingTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstreamWrapping_H
#define Foam_UPstreamWrapping_H

#include "UPstream.H"
#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// MPI_Bcast, using root=0
template<class Type>
void broadcast0
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    const label comm
);

// MPI_Reduce, using root=0
template<class Type>
void reduce0
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label comm
);

// MPI_Allreduce or MPI_Iallreduce
template<class Type>
void allReduce
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label comm,                   // Communicator
    label* requestID = nullptr          // Non-null for MPI_Iallreduce
);


// MPI_Alltoall or MPI_Ialltoall with one element per rank
template<class Type>
void allToAll
(
    const UList<Type>& sendData,
    UList<Type>& recvData,
    MPI_Datatype datatype,
    const label comm,                   // Communicator
    label* requestID = nullptr          // Non-null for MPI_Ialltoall
);


// MPI_Alltoallv or MPI_Ialltoallv
template<class Type>
void allToAllv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    MPI_Datatype datatype,
    const label comm,                   // Communicator
    label* requestID = nullptr          // Non-null for MPI_Ialltoallv
);


// MPI_Gather or MPI_Igather
template<class Type>
void gather
(
    const Type* sendData,
    int sendCount,

    Type* recvData,                     // Ignored on non-root rank
    int recvCount,                      // Ignored on non-root rank

    MPI_Datatype datatype,              // The send/recv data type
    const label comm,                   // Communicator
    label* requestID = nullptr          // Non-null for MPI_Igather
);


// MPI_Scatter or MPI_Iscatter
template<class Type>
void scatter
(
    const Type* sendData,               // Ignored on non-root rank
    int sendCount,                      // Ignored on non-root rank

    Type* recvData,
    int recvCount,

    MPI_Datatype datatype,              // The send/recv data type
    const label comm,                   // Communicator
    label* requestID = nullptr          // Non-null for MPI_Iscatter
);


// MPI_Gatherv or MPI_Igatherv
template<class Type>
void gatherv
(
    const Type* sendData,
    int sendCount,              // Ignored on master if recvCounts[0] == 0

    Type* recvData,                     // Ignored on non-root rank
    const UList<int>& recvCounts,       // Ignored on non-root rank
    const UList<int>& recvOffsets,      // Ignored on non-root rank

    MPI_Datatype datatype,              // The send/recv data type
    const label comm,                   // Communicator
    label* requestID = nullptr          // Non-null for MPI_Igatherv
);


// MPI_Scatterv or MPI_Iscatterv
template<class Type>
void scatterv
(
    const Type* sendData,               // Ignored on non-root rank
    const UList<int>& sendCounts,       // Ignored on non-root rank
    const UList<int>& sendOffsets,      // Ignored on non-root rank

    Type* recvData,
    int recvCount,

    MPI_Datatype datatype,              // The send/recv data type
    const label comm,                   // Communicator
    label* requestID = nullptr          // Non-null for MPI_Igatherv
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PstreamDetail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPstreamWrappingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
