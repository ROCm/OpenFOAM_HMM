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
    Functions to wrap MPI_Bcast, MPI_Allreduce, MPI_Iallreduce

SourceFiles
    allReduceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_allReduce_H
#define Foam_allReduce_H

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
    const label communicator
);

// MPI_Reduce, using root=0
template<class Type>
void reduce0
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label communicator
);

// MPI_Allreduce
template<class Type>
void allReduce
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label communicator
);


// MPI_Iallreduce
template<class Type>
void iallReduce
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label communicator,
    label& requestID
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PstreamDetail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "allReduceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
