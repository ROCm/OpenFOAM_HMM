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

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "UPstreamWrapping.H"

#include <mpi.h>
#include <cinttypes>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialisations for bool

void Foam::reduce
(
    bool& value,
    const andOp<bool>&,
    const int tag,  /* (unused) */
    const label comm
)
{
    // This can also work:
    // PstreamDetail::allReduce(&value, 1, MPI_BYTE, MPI_BAND, comm);
    PstreamDetail::allReduce(&value, 1, MPI_C_BOOL, MPI_LAND, comm);
}


void Foam::reduce
(
    bool& value,
    const orOp<bool>&,
    const int tag,  /* (unused) */
    const label comm
)
{
    // This can also work:
    // PstreamDetail::allReduce(&value, 1, MPI_BYTE, MPI_BOR, comm);
    PstreamDetail::allReduce(&value, 1, MPI_C_BOOL, MPI_LOR, comm);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialisations for common reduction types

#undef  Pstream_CommonReductions
#define Pstream_CommonReductions(Native, TaggedType)                          \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const minOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_MIN, comm                                  \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const maxOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_MAX, comm                                  \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_SUM, comm                                  \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        values, size, TaggedType, MPI_SUM, comm                               \
    );                                                                        \
}                                                                             \


Pstream_CommonReductions(int32_t, MPI_INT32_T);
Pstream_CommonReductions(int64_t, MPI_INT64_T);
Pstream_CommonReductions(uint32_t, MPI_UINT32_T);
Pstream_CommonReductions(uint64_t, MPI_UINT64_T);
Pstream_CommonReductions(float,   MPI_FLOAT);
Pstream_CommonReductions(double,  MPI_DOUBLE);

#undef Pstream_CommonReductions


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialisations for floating-point types

#undef  Pstream_FloatReductions
#define Pstream_FloatReductions(Native, TaggedType)                           \
                                                                              \
void Foam::sumReduce                                                          \
(                                                                             \
    Native& value,                                                            \
    label& count,                                                             \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    if (UPstream::parRun())                                                   \
    {                                                                         \
        Native values[2];                                                     \
        values[0] = value;                                                    \
        values[1] = static_cast<Native>(count);                               \
                                                                              \
        PstreamDetail::allReduce<Native>                                      \
        (                                                                     \
            values, 2, TaggedType, MPI_SUM, comm                              \
        );                                                                    \
                                                                              \
        value = values[0];                                                    \
        count = static_cast<label>(values[1]);                                \
    }                                                                         \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm,                                                         \
    label& requestID                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_SUM, comm, &requestID                      \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm,                                                         \
    label& requestID                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        values, size, TaggedType, MPI_SUM, comm, &requestID                   \
    );                                                                        \
}


Pstream_FloatReductions(float, MPI_FLOAT);
Pstream_FloatReductions(double, MPI_DOUBLE);

#undef Pstream_FloatReductions


// ************************************************************************* //
