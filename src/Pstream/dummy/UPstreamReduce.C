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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialisations for bool

void Foam::reduce
(
    bool& value,
    const andOp<bool>&,
    const int tag,
    const label comm
)
{}

void Foam::reduce
(
    bool& value,
    const orOp<bool>&,
    const int tag,
    const label comm
)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialisations for common reduction types

#undef  Pstream_CommonReductions
#define Pstream_CommonReductions(Native)                                      \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const minOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const maxOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}


Pstream_CommonReductions(int32_t);
Pstream_CommonReductions(int64_t);
Pstream_CommonReductions(uint32_t);
Pstream_CommonReductions(uint64_t);
Pstream_CommonReductions(float);
Pstream_CommonReductions(double);

#undef Pstream_CommonReductions


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialisations for floating-point types

#undef  Pstream_FloatReductions
#define Pstream_FloatReductions(Native)                                       \
                                                                              \
void Foam::sumReduce                                                          \
(                                                                             \
    Native& value,                                                            \
    label& count,                                                             \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm,                                                         \
    label& requestID                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm,                                                         \
    label& requestID                                                          \
)                                                                             \
{}


Pstream_FloatReductions(float);
Pstream_FloatReductions(double);

#undef Pstream_FloatReductions


// ************************************************************************* //
