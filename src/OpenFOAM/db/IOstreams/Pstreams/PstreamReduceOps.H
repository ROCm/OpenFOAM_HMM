/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Inter-processor communication reduction functions.

\*---------------------------------------------------------------------------*/

#ifndef Foam_PstreamReduceOps_H
#define Foam_PstreamReduceOps_H

#include "ops.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Reduce inplace (cf. MPI Allreduce)
//- using specified communication schedule.
template<class T, class BinaryOp>
void reduce
(
    const List<UPstream::commsStruct>& comms,
    T& value,
    const BinaryOp& bop,
    const int tag,
    const label comm
)
{
    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        Pout<< "** reducing:" << value << " with comm:" << comm << endl;
        error::printStack(Pout);
    }
    Pstream::gather(comms, value, bop, tag, comm);
    Pstream::broadcast(value, comm);
}


//- Reduce inplace (cf. MPI Allreduce)
//- using linear/tree communication schedule
template<class T, class BinaryOp>
void reduce
(
    T& value,
    const BinaryOp& bop,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    if (UPstream::parRun())
    {
        reduce(UPstream::whichCommunication(comm), value, bop, tag, comm);
    }
}


//- Reduce (copy) and return value
template<class T, class BinaryOp>
T returnReduce
(
    const T& value,
    const BinaryOp& bop,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    T work(value);
    reduce(work, bop, tag, comm);
    return work;
}


//- Reduce inplace (cf. MPI Allreduce)
//- the sum of both value and count (for averaging)
template<class T>
void sumReduce
(
    T& value,
    label& count,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    if (UPstream::parRun())
    {
        reduce(value, sumOp<T>(), tag, comm);
        reduce(count, sumOp<label>(), tag, comm);
    }
}


//- Reduce inplace (cf. MPI Allreduce)
//- multiple values (identical size on all processes!)
template<class T, class BinaryOp>
void reduce
(
    T values[],
    const int size,
    const BinaryOp&,
    const int tag,
    const label comm
)
{
    NotImplemented;
}

//- Non-blocking reduce inplace (cf. MPI Iallreduce)
//- single value. Sets request.
template<class T, class BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp&,
    const int tag,
    const label comm,
    label& request
)
{
    NotImplemented;
}

//- Non-blocking reduce inplace (cf. MPI Iallreduce)
//- of multiple values (same size on all processes!)
//  Sets request.
template<class T, class BinaryOp>
void reduce
(
    T values[],
    const int size,
    const BinaryOp&,
    const int tag,
    const label comm,
    label& request
)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialisations for bool

//- Logical (and) reduction (cf. MPI AllReduce)
void reduce
(
    bool& value,
    const andOp<bool>&,
    const int tag = UPstream::msgType(),  /*!< (ignored) */
    const label comm = UPstream::worldComm
);

//- Logical (or) reduction (cf. MPI AllReduce)
void reduce
(
    bool& value,
    const orOp<bool>&,
    const int tag = UPstream::msgType(),  /*!< (ignored) */
    const label comm = UPstream::worldComm
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialisations for common reduction types

#undef  Pstream_CommonReductions
#define Pstream_CommonReductions(Native)                                      \
                                                                              \
/*! \brief Reduce (min) single Native value */                                \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const minOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Reduce (max) single Native value */                                \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const maxOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Reduce (sum) single Native value */                                \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Reduce (sum) multiple Native values (identical size all procs!) */ \
void reduce                                                                   \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm                                                          \
);                                                                            \
                                                                              \
/*! \brief Reduce (sum) multiple Native values */                             \
template<unsigned N>                                                          \
inline void reduce                                                            \
(                                                                             \
    FixedList<Native, N>& values,                                             \
    const sumOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
)                                                                             \
{                                                                             \
    reduce(values.data(), int(values.size()), sumOp<Native>(), tag, comm);    \
}


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
/*! \brief Sum of both Native value and count (for averaging) */              \
void sumReduce                                                                \
(                                                                             \
    Native& value,                                                            \
    label& count,                                                             \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Non-blocking reduce (sum) single Native value. Sets request */     \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm,                                                         \
    label& requestID                                                          \
);                                                                            \
                                                                              \
/*! \brief Non-blocking reduce (sum) multiple Native values. Sets request */  \
void reduce                                                                   \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm,                                                         \
    label& requestID                                                          \
);


Pstream_FloatReductions(float);
Pstream_FloatReductions(double);

#undef Pstream_FloatReductions


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
