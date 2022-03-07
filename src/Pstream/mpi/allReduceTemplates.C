/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "allReduce.H"
#include "profilingPstream.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PstreamDetail::broadcast0
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    const label communicator
)
{
    if (!UPstream::parRun())
    {
        return;
    }

    profilingPstream::beginTiming();

    // const int retval =
    MPI_Bcast
    (
        values,
        count,
        datatype,
        0,  // (root process) is master == UPstream::masterNo()
        PstreamGlobals::MPICommunicators_[communicator]
    );

    profilingPstream::addBroadcastTime();
}


template<class Type>
void Foam::PstreamDetail::reduce0
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label communicator
)
{
    if (!UPstream::parRun())
    {
        return;
    }

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:";
        if (count == 1)
        {
            Pout<< (*values);
        }
        else
        {
            Pout<< UList<Type>(values, count);
        }
        Pout<< " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Pout);
    }

    profilingPstream::beginTiming();

    // const int retval =
    MPI_Reduce
    (
        MPI_IN_PLACE,
        values,
        count,
        datatype,
        optype,
        0,  // (root process) is master == UPstream::masterNo()
        PstreamGlobals::MPICommunicators_[communicator]
    );

    profilingPstream::addReduceTime();
}


template<class Type>
void Foam::PstreamDetail::allReduce
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label communicator
)
{
    if (!UPstream::parRun())
    {
        return;
    }

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:";
        if (count == 1)
        {
            Pout<< (*values);
        }
        else
        {
            Pout<< UList<Type>(values, count);
        }
        Pout<< " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Pout);
    }

    profilingPstream::beginTiming();

    // const int retval =
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        values,
        count,
        datatype,
        optype,
        PstreamGlobals::MPICommunicators_[communicator]
    );

    profilingPstream::addReduceTime();
}


template<class Type>
void Foam::PstreamDetail::iallReduce
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label communicator,
    label& requestID
)
{
    if (!UPstream::parRun())
    {
        return;
    }

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** non-blocking reducing:";
        if (count == 1)
        {
            Pout<< (*values);
        }
        else
        {
            Pout<< UList<Type>(values, count);
        }
        Pout<< " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Pout);
    }

    profilingPstream::beginTiming();

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    MPI_Request request;
    if
    (
        MPI_Iallreduce
        (
            MPI_IN_PLACE,
            values,
            count,
            datatype,
            optype,
            PstreamGlobals::MPICommunicators_[communicator],
            &request
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Iallreduce failed for "
            << UList<Type>(values, count)
            << Foam::abort(FatalError);
    }

    if (PstreamGlobals::freedRequests_.size())
    {
        requestID = PstreamGlobals::freedRequests_.remove();
        PstreamGlobals::outstandingRequests_[requestID] = request;
    }
    else
    {
        requestID = PstreamGlobals::outstandingRequests_.size();
        PstreamGlobals::outstandingRequests_.append(request);
    }

    if (UPstream::debug)
    {
        Pout<< "UPstream::allocateRequest for non-blocking reduce"
            << " : request:" << requestID << endl;
    }
#else
    // Non-blocking not yet implemented in mpi
    if
    (
        MPI_Allreduce
        (
            MPI_IN_PLACE,
            values,
            count,
            datatype,
            optype,
            PstreamGlobals::MPICommunicators_[communicator]
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Allreduce failed for "
            << UList<Type>(values, count)
            << Foam::abort(FatalError);
    }
    requestID = -1;
#endif

    profilingPstream::addReduceTime();
}


// ************************************************************************* //
