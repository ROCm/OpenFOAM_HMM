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

#include "UPstreamWrapping.H"
#include "profilingPstream.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PstreamDetail::broadcast0
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    const label comm
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
        0,  // (root rank) == UPstream::masterNo()
        PstreamGlobals::MPICommunicators_[comm]
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
    const label comm
)
{
    if (!UPstream::parRun())
    {
        return;
    }

    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
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
        Pout<< " with comm:" << comm
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Pout);
    }

    profilingPstream::beginTiming();

    // const int retval =
    MPI_Reduce
    (
        MPI_IN_PLACE,  // recv is also send
        values,
        count,
        datatype,
        optype,
        0,  // (root rank) == UPstream::masterNo()
        PstreamGlobals::MPICommunicators_[comm]
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
    const label comm,
    label* requestID
)
{
    if (!UPstream::parRun())
    {
        return;
    }

    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        if (requestID != nullptr)
        {
            Pout<< "** MPI_Iallreduce (non-blocking):";
        }
        else
        {
            Pout<< "** MPI_Allreduce (blocking):";
        }
        if (count == 1)
        {
            Pout<< (*values);
        }
        else
        {
            Pout<< UList<Type>(values, count);
        }
        Pout<< " with comm:" << comm
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Pout);
    }

    profilingPstream::beginTiming();

    bool handled(false);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (requestID != nullptr)
    {
        handled = true;
        MPI_Request request;
        if
        (
            MPI_Iallreduce
            (
                MPI_IN_PLACE,  // recv is also send
                values,
                count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[comm],
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
            *requestID = PstreamGlobals::freedRequests_.remove();
            PstreamGlobals::outstandingRequests_[*requestID] = request;
        }
        else
        {
            *requestID = PstreamGlobals::outstandingRequests_.size();
            PstreamGlobals::outstandingRequests_.append(request);
        }
    }
#endif

    if (!handled)
    {
        if (requestID != nullptr)
        {
            *requestID = -1;
        }
        if
        (
            MPI_Allreduce
            (
                MPI_IN_PLACE,  // recv is also send
                values,
                count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Allreduce failed for "
                << UList<Type>(values, count)
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addReduceTime();
}


template<class Type>
void Foam::PstreamDetail::allToAll
(
    const UList<Type>& sendData,
    UList<Type>& recvData,
    MPI_Datatype datatype,
    const label comm,
    label* requestID
)
{
    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        if (requestID != nullptr)
        {
            Pout<< "** MPI_Ialltoall (non-blocking):";
        }
        else
        {
            Pout<< "** MPI_Alltoall (blocking):";
        }
        Pout<< " np:" << np
            << " sendData:" << sendData.size()
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if (sendData.size() != np || recvData.size() != np)
    {
        FatalErrorInFunction
            << "Have " << np << " ranks, but size of sendData:"
            << sendData.size() << " or recvData:" << recvData.size()
            << " is different!"
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        recvData.deepCopy(sendData);
        return;
    }

    profilingPstream::beginTiming();

    bool handled(false);

    #if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (requestID != nullptr)
    {
        handled = true;
        MPI_Request request;
        if
        (
            MPI_Ialltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                const_cast<Type*>(sendData.cdata()),
                1,                      // one element per rank
                datatype,
                recvData.data(),
                1,                      // one element per rank
                datatype,
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Ialltoall [comm: " << comm << "] failed."
                << " For " << sendData
                << Foam::abort(FatalError);
        }

        if (PstreamGlobals::freedRequests_.size())
        {
            *requestID = PstreamGlobals::freedRequests_.remove();
            PstreamGlobals::outstandingRequests_[*requestID] = request;
        }
        else
        {
            *requestID = PstreamGlobals::outstandingRequests_.size();
            PstreamGlobals::outstandingRequests_.append(request);
        }
    }
#endif


    if (!handled)
    {
        if (requestID != nullptr)
        {
            *requestID = -1;
        }
        if
        (
            MPI_Alltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                const_cast<Type*>(sendData.cdata()),
                1,                      // one element per rank
                datatype,
                recvData.data(),
                1,                      // one element per rank
                datatype,
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoall [comm: " << comm << "] failed."
                << " For " << sendData
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addAllToAllTime();
}


template<class Type>
void Foam::PstreamDetail::allToAllv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    MPI_Datatype datatype,
    const label comm,
    label* requestID
)
{
    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        if (requestID != nullptr)
        {
            Pout<< "** MPI_Ialltoallv (non-blocking):";
        }
        else
        {
            Pout<< "** MPI_Alltoallv (blocking):";
        }
        Pout<< " sendCounts:" << sendCounts
            << " sendOffsets:" << sendOffsets
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        (sendCounts.size() != np || sendOffsets.size() < np)
     || (recvCounts.size() != np || recvOffsets.size() < np)
    )
    {
        FatalErrorInFunction
            << "Have " << np << " ranks, but sendCounts:" << sendCounts.size()
            << ", sendOffsets:" << sendOffsets.size()
            << ", recvCounts:" << recvCounts.size()
            << " or recvOffsets:" << recvOffsets.size()
            << " is different!"
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        if (recvCounts[0] != sendCounts[0])
        {
            FatalErrorInFunction
                << "Bytes to send " << sendCounts[0]
                << " does not equal bytes to receive " << recvCounts[0]
                << Foam::abort(FatalError);
        }
        std::memmove
        (
            recvData,
            (sendData + sendOffsets[0]),
            recvCounts[0]*sizeof(Type)
        );
        return;
    }

    profilingPstream::beginTiming();

    bool handled(false);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (requestID != nullptr)
    {
        handled = true;
        MPI_Request request;

        if
        (
            MPI_Ialltoallv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Ialltoallv [comm: " << comm << "] failed."
                << " For sendCounts " << sendCounts
                << " recvCounts " << recvCounts
                << Foam::abort(FatalError);
        }

        if (PstreamGlobals::freedRequests_.size())
        {
            *requestID = PstreamGlobals::freedRequests_.remove();
            PstreamGlobals::outstandingRequests_[*requestID] = request;
        }
        else
        {
            *requestID = PstreamGlobals::outstandingRequests_.size();
            PstreamGlobals::outstandingRequests_.append(request);
        }
    }
#endif

    if (!handled)
    {
        if (requestID != nullptr)
        {
            *requestID = -1;
        }
        if
        (
            MPI_Alltoallv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoallv [comm: " << comm << "] failed."
                << " For sendCounts " << sendCounts
                << " recvCounts " << recvCounts
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addAllToAllTime();
}


template<class Type>
void Foam::PstreamDetail::gather
(
    const Type* sendData,
    int sendCount,

    Type* recvData,
    int recvCount,

    MPI_Datatype datatype,
    const label comm,
    label* requestID
)
{
    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvCount*sizeof(Type));
        return;
    }

    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        if (requestID != nullptr)
        {
            Pout<< "** MPI_Igather (non-blocking):";
        }
        else
        {
            Pout<< "** MPI_Gather (blocking):";
        }
        Pout<< " np:" << np
            << " recvCount:" << recvCount
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    profilingPstream::beginTiming();

    bool handled(false);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (requestID != nullptr)
    {
        handled = true;
        MPI_Request request;
        if
        (
            MPI_Igather
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Igather [comm: " << comm << "] failed."
                << " sendCount " << sendCount
                << " recvCount " << recvCount
                << Foam::abort(FatalError);
        }

        if (PstreamGlobals::freedRequests_.size())
        {
            *requestID = PstreamGlobals::freedRequests_.remove();
            PstreamGlobals::outstandingRequests_[*requestID] = request;
        }
        else
        {
            *requestID = PstreamGlobals::outstandingRequests_.size();
            PstreamGlobals::outstandingRequests_.append(request);
        }
    }
#endif

    if (!handled)
    {
        if (requestID != nullptr)
        {
            *requestID = -1;
        }
        if
        (
            MPI_Gather
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gather [comm: " << comm << "] failed."
                << " sendCount " << sendCount
                << " recvCount " << recvCount
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addGatherTime();
}


template<class Type>
void Foam::PstreamDetail::scatter
(
    const Type* sendData,
    int sendCount,

    Type* recvData,
    int recvCount,

    MPI_Datatype datatype,
    const label comm,
    label* requestID
)
{
    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvCount*sizeof(Type));
        return;
    }

    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        if (requestID != nullptr)
        {
            Pout<< "** MPI_Iscatter (non-blocking):";
        }
        else
        {
            Pout<< "** MPI_Scatter (blocking):";
        }
        Pout<< " np:" << np
            << " recvCount:" << recvCount
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    profilingPstream::beginTiming();

    bool handled(false);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (requestID != nullptr)
    {
        handled = true;
        MPI_Request request;
        if
        (
            MPI_Iscatter
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iscatter [comm: " << comm << "] failed."
                << " sendCount " << sendCount
                << " recvCount " << recvCount
                << Foam::abort(FatalError);
        }

        if (PstreamGlobals::freedRequests_.size())
        {
            *requestID = PstreamGlobals::freedRequests_.remove();
            PstreamGlobals::outstandingRequests_[*requestID] = request;
        }
        else
        {
            *requestID = PstreamGlobals::outstandingRequests_.size();
            PstreamGlobals::outstandingRequests_.append(request);
        }
    }
#endif

    if (!handled)
    {
        if (requestID != nullptr)
        {
            *requestID = -1;
        }
        if
        (
            MPI_Scatter
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iscatter [comm: " << comm << "] failed."
                << " sendCount " << sendCount
                << " recvCount " << recvCount
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addScatterTime();
}


template<class Type>
void Foam::PstreamDetail::gatherv
(
    const Type* sendData,
    int sendCount,

    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    MPI_Datatype datatype,
    const label comm,
    label* requestID
)
{
    if (!UPstream::parRun())
    {
        // recvCounts[0] may be invalid - use sendCount instead
        std::memmove(recvData, sendData, sendCount*sizeof(Type));
        return;
    }

    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        if (requestID != nullptr)
        {
            Pout<< "** MPI_Igatherv (non-blocking):";
        }
        else
        {
            Pout<< "** MPI_Gatherv (blocking):";
        }
        Pout<< " np:" << np
            << " recvCounts:" << recvCounts
            << " recvOffsets:" << recvOffsets
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        UPstream::master(comm)
     && (recvCounts.size() != np || recvOffsets.size() < np)
    )
    {
        // Note: allow offsets to be e.g. 1 larger than nProc so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Have " << np << " ranks, but recvCounts:" << recvCounts.size()
            << " or recvOffsets:" << recvOffsets.size()
            << " is too small!"
            << Foam::abort(FatalError);
    }

    profilingPstream::beginTiming();

    // Ensure send/recv consistency on master
    if (UPstream::master(comm) && !recvCounts[0])
    {
        sendCount = 0;
    }

    bool handled(false);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (requestID != nullptr)
    {
        handled = true;
        MPI_Request request;
        if
        (
            MPI_Igatherv
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Igatherv failed [comm: " << comm << ']'
                << " sendCount " << sendCount
                << " recvCounts " << recvCounts
                << Foam::abort(FatalError);
        }

        if (PstreamGlobals::freedRequests_.size())
        {
            *requestID = PstreamGlobals::freedRequests_.remove();
            PstreamGlobals::outstandingRequests_[*requestID] = request;
        }
        else
        {
            *requestID = PstreamGlobals::outstandingRequests_.size();
            PstreamGlobals::outstandingRequests_.append(request);
        }
    }
#endif

    if (!handled)
    {
        if (requestID != nullptr)
        {
            *requestID = -1;
        }
        if
        (
            MPI_Gatherv
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gatherv failed [comm: " << comm << ']'
                << " sendCount " << sendCount
                << " recvCounts " << recvCounts
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addGatherTime();
}


template<class Type>
void Foam::PstreamDetail::scatterv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    Type* recvData,
    int recvCount,

    MPI_Datatype datatype,
    const label comm,
    label* requestID
)
{
    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvCount*sizeof(Type));
        return;
    }

    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        if (requestID != nullptr)
        {
            Pout<< "** MPI_Iscatterv (non-blocking):";
        }
        else
        {
            Pout<< "** MPI_Scatterv (blocking):";
        }
        Pout<< " np:" << np
            << " sendCounts:" << sendCounts
            << " sendOffsets:" << sendOffsets
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        UPstream::master(comm)
     && (sendCounts.size() != np || sendOffsets.size() < np)
    )
    {
        // Note: allow offsets to be e.g. 1 larger than nProc so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Have " << np << " ranks, but sendCounts:" << sendCounts.size()
            << " or sendOffsets:" << sendOffsets.size()
            << " is too small!"
            << Foam::abort(FatalError);
    }

    profilingPstream::beginTiming();

    bool handled(false);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (requestID != nullptr)
    {
        handled = true;
        MPI_Request request;
        if
        (
            MPI_Iscatterv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iscatterv [comm: " << comm << "] failed."
                << " sendCounts " << sendCounts
                << " sendOffsets " << sendOffsets
                << Foam::abort(FatalError);
        }

        if (PstreamGlobals::freedRequests_.size())
        {
            *requestID = PstreamGlobals::freedRequests_.remove();
            PstreamGlobals::outstandingRequests_[*requestID] = request;
        }
        else
        {
            *requestID = PstreamGlobals::outstandingRequests_.size();
            PstreamGlobals::outstandingRequests_.append(request);
        }
    }
#endif

    if (!handled)
    {
        if (requestID != nullptr)
        {
            *requestID = -1;
        }
        if
        (
            MPI_Scatterv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatterv [comm: " << comm << "] failed."
                << " sendCounts " << sendCounts
                << " sendOffsets " << sendOffsets
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addScatterTime();
}


// ************************************************************************* //
