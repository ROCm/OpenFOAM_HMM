/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
#include "PstreamGlobals.H"
#include "profilingPstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UPstream::nRequests() noexcept
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label n)
{
    if (n >= 0 && n < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.resize(n);
    }
}


void Foam::UPstream::waitRequests(const label start)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || start < 0
     || start >= PstreamGlobals::outstandingRequests_.size()
    )
    {
        return;
    }

    const label count = (PstreamGlobals::outstandingRequests_.size() - start);
    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + start);

    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << count << " requests starting at " << start << endl;
    }

    profilingPstream::beginTiming();

    // On success: sets each request to MPI_REQUEST_NULL
    if (MPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Waitall returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    // ie, resetRequests(start)
    PstreamGlobals::outstandingRequests_.resize(start);

    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequest(const label i)
{
    // No-op for non-parallel, or out-of-range (eg, placeholder indices)
    if
    (
        !UPstream::parRun()
     || i < 0
     || i >= PstreamGlobals::outstandingRequests_.size()
    )
    {
        return;
    }

    // Push index onto free cache (for later reuse)
    PstreamGlobals::freedRequests_.push_back(i);

    auto& request = PstreamGlobals::outstandingRequests_[i];

    // No-op for null request
    if (MPI_REQUEST_NULL == request)
    {
        return;
    }

    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequest : starting wait for request:"
            << i << endl;
    }

    profilingPstream::beginTiming();

    // On success: sets request to MPI_REQUEST_NULL
    if (MPI_Wait(&request, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequest : finished wait for request:"
            << i << endl;
    }
}


bool Foam::UPstream::finishedRequest(const label i)
{
    // No-op for non-parallel, or out-of-range (eg, placeholder indices)
    if
    (
        !UPstream::parRun()
     || i < 0
     || i >= PstreamGlobals::outstandingRequests_.size()
    )
    {
        return true;
    }

    auto& request = PstreamGlobals::outstandingRequests_[i];

    // No-op for null request
    if (MPI_REQUEST_NULL == request)
    {
        return true;
    }

    if (UPstream::debug)
    {
        Pout<< "UPstream::finishedRequest : checking request:"
            << i << endl;
    }

    // On success: sets request to MPI_REQUEST_NULL
    int flag = 0;
    MPI_Test(&request, &flag, MPI_STATUS_IGNORE);

    if (UPstream::debug)
    {
        Pout<< "UPstream::finishedRequest : finished request:" << i
            << endl;
    }

    return flag != 0;
}


// ************************************************************************* //
