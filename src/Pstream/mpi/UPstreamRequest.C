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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Request::Request() noexcept
:
    UPstream::Request(MPI_REQUEST_NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Request::good() const noexcept
{
    return MPI_REQUEST_NULL != PstreamDetail::Request::get(*this);
}


void Foam::UPstream::Request::reset() noexcept
{
    *this = UPstream::Request(MPI_REQUEST_NULL);
}


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


void Foam::UPstream::waitRequests(const label pos)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || pos < 0
     || pos >= PstreamGlobals::outstandingRequests_.size()
     /// || !len
    )
    {
        return;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    /// // Treat len < 0 like npos (ie, the rest of the list) but also
    /// // apply range checking to avoid bad slices
    /// if (len > 0 && len < count)
    /// {
    ///     count = len;
    /// }

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << count << " requests starting at " << pos << endl;
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

    // ie, resetRequests(pos)
    PstreamGlobals::outstandingRequests_.resize(pos);

    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    for (auto& req : requests)
    {
        MPI_Request request = PstreamDetail::Request::get(req);

        if (MPI_REQUEST_NULL != request)
        {
            waitRequests[count] = request;
            ++count;
        }
    }

    if (!count)
    {
        // Early exit: non-NULL requests found
        return;
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

    // Everything handled, reset all to MPI_REQUEST_NULL
    requests = UPstream::Request(MPI_REQUEST_NULL);
}


Foam::label Foam::UPstream::waitAnyRequest(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return -1;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    // Transcribe UPstream::Request into MPI_Request
    // - do not change locations within the list since these are relevant
    //   for the return index.
    for (auto& req : requests)
    {
        waitRequests[count] = PstreamDetail::Request::get(req);
        ++count;
    }

    profilingPstream::beginTiming();

    // On success: sets request to MPI_REQUEST_NULL
    int index = -1;
    if (MPI_Waitany(count, waitRequests, &index, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Waitany returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (index == MPI_UNDEFINED)
    {
        index = -1;  // No outstanding requests
    }

    // Transcribe MPI_Request back into UPstream::Request
    while (--count >= 0)
    {
        requests[count] = UPstream::Request(waitRequests[count]);
    }

    return index;
}


// FUTURE?
//
/// void Foam::UPstream::waitRequests
/// (
///     UPstream::Request& req1,
///     UPstream::Request& req2
/// )
/// {
///     // No-op for non-parallel
///     if (!UPstream::parRun())
///     {
///         return;
///     }
///
///     int count = 0;
///     MPI_Request waitRequests[2];
///
///     waitRequests[count] = PstreamDetail::Request::get(req1);
///     if (MPI_REQUEST_NULL != waitRequests[count])
///     {
///         // Flag in advance as being handled
///         req1 = UPstream::Request(MPI_REQUEST_NULL);
///         ++count;
///     }
///
///     waitRequests[count] = PstreamDetail::Request::get(req2);
///     if (MPI_REQUEST_NULL != waitRequests[count])
///     {
///         // Flag in advance as being handled
///         req2 = UPstream::Request(MPI_REQUEST_NULL);
///         ++count;
///     }
///
///     if (!count)
///     {
///         return;
///     }
///
///     profilingPstream::beginTiming();
///
///     // On success: sets each request to MPI_REQUEST_NULL
///     if (MPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
///     {
///         FatalErrorInFunction
///             << "MPI_Waitall returned with error"
///             << Foam::abort(FatalError);
///     }
///
///     profilingPstream::addWaitTime();
/// }


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


void Foam::UPstream::waitRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    MPI_Request request = PstreamDetail::Request::get(req);

    // No-op for null request
    if (MPI_REQUEST_NULL == request)
    {
        return;
    }

    profilingPstream::beginTiming();

    if (MPI_Wait(&request, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    // Handled, reset to MPI_REQUEST_NULL
    req = UPstream::Request(MPI_REQUEST_NULL);
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


bool Foam::UPstream::finishedRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return true;
    }

    MPI_Request request = PstreamDetail::Request::get(req);

    // No-op for null request
    if (MPI_REQUEST_NULL == request)
    {
        return true;
    }

    int flag = 0;
    MPI_Test(&request, &flag, MPI_STATUS_IGNORE);

    if (flag)
    {
        // Success: reset request to MPI_REQUEST_NULL
        req = UPstream::Request(MPI_REQUEST_NULL);
    }

    return flag != 0;
}


bool Foam::UPstream::finishedRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return true;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    for (auto& req : requests)
    {
        MPI_Request request = PstreamDetail::Request::get(req);

        if (MPI_REQUEST_NULL != request)
        {
            waitRequests[count] = request;
            ++count;
        }
    }

    if (!count)
    {
        // Early exit: non-NULL requests found
        return true;
    }

    // On success: sets each request to MPI_REQUEST_NULL
    // On failure: no request is modified
    int flag = 0;
    MPI_Testall(count, waitRequests, &flag, MPI_STATUSES_IGNORE);

    if (flag)
    {
        // Success: reset all requests to MPI_REQUEST_NULL
        requests = UPstream::Request(MPI_REQUEST_NULL);
    }
    else
    {
        // Not all done. Recover wrapped representation but in reverse order
        // since sizeof(MPI_Request) can be smaller than
        // sizeof(UPstream::Request::value_type)
        // eg, mpich has MPI_Request as 'int'
        //
        // This is uglier that we'd like, but much better than allocating
        // and freeing a scratch buffer each time we query things.

        while (--count >= 0)
        {
            requests[count] = UPstream::Request(waitRequests[count]);
        }
    }

    return flag != 0;
}


// ************************************************************************* //
