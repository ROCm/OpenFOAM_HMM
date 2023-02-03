/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "PstreamBuffers.H"
#include "bitSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PstreamBuffers::finalExchange
(
    labelList& recvSizes,
    const bool wait
)
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        if
        (
            wait
         && UPstream::parRun()
         && UPstream::nProcsNonblockingExchange > 1
         && UPstream::nProcsNonblockingExchange <= nProcs()
        )
        {
            Pstream::exchangeConsensus<DynamicList<char>, char>
            (
                sendBuffers_,
                recvBuffers_,
                (tag_ + 314159),  // some unique tag?
                comm_
            );

            // Copy back out
            recvSizes.resize_nocopy(recvBuffers_.size());
            forAll(recvBuffers_, proci)
            {
                recvSizes[proci] = recvBuffers_[proci].size();
            }

            return;
        }

        // all-to-all
        Pstream::exchangeSizes(sendBuffers_, recvSizes, comm_);

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuffers_,
            recvSizes,
            recvBuffers_,
            tag_,
            comm_,
            wait
        );
    }
}


void Foam::PstreamBuffers::finalExchange
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    labelList& recvSizes,
    const bool wait
)
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        Pstream::exchangeSizes
        (
            sendProcs,
            recvProcs,
            sendBuffers_,
            recvSizes,
            tag_,
            comm_
        );

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuffers_,
            recvSizes,
            recvBuffers_,
            tag_,
            comm_,
            wait
        );
    }
}


void Foam::PstreamBuffers::finalExchangeGatherScatter
(
    const bool isGather,
    const bool wait
)
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        labelList recvSizes;

        if (isGather)
        {
            // gather mode (all-to-one): master [0] <- everyone

            recvSizes =
                UPstream::listGatherValues(sendBuffers_[0].size(), comm_);

            if (!UPstream::master(comm_))
            {
                recvSizes.resize_nocopy(recvBuffers_.size());
                recvSizes = Zero;
            }
        }
        else
        {
            // scatter mode (one-to-all): master [0] -> everyone

            recvSizes.resize_nocopy(sendBuffers_.size());

            if (UPstream::master(comm_))
            {
                forAll(sendBuffers_, proci)
                {
                    recvSizes[proci] = sendBuffers_[proci].size();
                }
            }

            const label myRecv(UPstream::listScatterValues(recvSizes, comm_));

            recvSizes = Zero;
            recvSizes[0] = myRecv;
        }


        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuffers_,
            recvSizes,
            recvBuffers_,
            tag_,
            comm_,
            wait
        );
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::PstreamBuffers
(
    const UPstream::commsTypes commsType,
    const int tag,
    const label comm,
    IOstreamOption::streamFormat fmt
)
:
    finishedSendsCalled_(false),
    allowClearRecv_(true),
    format_(fmt),
    commsType_(commsType),
    tag_(tag),
    comm_(comm),
    sendBuffers_(UPstream::nProcs(comm)),
    recvBuffers_(UPstream::nProcs(comm)),
    recvPositions_(UPstream::nProcs(comm), Zero)
{}


Foam::PstreamBuffers::PstreamBuffers
(
    const label comm,
    const UPstream::commsTypes commsType,
    const int tag,
    IOstreamOption::streamFormat fmt
)
:
    finishedSendsCalled_(false),
    allowClearRecv_(true),
    format_(fmt),
    commsType_(commsType),
    tag_(tag),
    comm_(comm),
    sendBuffers_(UPstream::nProcs(comm)),
    recvBuffers_(UPstream::nProcs(comm)),
    recvPositions_(UPstream::nProcs(comm), Zero)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::~PstreamBuffers()
{
    // Check that all data has been consumed.
    forAll(recvPositions_, proci)
    {
        const label pos = recvPositions_[proci];
        const label len = recvBuffers_[proci].size();

        if (pos < len)
        {
            FatalErrorInFunction
                << "Message from processor " << proci
                << " Only consumed " << pos << " of " << len << " bytes" << nl
                << Foam::abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::DynamicList<char>& Foam::PstreamBuffers::accessSendBuffer
(
    const label proci
)
{
    return sendBuffers_[proci];
}


Foam::DynamicList<char>& Foam::PstreamBuffers::accessRecvBuffer
(
    const label proci
)
{
    return recvBuffers_[proci];
}


Foam::label& Foam::PstreamBuffers::accessRecvPosition(const label proci)
{
    return recvPositions_[proci];
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PstreamBuffers::clear()
{
    for (DynamicList<char>& buf : sendBuffers_)
    {
        buf.clear();
    }
    for (DynamicList<char>& buf : recvBuffers_)
    {
        buf.clear();
    }
    recvPositions_ = 0;

    finishedSendsCalled_ = false;
}


void Foam::PstreamBuffers::clearRecv(const label proci)
{
    recvBuffers_[proci].clear();
    recvPositions_[proci] = 0;
}


void Foam::PstreamBuffers::clearStorage()
{
    // Could also clear out entire sendBuffers_, recvBuffers_ and reallocate.
    // Not sure if it makes much difference
    for (DynamicList<char>& buf : sendBuffers_)
    {
        buf.clearStorage();
    }
    for (DynamicList<char>& buf : recvBuffers_)
    {
        buf.clearStorage();
    }
    recvPositions_ = 0;

    finishedSendsCalled_ = false;
}


bool Foam::PstreamBuffers::hasSendData() const
{
    for (const DynamicList<char>& buf : sendBuffers_)
    {
        if (!buf.empty())
        {
            return true;
        }
    }
    return false;
}


bool Foam::PstreamBuffers::hasRecvData() const
{
    if (finishedSendsCalled_)
    {
        forAll(recvPositions_, proci)
        {
            if (recvPositions_[proci] < recvBuffers_[proci].size())
            {
                return true;
            }
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return false;
}


Foam::label Foam::PstreamBuffers::sendDataCount(const label proci) const
{
    return sendBuffers_[proci].size();
}


Foam::label Foam::PstreamBuffers::recvDataCount(const label proci) const
{
    if (finishedSendsCalled_)
    {
        const label len(recvBuffers_[proci].size() - recvPositions_[proci]);

        if (len > 0)
        {
            return len;
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return 0;
}


Foam::labelList Foam::PstreamBuffers::recvDataCounts() const
{
    labelList counts(recvPositions_.size(), Zero);

    if (finishedSendsCalled_)
    {
        forAll(recvPositions_, proci)
        {
            const label len(recvBuffers_[proci].size() - recvPositions_[proci]);

            if (len > 0)
            {
                counts[proci] = len;
            }
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return counts;
}


const Foam::UList<char>
Foam::PstreamBuffers::peekRecvData(const label proci) const
{
    if (finishedSendsCalled_)
    {
        const label pos = recvPositions_[proci];
        const label len = recvBuffers_[proci].size();

        if (pos < len)
        {
            return UList<char>
            (
                const_cast<char*>(recvBuffers_[proci].cdata()) + pos,
                (len - pos)
            );
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return UList<char>();
}


bool Foam::PstreamBuffers::allowClearRecv(bool on) noexcept
{
    bool old(allowClearRecv_);
    allowClearRecv_ = on;
    return old;
}


void Foam::PstreamBuffers::finishedSends(const bool wait)
{
    labelList recvSizes;
    finalExchange(recvSizes, wait);
}


void Foam::PstreamBuffers::finishedSends
(
    labelList& recvSizes,
    const bool wait
)
{
    finalExchange(recvSizes, wait);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


void Foam::PstreamBuffers::finishedSends
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const bool wait
)
{
    labelList recvSizes;
    finalExchange(sendProcs, recvProcs, recvSizes, wait);
}


void Foam::PstreamBuffers::finishedSends
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    labelList& recvSizes,
    const bool wait
)
{
    finalExchange(sendProcs, recvProcs, recvSizes, wait);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


bool Foam::PstreamBuffers::finishedSends
(
    bitSet& sendConnections,
    DynamicList<label>& sendProcs,
    DynamicList<label>& recvProcs,
    const bool wait
)
{
    bool changed = (sendConnections.size() != nProcs());

    if (changed)
    {
        sendConnections.resize(nProcs());
    }

    // Update send connections
    // - reasonable to assume there are no self-sends on UPstream::myProcNo
    forAll(sendBuffers_, proci)
    {
        // ie, sendDataCount(proci) != 0
        if (sendConnections.set(proci, !sendBuffers_[proci].empty()))
        {
            // The state changed
            changed = true;
        }
    }

    UPstream::reduceOr(changed, comm_);

    if (changed)
    {
        // Create send/recv topology

        // The send ranks
        sendProcs.clear();
        forAll(sendBuffers_, proci)
        {
            // ie, sendDataCount(proci) != 0
            if (!sendBuffers_[proci].empty())
            {
                sendProcs.push_back(proci);
            }
        }

        finishedSends(wait);  // All-to-all

        // The recv ranks
        recvProcs.clear();
        forAll(recvBuffers_, proci)
        {
            // ie, recvDataCount(proci)
            if (!recvBuffers_[proci].empty())
            {
                recvProcs.push_back(proci);
            }
        }
    }
    else
    {
        // Use existing send/recv ranks

        finishedSends(sendProcs, recvProcs, wait);
    }

    return changed;
}


void Foam::PstreamBuffers::finishedGathers(const bool wait)
{
    finalExchangeGatherScatter(true, wait);
}


void Foam::PstreamBuffers::finishedScatters(const bool wait)
{
    finalExchangeGatherScatter(false, wait);
}


void Foam::PstreamBuffers::finishedGathers
(
    labelList& recvSizes,
    const bool wait
)
{
    finalExchangeGatherScatter(true, wait);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }

    // For nonBlocking mode, simply recover received sizes
    // from the buffers themselves.

    recvSizes = recvDataCounts();
}


void Foam::PstreamBuffers::finishedScatters
(
    labelList& recvSizes,
    const bool wait
)
{
    finalExchangeGatherScatter(false, wait);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }

    // For nonBlocking mode, simply recover received sizes
    // from the buffers themselves.

    recvSizes = recvDataCounts();
}


// ************************************************************************* //
