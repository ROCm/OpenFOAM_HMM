/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
        // all-to-all
        Pstream::exchangeSizes(sendBuf_, recvSizes, comm_);

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvSizes,
            recvBuf_,
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
            sendBuf_,
            recvSizes,
            tag_,
            comm_
        );

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvSizes,
            recvBuf_,
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

            recvSizes = UPstream::listGatherValues(sendBuf_[0].size(), comm_);

            if (!UPstream::master(comm_))
            {
                recvSizes.resize_nocopy(recvBuf_.size());
                recvSizes = Zero;
            }
        }
        else
        {
            // scatter mode (one-to-all): master [0] -> everyone

            recvSizes.resize_nocopy(sendBuf_.size());

            if (UPstream::master(comm_))
            {
                forAll(sendBuf_, proci)
                {
                    recvSizes[proci] = sendBuf_[proci].size();
                }
            }

            const label myRecv(UPstream::listScatterValues(recvSizes, comm_));

            recvSizes = Zero;
            recvSizes[0] = myRecv;
        }


        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvSizes,
            recvBuf_,
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
    sendBuf_(UPstream::nProcs(comm_)),
    recvBuf_(UPstream::nProcs(comm_)),
    recvBufPos_(UPstream::nProcs(comm_), Zero)
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
    sendBuf_(UPstream::nProcs(comm_)),
    recvBuf_(UPstream::nProcs(comm_)),
    recvBufPos_(UPstream::nProcs(comm_), Zero)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::~PstreamBuffers()
{
    // Check that all data has been consumed.
    forAll(recvBufPos_, proci)
    {
        if (recvBufPos_[proci] < recvBuf_[proci].size())
        {
            FatalErrorInFunction
                << "Message from processor " << proci
                << " Only consumed " << recvBufPos_[proci] << " of "
                << recvBuf_[proci].size() << " bytes" << nl
                << Foam::abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PstreamBuffers::clear()
{
    for (DynamicList<char>& buf : sendBuf_)
    {
        buf.clear();
    }
    for (DynamicList<char>& buf : recvBuf_)
    {
        buf.clear();
    }
    recvBufPos_ = 0;

    finishedSendsCalled_ = false;
}


void Foam::PstreamBuffers::clearRecv(const label proci)
{
    recvBuf_[proci].clear();
    recvBufPos_[proci] = 0;
}


void Foam::PstreamBuffers::clearStorage()
{
    // Could also clear out entire sendBuf_, recvBuf_ and reallocate.
    // Not sure if it makes much difference
    for (DynamicList<char>& buf : sendBuf_)
    {
        buf.clearStorage();
    }
    for (DynamicList<char>& buf : recvBuf_)
    {
        buf.clearStorage();
    }
    recvBufPos_ = 0;

    finishedSendsCalled_ = false;
}


bool Foam::PstreamBuffers::hasSendData() const
{
    for (const DynamicList<char>& buf : sendBuf_)
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
        forAll(recvBufPos_, proci)
        {
            if (recvBuf_[proci].size() > recvBufPos_[proci])
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
    return sendBuf_[proci].size();
}


Foam::label Foam::PstreamBuffers::recvDataCount(const label proci) const
{
    if (finishedSendsCalled_)
    {
        const label len(recvBuf_[proci].size() > recvBufPos_[proci]);

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
    labelList counts(recvBuf_.size(), Zero);

    if (finishedSendsCalled_)
    {
        forAll(recvBufPos_, proci)
        {
            const label len(recvBuf_[proci].size() - recvBufPos_[proci]);

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
        const label len(recvBuf_[proci].size() - recvBufPos_[proci]);

        if (len > 0)
        {
            return UList<char>
            (
                const_cast<char*>(&recvBuf_[proci][recvBufPos_[proci]]),
                len
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
    forAll(sendBuf_, proci)
    {
        // ie, sendDataCount(proci) != 0
        if (sendConnections.set(proci, !sendBuf_[proci].empty()))
        {
            // The state changed
            changed = true;
        }
    }

    reduce(changed, orOp<bool>());

    if (changed)
    {
        // Create send/recv topology

        // The send ranks
        sendProcs.clear();
        forAll(sendBuf_, proci)
        {
            // ie, sendDataCount(proci) != 0
            if (!sendBuf_[proci].empty())
            {
                sendProcs.append(proci);
            }
        }

        finishedSends(wait);  // All-to-all

        // The recv ranks
        recvProcs.clear();
        forAll(recvBuf_, proci)
        {
            // ie, recvDataCount(proci)
            if (!recvBuf_[proci].empty())
            {
                recvProcs.append(proci);
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
