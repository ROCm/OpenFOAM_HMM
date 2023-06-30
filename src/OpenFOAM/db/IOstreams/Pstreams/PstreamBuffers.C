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
#include "debug.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::PstreamBuffers::algorithm
(
    // Name may change in the future (JUN-2023)
    Foam::debug::optimisationSwitch("pbufs.tuning", -1)
);
registerOptSwitch
(
    "pbufs.tuning",
    int,
    Foam::PstreamBuffers::algorithm
);


// Simple enumerations
// -------------------
static constexpr int algorithm_PEX_allToAll = -1;  // Traditional PEX
//static constexpr int algorithm_PEX_hybrid = 0;   // Possible new default?
static constexpr int algorithm_full_NBX = 1;       // Very experimental


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PstreamBuffers::finalExchange
(
    const bool wait,
    labelList& recvSizes
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
         && (algorithm >= algorithm_full_NBX)
         && (UPstream::maxCommsSize <= 0)
        )
        {
            // NBX algorithm (nonblocking exchange)
            // - when requested and waiting, no data chunking etc

            PstreamDetail::exchangeConsensus<DynamicList<char>, char>
            (
                sendBuffers_,
                recvBuffers_,
                recvSizes,
                (tag_ + 271828),  // some unique tag?
                comm_,
                wait
            );

            return;
        }


        // PEX algorithm with two different flavours of exchanging sizes

        // Assemble the send sizes (cf. Pstream::exchangeSizes)
        labelList sendSizes(nProcs_);
        forAll(sendBuffers_, proci)
        {
            sendSizes[proci] = sendBuffers_[proci].size();
        }
        recvSizes.resize_nocopy(nProcs_);

        if (algorithm == algorithm_PEX_allToAll)
        {
            // PEX stage 1: exchange sizes (all-to-all)
            UPstream::allToAll(sendSizes, recvSizes, comm_);
        }
        else
        {
            // PEX stage 1: exchange sizes (non-blocking consensus)
            UPstream::allToAllConsensus
            (
                sendSizes,
                recvSizes,
                (tag_ + 314159),  // some unique tag?
                comm_
            );
        }

        // PEX stage 2: point-to-point data exchange
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
    const bool wait,
    labelList& recvSizes
)
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        // Preparation. Temporarily abuse recvSizes as logic to clear
        // send buffers that are not in the neighbourhood connection
        {
            recvSizes.resize_nocopy(nProcs_);
            recvSizes = 0;

            // Preserve self-send, even if not described by neighbourhood
            recvSizes[UPstream::myProcNo(comm_)] = 1;

            for (const label proci : sendProcs)
            {
                recvSizes[proci] = 1;  // Connected
            }

            for (label proci=0; proci < nProcs_; ++proci)
            {
                if (!recvSizes[proci])  // Not connected
                {
                    sendBuffers_[proci].clear();
                }
            }
        }

        // PEX stage 1: exchange sizes (limited neighbourhood)
        Pstream::exchangeSizes
        (
            sendProcs,
            recvProcs,
            sendBuffers_,
            recvSizes,
            tag_,
            comm_
        );

        // PEX stage 2: point-to-point data exchange
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


void Foam::PstreamBuffers::finalGatherScatter
(
    const bool isGather,
    const bool wait,
    labelList& recvSizes
)
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    finishedSendsCalled_ = true;

    if (isGather)
    {
        // gather mode (all-to-one)

        // Only send to master [0]. Master is also allowed to 'send' to itself

        for (label proci=1; proci < sendBuffers_.size(); ++proci)
        {
            sendBuffers_[proci].clear();
        }
    }
    else
    {
        // scatter mode (one-to-all)

        if (!UPstream::master(comm_))
        {
            // Non-master: has no sends
            clearSends();
        }
    }


    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        // Use PEX algorithm
        // - for a non-sparse gather/scatter, it is presumed that
        // MPI_Gather/MPI_Scatter will be the most efficient way to
        // communicate the sizes.

        // PEX stage 1: exchange sizes (gather or scatter)
        if (isGather)
        {
            // gather mode (all-to-one): master [0] <- everyone

            recvSizes =
                UPstream::listGatherValues(sendBuffers_[0].size(), comm_);

            if (!UPstream::master(comm_))
            {
                recvSizes.resize_nocopy(nProcs_);
                recvSizes = Zero;
            }
        }
        else
        {
            // scatter mode (one-to-all): master [0] -> everyone

            recvSizes.resize_nocopy(nProcs_);

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

        // PEX stage 2: point-to-point data exchange
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
    UPstream::commsTypes commsType,
    int tag,
    label communicator,
    IOstreamOption::streamFormat fmt
)
:
    finishedSendsCalled_(false),
    allowClearRecv_(true),
    format_(fmt),
    commsType_(commsType),
    tag_(tag),
    comm_(communicator),
    nProcs_(UPstream::nProcs(comm_)),
    sendBuffers_(nProcs_),
    recvBuffers_(nProcs_),
    recvPositions_(nProcs_, Zero)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::~PstreamBuffers()
{
    // Check that all data has been consumed.
    forAll(recvBuffers_, proci)
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

void Foam::PstreamBuffers::clearSends()
{
    for (DynamicList<char>& buf : sendBuffers_)
    {
        buf.clear();
    }
}


void Foam::PstreamBuffers::clearRecvs()
{
    for (DynamicList<char>& buf : recvBuffers_)
    {
        buf.clear();
    }
    recvPositions_ = Zero;
}


void Foam::PstreamBuffers::clear()
{
    clearSends();
    clearRecvs();
    finishedSendsCalled_ = false;
}


void Foam::PstreamBuffers::clearSend(const label proci)
{
    sendBuffers_[proci].clear();
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
    recvPositions_ = Zero;

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
        forAll(recvBuffers_, proci)
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
    labelList counts(nProcs_, Zero);

    if (finishedSendsCalled_)
    {
        forAll(recvBuffers_, proci)
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


Foam::label Foam::PstreamBuffers::maxNonLocalRecvCount
(
    const label excludeProci
) const
{
    label maxLen = 0;

    if (finishedSendsCalled_)
    {
        forAll(recvBuffers_, proci)
        {
            if (excludeProci != proci)
            {
                label len(recvBuffers_[proci].size() - recvPositions_[proci]);
                maxLen = max(maxLen, len);
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

    return maxLen;
}


Foam::label Foam::PstreamBuffers::maxRecvCount() const
{
    // Use out-of-range proci to avoid excluding any processor
    return maxNonLocalRecvCount(-1);
}


Foam::label Foam::PstreamBuffers::maxNonLocalRecvCount() const
{
    return maxNonLocalRecvCount(UPstream::myProcNo(comm_));
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
                const_cast<char*>(recvBuffers_[proci].cbegin(pos)),
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
    finalExchange(wait, recvSizes);
}


void Foam::PstreamBuffers::finishedSends
(
    labelList& recvSizes,
    const bool wait
)
{
    // Resize for copying back
    recvSizes.resize_nocopy(sendBuffers_.size());

    finalExchange(wait, recvSizes);

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


void Foam::PstreamBuffers::finishedNeighbourSends
(
    const labelUList& neighProcs,
    labelList& recvSizes,
    const bool wait
)
{
    finalExchange(neighProcs, neighProcs, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedNeighbourSends
(
    const labelUList& neighProcs,
    const bool wait
)
{
    labelList recvSizes;
    finalExchange(neighProcs, neighProcs, wait, recvSizes);
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
    forAll(sendBuffers_, proci)
    {
        if (sendConnections.set(proci, !sendBuffers_[proci].empty()))
        {
            // The state changed
            changed = true;
        }
    }

    UPstream::reduceOr(changed, comm_);

    if (changed)
    {
        // Update send/recv topology
        labelList recvSizes;
        finishedSends(recvSizes, wait);  // eg, using all-to-all

        // The send ranks
        sendProcs.clear();
        forAll(sendBuffers_, proci)
        {
            if (!sendBuffers_[proci].empty())
            {
                sendProcs.push_back(proci);
            }
        }

        // The recv ranks
        recvProcs.clear();
        forAll(recvSizes, proci)
        {
            if (recvSizes[proci] > 0)
            {
                recvProcs.push_back(proci);
            }
        }
    }
    else
    {
        // Use existing send/recv ranks
        labelList recvSizes;
        finalExchange(sendProcs, recvProcs, wait, recvSizes);
    }

    return changed;
}


void Foam::PstreamBuffers::finishedGathers(const bool wait)
{
    labelList recvSizes;
    finalGatherScatter(true, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedScatters(const bool wait)
{
    labelList recvSizes;
    finalGatherScatter(false, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedGathers
(
    labelList& recvSizes,
    const bool wait
)
{
    finalGatherScatter(true, wait, recvSizes);

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


void Foam::PstreamBuffers::finishedScatters
(
    labelList& recvSizes,
    const bool wait
)
{
    finalGatherScatter(false, wait, recvSizes);

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


// ************************************************************************* //
