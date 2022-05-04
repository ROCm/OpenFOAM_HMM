/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "mapDistributeBase.H"
#include "bitSet.H"
#include "commSchedule.H"
#include "labelPairHashes.H"
#include "globalIndex.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mapDistributeBase, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::mapDistributeBase::hasFlipAddressing(const labelUList& map)
{
    for (const label val : map)
    {
        if (!val)
        {
            // Cannot be flipped addressing if it contains zero.
            return false;
        }
        else if (val < 0)
        {
            // Must be flipped addressing if it contains negatives.
            return true;
        }
    }

    return false;
}


bool Foam::mapDistributeBase::hasFlipAddressing(const labelListList& maps)
{
    for (const labelList& map : maps)
    {
        for (const label val : map)
        {
            if (!val)
            {
                // Cannot be flipped addressing if it contains zero.
                return false;
            }
            else if (val < 0)
            {
                // Must be flipped addressing if it contains negatives.
                return true;
            }
        }
    }

    return false;
}


Foam::label Foam::mapDistributeBase::getMappedSize
(
    const labelListList& maps,
    const bool hasFlip
)
{
    label maxIndex = -1;

    for (const labelList& map : maps)
    {
        for (label index : map)
        {
            if (hasFlip)
            {
                index = mag(index)-1;
            }

            maxIndex = max(maxIndex, index);
        }
    }

    return (maxIndex+1);
}


Foam::label Foam::mapDistributeBase::countUnmapped
(
    const labelUList& elements,
    const labelListList& maps,
    const bool hasFlip
)
{
    if (elements.empty())
    {
        return 0;
    }

    // Moderately efficient markup/search

    bitSet unvisited(elements);
    label nUnmapped = unvisited.count();

    if (hasFlip)
    {
        for (const labelList& map : maps)
        {
            for (label index : map)
            {
                index = mag(index)-1;

                if (unvisited.unset(index))
                {
                    --nUnmapped;
                    if (!nUnmapped) break;
                }
            }
        }
    }
    else
    {
        for (const labelList& map : maps)
        {
            for (label index : map)
            {
                if (unvisited.unset(index))
                {
                    --nUnmapped;
                    if (!nUnmapped) break;
                }
            }
        }
    }

    return nUnmapped;
}


void Foam::mapDistributeBase::checkReceivedSize
(
    const label proci,
    const label expectedSize,
    const label receivedSize
)
{
    if (receivedSize != expectedSize)
    {
        FatalErrorInFunction
            << "Expected from processor " << proci
            << " " << expectedSize << " but received "
            << receivedSize << " elements."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::List<Foam::labelPair> Foam::mapDistributeBase::schedule
(
    const labelListList& subMap,
    const labelListList& constructMap,
    const int tag,
    const label comm
)
{
    const label myRank = Pstream::myProcNo(comm);
    const label nProcs = Pstream::nProcs(comm);

    // Communications: send and receive processor
    List<labelPair> allComms;

    {
        labelPairHashSet commsSet(nProcs);

        // Find what communication is required
        forAll(subMap, proci)
        {
            if (proci != myRank)
            {
                if (subMap[proci].size())
                {
                    // I need to send to proci
                    commsSet.insert(labelPair(myRank, proci));
                }
                if (constructMap[proci].size())
                {
                    // I need to receive from proci
                    commsSet.insert(labelPair(proci, myRank));
                }
            }
        }
        allComms = commsSet.toc();
    }


    // Gather/reduce
    if (Pstream::master(comm))
    {
        // Receive and merge
        for (const int proci : Pstream::subProcs(comm))
        {
            IPstream fromProc
            (
                Pstream::commsTypes::scheduled,
                proci,
                0,
                tag,
                comm
            );
            List<labelPair> nbrData(fromProc);

            for (const labelPair& connection : nbrData)
            {
                allComms.appendUniq(connection);
            }
        }
    }
    else
    {
        if (Pstream::parRun())
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                0,
                tag,
                comm
            );
            toMaster << allComms;
        }
    }

    // Broadcast: send comms information to all
    Pstream::broadcast(allComms, comm);

    // Determine my schedule.
    labelList mySchedule
    (
        commSchedule
        (
            nProcs,
            allComms
        ).procSchedule()[myRank]
    );

    // Processors involved in my schedule
    return List<labelPair>(allComms, mySchedule);
}


const Foam::List<Foam::labelPair>& Foam::mapDistributeBase::schedule() const
{
    if (!schedulePtr_)
    {
        schedulePtr_.reset
        (
            new List<labelPair>
            (
                schedule(subMap_, constructMap_, UPstream::msgType(), comm_)
            )
        );
    }

    return *schedulePtr_;
}


const Foam::List<Foam::labelPair>& Foam::mapDistributeBase::whichSchedule
(
    const UPstream::commsTypes commsType
) const
{
    if (commsType == UPstream::commsTypes::scheduled)
    {
        return schedule();
    }

    return List<labelPair>::null();
}


void Foam::mapDistributeBase::printLayout(Ostream& os) const
{
    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    // Determine offsets of remote data.
    labelList minIndex(nProcs, labelMax);
    labelList maxIndex(nProcs, labelMin);
    forAll(constructMap_, proci)
    {
        const labelList& construct = constructMap_[proci];
        if (constructHasFlip_)
        {
            forAll(construct, i)
            {
                label index = mag(construct[i])-1;
                minIndex[proci] = min(minIndex[proci], index);
                maxIndex[proci] = max(maxIndex[proci], index);
            }
        }
        else
        {
            forAll(construct, i)
            {
                label index = construct[i];
                minIndex[proci] = min(minIndex[proci], index);
                maxIndex[proci] = max(maxIndex[proci], index);
            }
        }
    }

    label localSize(0);

    if (maxIndex[myRank] != labelMin)
    {
        localSize = maxIndex[myRank]+1;
    }

    os  << "Layout: (constructSize:" << constructSize_
        << " subHasFlip:" << subHasFlip_
        << " constructHasFlip:" << constructHasFlip_
        << ")" << nl
        << "local (processor " << myRank << "):" << nl
        << "    start : 0" << nl
        << "    size  : " << localSize << endl;

    label offset = localSize;
    forAll(minIndex, proci)
    {
        if (proci != myRank && !constructMap_[proci].empty())
        {
            label size(0);

            if (maxIndex[proci] != labelMin)
            {
                size = maxIndex[proci]-minIndex[proci]+1;
                if (minIndex[proci] != offset)
                {
                    FatalErrorInFunction
                        << "offset:" << offset
                        << " proci:" << proci
                        << " minIndex:" << minIndex[proci]
                        << abort(FatalError);
                }
            }

            os  << "processor " << proci << ':' << nl
                << "    start : " << offset << nl
                << "    size  : " << size << endl;

            offset += size;
        }
    }
}


void Foam::mapDistributeBase::calcCompactAddressing
(
    const globalIndex& globalNumbering,
    const labelUList& elements,
    List<Map<label>>& compactMap
) const
{
    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    // Count all (non-local) elements needed. Just for presizing map.
    labelList nNonLocal(nProcs, Zero);

    for (const label globalIdx : elements)
    {
        if (globalIdx != -1 && !globalNumbering.isLocal(globalIdx))
        {
            label proci = globalNumbering.whichProcID(globalIdx);
            nNonLocal[proci]++;
        }
    }

    compactMap.resize_nocopy(nProcs);

    forAll(compactMap, proci)
    {
        compactMap[proci].clear();
        if (proci != myRank)
        {
            compactMap[proci].resize(2*nNonLocal[proci]);
        }
    }


    // Collect all (non-local) elements needed.
    for (const label globalIdx : elements)
    {
        if (globalIdx != -1 && !globalNumbering.isLocal(globalIdx))
        {
            label proci = globalNumbering.whichProcID(globalIdx);
            label index = globalNumbering.toLocal(proci, globalIdx);
            label nCompact = compactMap[proci].size();
            compactMap[proci].insert(index, nCompact);
        }
    }
}


void Foam::mapDistributeBase::calcCompactAddressing
(
    const globalIndex& globalNumbering,
    const labelListList& cellCells,
    List<Map<label>>& compactMap
) const
{
    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    // Count all (non-local) elements needed. Just for presizing map.
    labelList nNonLocal(nProcs, Zero);

    for (const labelList& cCells : cellCells)
    {
        for (const label globalIdx : cCells)
        {
            if (globalIdx != -1 && !globalNumbering.isLocal(globalIdx))
            {
                label proci = globalNumbering.whichProcID(globalIdx);
                nNonLocal[proci]++;
            }
        }
    }

    compactMap.resize_nocopy(nProcs);

    forAll(compactMap, proci)
    {
        compactMap[proci].clear();
        if (proci != myRank)
        {
            compactMap[proci].resize(2*nNonLocal[proci]);
        }
    }


    // Collect all (non-local) elements needed.
    for (const labelList& cCells : cellCells)
    {
        for (const label globalIdx : cCells)
        {
            if (globalIdx != -1 && !globalNumbering.isLocal(globalIdx))
            {
                label proci = globalNumbering.whichProcID(globalIdx);
                label index = globalNumbering.toLocal(proci, globalIdx);
                label nCompact = compactMap[proci].size();
                compactMap[proci].insert(index, nCompact);
            }
        }
    }
}


void Foam::mapDistributeBase::exchangeAddressing
(
    const int tag,
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label>>& compactMap,
    labelList& compactStart
)
{
    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    // The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    compactStart.setSize(nProcs);
    compactStart[myRank] = 0;
    constructSize_ = globalNumbering.localSize();
    forAll(compactStart, proci)
    {
        if (proci != myRank)
        {
            compactStart[proci] = constructSize_;
            constructSize_ += compactMap[proci].size();
        }
    }


    // Find out what to receive/send in compact addressing.

    // What I want to receive is what others have to send
    labelListList wantedRemoteElements(nProcs);
    // Compact addressing for received data
    constructMap_.setSize(nProcs);
    forAll(compactMap, proci)
    {
        if (proci == myRank)
        {
            // All my own elements are used
            label nLocal = globalNumbering.localSize();
            wantedRemoteElements[proci] = identity(nLocal);
            constructMap_[proci] = identity(nLocal);
        }
        else
        {
            // Remote elements wanted from processor proci
            labelList& remoteElem = wantedRemoteElements[proci];
            labelList& localElem = constructMap_[proci];
            remoteElem.setSize(compactMap[proci].size());
            localElem.setSize(compactMap[proci].size());
            label i = 0;
            forAllIters(compactMap[proci], iter)
            {
                const label compactI = compactStart[proci] + iter.val();
                remoteElem[i] = iter.key();
                localElem[i]  = compactI;
                iter.val() = compactI;
                i++;
            }
        }
    }

    subMap_.setSize(nProcs);
    Pstream::exchange<labelList, label>
    (
        wantedRemoteElements,
        subMap_,
        tag,
        comm_
    );

    // Renumber elements
    for (label& elem : elements)
    {
        elem = renumber(globalNumbering, compactMap, elem);
    }
}


void Foam::mapDistributeBase::exchangeAddressing
(
    const int tag,
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label>>& compactMap,
    labelList& compactStart
)
{
    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    // The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    compactStart.setSize(nProcs);
    compactStart[myRank] = 0;
    constructSize_ = globalNumbering.localSize();
    forAll(compactStart, proci)
    {
        if (proci != myRank)
        {
            compactStart[proci] = constructSize_;
            constructSize_ += compactMap[proci].size();
        }
    }


    // Find out what to receive/send in compact addressing.

    // What I want to receive is what others have to send
    labelListList wantedRemoteElements(nProcs);
    // Compact addressing for received data
    constructMap_.setSize(nProcs);
    forAll(compactMap, proci)
    {
        if (proci == myRank)
        {
            // All my own elements are used
            label nLocal = globalNumbering.localSize();
            wantedRemoteElements[proci] = identity(nLocal);
            constructMap_[proci] = identity(nLocal);
        }
        else
        {
            // Remote elements wanted from processor proci
            labelList& remoteElem = wantedRemoteElements[proci];
            labelList& localElem = constructMap_[proci];
            remoteElem.setSize(compactMap[proci].size());
            localElem.setSize(compactMap[proci].size());
            label i = 0;
            forAllIters(compactMap[proci], iter)
            {
                const label compactI = compactStart[proci] + iter.val();
                remoteElem[i] = iter.key();
                localElem[i]  = compactI;
                iter.val() = compactI;
                i++;
            }
        }
    }

    subMap_.setSize(nProcs);
    Pstream::exchange<labelList, label>
    (
        wantedRemoteElements,
        subMap_,
        tag,
        comm_
    );

    // Renumber elements
    for (labelList& cCells : cellCells)
    {
        for (label& celli : cCells)
        {
            celli = renumber(globalNumbering, compactMap, celli);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapDistributeBase::mapDistributeBase()
:
    mapDistributeBase(UPstream::worldComm)
{}


Foam::mapDistributeBase::mapDistributeBase(const label comm)
:
    constructSize_(0),
    subMap_(),
    constructMap_(),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(comm),
    schedulePtr_(nullptr)
{}


Foam::mapDistributeBase::mapDistributeBase(const mapDistributeBase& map)
:
    constructSize_(map.constructSize_),
    subMap_(map.subMap_),
    constructMap_(map.constructMap_),
    subHasFlip_(map.subHasFlip_),
    constructHasFlip_(map.constructHasFlip_),
    comm_(map.comm_),
    schedulePtr_(nullptr)
{}


Foam::mapDistributeBase::mapDistributeBase(mapDistributeBase&& map)
:
    mapDistributeBase(map.comm())
{
    transfer(map);
}


Foam::mapDistributeBase::mapDistributeBase
(
    const label constructSize,
    labelListList&& subMap,
    labelListList&& constructMap,
    const bool subHasFlip,
    const bool constructHasFlip,
    const label comm
)
:
    constructSize_(constructSize),
    subMap_(std::move(subMap)),
    constructMap_(std::move(constructMap)),
    subHasFlip_(subHasFlip),
    constructHasFlip_(constructHasFlip),
    comm_(comm),
    schedulePtr_(nullptr)
{}


Foam::mapDistributeBase::mapDistributeBase
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const label comm
)
:
    constructSize_(0),
    subMap_(),
    constructMap_(),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(comm),
    schedulePtr_(nullptr)
{
    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    if (sendProcs.size() != recvProcs.size())
    {
        FatalErrorInFunction
            << "The send and receive data is not the same length. sendProcs:"
            << sendProcs.size() << " recvProcs:" << recvProcs.size()
            << abort(FatalError);
    }

    // Per processor the number of samples we have to send/receive.
    labelList nSend(nProcs, Zero);
    labelList nRecv(nProcs, Zero);

    forAll(sendProcs, sampleI)
    {
        const label sendProc = sendProcs[sampleI];
        const label recvProc = recvProcs[sampleI];

        // Note that also need to include local communication (both
        // RecvProc and sendProc on local processor)

        if (myRank == sendProc)
        {
            // I am the sender.
            nSend[recvProc]++;
        }
        if (myRank == recvProc)
        {
            // I am the receiver.
            nRecv[sendProc]++;
        }
    }

    subMap_.setSize(nProcs);
    constructMap_.setSize(nProcs);
    forAll(nSend, proci)
    {
        subMap_[proci].setSize(nSend[proci]);
        constructMap_[proci].setSize(nRecv[proci]);
    }
    nSend = 0;
    nRecv = 0;

    // Largest entry inside constructMap
    label maxRecvIndex = -1;

    forAll(sendProcs, sampleI)
    {
        const label sendProc = sendProcs[sampleI];
        const label recvProc = recvProcs[sampleI];

        if (myRank == sendProc)
        {
            // I am the sender. Store index I need to send.
            subMap_[recvProc][nSend[recvProc]++] = sampleI;
        }
        if (myRank == recvProc)
        {
            // I am the receiver.
            constructMap_[sendProc][nRecv[sendProc]++] = sampleI;
            maxRecvIndex = sampleI;
        }
    }

    constructSize_ = maxRecvIndex+1;
}


Foam::mapDistributeBase::mapDistributeBase
(
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label>>& compactMap,
    const int tag,
    const label comm
)
:
    constructSize_(0),
    subMap_(),
    constructMap_(),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(comm),
    schedulePtr_(nullptr)
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        elements,
        compactMap
    );

    //// Sort remote elements needed (not really necessary)
    //forAll(compactMap, proci)
    //{
    //    if (proci != myRank)
    //    {
    //        Map<label>& globalMap = compactMap[proci];
    //
    //        const List<label> sorted(globalMap.sortedToc());
    //
    //        forAll(sorted, i)
    //        {
    //            globalMap(sorted[i]) = i;
    //        }
    //    }
    //}


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        elements,
        compactMap,
        compactStart
    );

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistributeBase::mapDistributeBase
(
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label>>& compactMap,
    const int tag,
    const label comm
)
:
    constructSize_(0),
    subMap_(),
    constructMap_(),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(comm),
    schedulePtr_(nullptr)
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        cellCells,
        compactMap
    );

    //// Sort remote elements needed (not really necessary)
    //forAll(compactMap, proci)
    //{
    //    if (proci != myRank)
    //    {
    //        Map<label>& globalMap = compactMap[proci];
    //
    //        const List<label> sorted(globalMap.sortedToc());
    //
    //        forAll(sorted, i)
    //        {
    //            globalMap(sorted[i]) = i;
    //        }
    //    }
    //}


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        cellCells,
        compactMap,
        compactStart
    );

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistributeBase::mapDistributeBase
(
    labelListList&& subMap,
    const bool subHasFlip,
    const bool constructHasFlip,
    const label comm
)
:
    constructSize_(0),
    subMap_(std::move(subMap)),
    constructMap_(),
    subHasFlip_(subHasFlip),
    constructHasFlip_(constructHasFlip),
    comm_(comm),
    schedulePtr_(nullptr)
{
    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    // Send over how many i need to receive.
    labelList recvSizes;
    Pstream::exchangeSizes(subMap_, recvSizes, comm_);

    // Determine order of receiving
    constructSize_ = 0;
    constructMap_.resize(nProcs);


    // My data first
    {
        const label len = recvSizes[myRank];

        constructMap_[myRank] = identity(len, constructSize_);
        constructSize_ += len;
    }

    // What the other processors are sending to me
    forAll(constructMap_, proci)
    {
        if (proci != myRank)
        {
            const label len = recvSizes[proci];

            constructMap_[proci] = identity(len, constructSize_);
            constructSize_ += len;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::mapDistributeBase::subMapSizes() const
{
    labelList sizes(subMap_.size());
    forAll(subMap_, i)
    {
        sizes[i] = subMap_[i].size();
    }

    return sizes;
}


Foam::labelList Foam::mapDistributeBase::constructMapSizes() const
{
    labelList sizes(constructMap_.size());
    forAll(constructMap_, i)
    {
        sizes[i] = constructMap_[i].size();
    }

    return sizes;
}


void Foam::mapDistributeBase::clear()
{
    constructSize_ = 0;
    subMap_.clear();
    constructMap_.clear();
    subHasFlip_ = false;
    constructHasFlip_ = false;
    // Leave comm_ intact
    schedulePtr_.reset(nullptr);
}


void Foam::mapDistributeBase::transfer(mapDistributeBase& rhs)
{
    if (this == &rhs)
    {
        // Self-assignment is a no-op
        return;
    }

    constructSize_ = rhs.constructSize_;
    subMap_.transfer(rhs.subMap_);
    constructMap_.transfer(rhs.constructMap_);
    subHasFlip_ = rhs.subHasFlip_;
    constructHasFlip_ = rhs.constructHasFlip_;
    comm_ = rhs.comm_;
    schedulePtr_.reset(nullptr);

    rhs.constructSize_ = 0;
    rhs.subHasFlip_ = false;
    rhs.constructHasFlip_ = false;
}


Foam::label Foam::mapDistributeBase::renumber
(
    const globalIndex& globalNumbering,
    const List<Map<label>>& compactMap,
    const label globalI
)
{
    if (globalI == -1)
    {
        return globalI;
    }
    if (globalNumbering.isLocal(globalI))
    {
        return globalNumbering.toLocal(globalI);
    }
    else
    {
        label proci = globalNumbering.whichProcID(globalI);
        label index = globalNumbering.toLocal(proci, globalI);
        return compactMap[proci][index];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::mapDistributeBase::operator=(const mapDistributeBase& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    constructSize_ = rhs.constructSize_;
    subMap_ = rhs.subMap_;
    constructMap_ = rhs.constructMap_;
    subHasFlip_ = rhs.subHasFlip_;
    constructHasFlip_ = rhs.constructHasFlip_;
    comm_ = rhs.comm_;
    schedulePtr_.reset(nullptr);
}


void Foam::mapDistributeBase::operator=(mapDistributeBase&& rhs)
{
    if (this != &rhs)
    {
        // Avoid self assignment
        transfer(rhs);
    }
}


// ************************************************************************* //
