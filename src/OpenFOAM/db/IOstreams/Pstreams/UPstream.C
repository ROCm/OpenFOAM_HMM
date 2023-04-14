/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

Note
    Included by global/globals.C

\*---------------------------------------------------------------------------*/

#include "UPstream.H"
#include "debug.H"
#include "registerSwitch.H"
#include "dictionary.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPstream, 0);
}

const Foam::Enum
<
    Foam::UPstream::commsTypes
>
Foam::UPstream::commsTypeNames
({
    { commsTypes::blocking, "blocking" },
    { commsTypes::scheduled, "scheduled" },
    { commsTypes::nonBlocking, "nonBlocking" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UPstream::setParRun(const label nProcs, const bool haveThreads)
{
    parRun_ = (nProcs > 0);
    haveThreads_ = haveThreads;

    label comm = -1;

    if (!parRun_)
    {
        // These are already correct from the static initialisation,
        // but just in case of future changes

        // Using (world, self) ordering
        freeCommunicator(UPstream::selfComm);
        freeCommunicator(UPstream::globalComm);

        // 0: worldComm
        comm = allocateCommunicator(-1, Foam::labelList(Foam::one{}, 0), false);
        if (comm != UPstream::globalComm)
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::globalComm:" << UPstream::globalComm
                << Foam::exit(FatalError);
        }

        // 1: selfComm
        comm = allocateCommunicator(-2, Foam::labelList(Foam::one{}, 0), false);
        if (comm != UPstream::selfComm)
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::selfComm:" << UPstream::selfComm
                << Foam::exit(FatalError);
        }

        Pout.prefix().clear();
        Perr.prefix().clear();
    }
    else
    {
        // Redo communicators that were created during static initialisation
        // but this time with Pstream components

        // Using (world, self) ordering
        freeCommunicator(UPstream::selfComm);
        freeCommunicator(UPstream::globalComm);

        // 0: worldComm
        comm = allocateCommunicator(-1, identity(nProcs), true);
        if (comm != UPstream::globalComm)
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::globalComm:" << UPstream::globalComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = '[' +  Foam::name(myProcNo(comm)) + "] ";
        Perr.prefix() = Pout.prefix();

        // 1: selfComm
        comm = allocateCommunicator(-2, Foam::labelList(Foam::one{}, 0), true);
        if (comm != UPstream::selfComm)
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::selfComm:" << UPstream::selfComm
                << Foam::exit(FatalError);
        }
    }

    if (debug)
    {
        Pout<< "UPstream::setParRun :"
            << " nProcs:" << nProcs
            << " haveThreads:" << haveThreads
            << endl;
    }
}


Foam::label Foam::UPstream::getAvailableCommIndex(const label parentIndex)
{
    label index;
    if (!freeComms_.empty())
    {
        // LIFO pop
        index = freeComms_.back();
        freeComms_.pop_back();

        // Reset existing
        myProcNo_[index] = -1;
        parentComm_[index] = parentIndex;

        procIDs_[index].clear();
        linearCommunication_[index].clear();
        treeCommunication_[index].clear();
    }
    else
    {
        // Extend storage
        index = parentComm_.size();

        myProcNo_.push_back(-1);
        parentComm_.push_back(parentIndex);

        procIDs_.emplace_back();
        linearCommunication_.emplace_back();
        treeCommunication_.emplace_back();
    }

    return index;
}


Foam::label Foam::UPstream::allocateCommunicator
(
    const label parentIndex,
    const labelUList& subRanks,
    const bool doPstream
)
{
    const label index = getAvailableCommIndex(parentIndex);

    if (debug)
    {
        Pout<< "Allocating communicator " << index << nl
            << "    parent : " << parentIndex << nl
            << "    procs  : " << subRanks << nl
            << endl;
    }

    // Initially treat as master, overwritten by allocatePstreamCommunicator
    myProcNo_[index] = UPstream::masterNo();

    // The selected sub-ranks.
    // - transcribe from label to int. Treat negative values as 'ignore'
    // - enforce incremental order (so index is rank in next communicator)

    auto& procIds = procIDs_[index];
    procIds.resize_nocopy(subRanks.size());

    label numSubRanks = 0;
    bool monotonicOrder = true;
    for (const label subRanki : subRanks)
    {
        if (subRanki < 0)
        {
            continue;
        }
        if (monotonicOrder && numSubRanks)
        {
            monotonicOrder = (procIds[numSubRanks-1] < subRanki);
        }

        procIds[numSubRanks] = subRanki;
        ++numSubRanks;
    }

    if (!monotonicOrder)
    {
        auto last = procIds.begin() + numSubRanks;
        std::sort(procIds.begin(), last);
        last = std::unique(procIds.begin(), last);
        numSubRanks = label(last - procIds.begin());
    }

    procIds.resize(numSubRanks);

    // Sizing and filling are demand-driven
    linearCommunication_[index].clear();
    treeCommunication_[index].clear();

    if (doPstream && parRun())
    {
        allocatePstreamCommunicator(parentIndex, index);

        // Could 'remember' locations of uninvolved ranks
        /// if (myProcNo_[index] < 0 && parentIndex >= 0)
        /// {
        ///     // As global rank
        ///     myProcNo_[index] = -(myProcNo_[worldComm]+1);
        ///
        /// OR:
        ///     // As parent rank number
        ///     if (myProcNo_[parentIndex] >= 0)
        ///     {
        ///         myProcNo_[index] = -(myProcNo_[parentIndex]+1);
        ///     }
        /// }
    }

    return index;
}


void Foam::UPstream::freeCommunicator
(
    const label communicator,
    const bool doPstream
)
{
    // Filter out any placeholders
    if (communicator < 0)
    {
        return;
    }

    if (debug)
    {
        Pout<< "Communicators : Freeing communicator " << communicator
            << " parent: " << parentComm_[communicator]
            << " myProcNo: " << myProcNo_[communicator]
            << endl;
    }

    if (doPstream && parRun())
    {
        freePstreamCommunicator(communicator);
    }

    myProcNo_[communicator] = -1;
    parentComm_[communicator] = -1;
    //procIDs_[communicator].clear();
    linearCommunication_[communicator].clear();
    treeCommunication_[communicator].clear();

    // LIFO push
    freeComms_.push_back(communicator);
}


void Foam::UPstream::freeCommunicators(const bool doPstream)
{
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] >= 0)
        {
            freeCommunicator(communicator, doPstream);
        }
    }
}


int Foam::UPstream::baseProcNo(label comm, int procID)
{
    while (UPstream::parent(comm) >= 0 && procID >= 0)
    {
        const auto& parentRanks = UPstream::procID(comm);
        procID = parentRanks[procID];
        comm = parent(comm);
    }

    return procID;
}


Foam::label Foam::UPstream::procNo(const label comm, const int baseProcID)
{
    const auto& parentRanks = UPstream::procID(comm);
    label parentComm = UPstream::parent(comm);

    int procID = baseProcID;

    if (parentComm >= 0)
    {
        procID = UPstream::procNo(parentComm, baseProcID);
    }

    return parentRanks.find(procID);
}


Foam::label Foam::UPstream::procNo
(
    const label comm,
    const label currentComm,
    const int currentProcID
)
{
    label physProcID = UPstream::baseProcNo(currentComm, currentProcID);
    return UPstream::procNo(comm, physProcID);
}


const Foam::List<Foam::UPstream::commsStruct>&
Foam::UPstream::linearCommunication(const label communicator)
{
    if (linearCommunication_[communicator].empty())
    {
        linearCommunication_[communicator] =
            List<commsStruct>(UPstream::nProcs(communicator));
    }

    return linearCommunication_[communicator];
}


const Foam::List<Foam::UPstream::commsStruct>&
Foam::UPstream::treeCommunication(const label communicator)
{
    if (treeCommunication_[communicator].empty())
    {
        treeCommunication_[communicator] =
            List<commsStruct>(UPstream::nProcs(communicator));
    }

    return treeCommunication_[communicator];
}


void Foam::UPstream::printCommTree(const label communicator)
{
    const auto& comms = UPstream::whichCommunication(communicator);

    if (UPstream::master(communicator))
    {
        commsStruct::printGraph(Info(), comms);
    }
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::UPstream::parRun_(false);

bool Foam::UPstream::haveThreads_(false);

int Foam::UPstream::msgType_(1);


Foam::wordList Foam::UPstream::allWorlds_(Foam::one{}, "");
Foam::labelList Foam::UPstream::worldIDs_(Foam::one{}, 0);


Foam::DynamicList<int> Foam::UPstream::myProcNo_(16);
Foam::DynamicList<Foam::List<int>> Foam::UPstream::procIDs_(16);

Foam::DynamicList<Foam::label> Foam::UPstream::parentComm_(16);
Foam::DynamicList<Foam::label> Foam::UPstream::freeComms_;

Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::linearCommunication_(16);

Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::treeCommunication_(16);


Foam::label Foam::UPstream::worldComm(0);
Foam::label Foam::UPstream::warnComm(-1);


// Predefine worldComm, selfComm slots.
// These are overwritten in parallel mode (by UPstream::setParRun())
const Foam::label nPredefinedComm = []()
{
    const Foam::labelList singleProc(Foam::one{}, 0);

    // 0: worldComm
    (void) Foam::UPstream::allocateCommunicator(-1, singleProc, false);

    // 1: selfComm
    (void) Foam::UPstream::allocateCommunicator(-2, singleProc, false);

    return Foam::UPstream::nComms();
}();


bool Foam::UPstream::floatTransfer
(
    Foam::debug::optimisationSwitch("floatTransfer", 0)
);
registerOptSwitch
(
    "floatTransfer",
    bool,
    Foam::UPstream::floatTransfer
);

int Foam::UPstream::nProcsSimpleSum
(
    Foam::debug::optimisationSwitch("nProcsSimpleSum", 0)
);
registerOptSwitch
(
    "nProcsSimpleSum",
    int,
    Foam::UPstream::nProcsSimpleSum
);

int Foam::UPstream::nProcsNonblockingExchange
(
    Foam::debug::optimisationSwitch("nonBlockingExchange", 0)
);
registerOptSwitch
(
    "nonBlockingExchange",
    int,
    Foam::UPstream::nProcsNonblockingExchange
);


int Foam::UPstream::nPollProcInterfaces
(
    Foam::debug::optimisationSwitch("nPollProcInterfaces", 0)
);
registerOptSwitch
(
    "nPollProcInterfaces",
    int,
    Foam::UPstream::nPollProcInterfaces
);


Foam::UPstream::commsTypes Foam::UPstream::defaultCommsType
(
    commsTypeNames.get
    (
        "commsType",
        Foam::debug::optimisationSwitches()
    )
);


//! \cond file_scope
namespace Foam
{
    //- Registered reader for UPstream::defaultCommsType
    class addcommsTypeToOpt
    :
        public ::Foam::simpleRegIOobject
    {
    public:

        addcommsTypeToOpt(const char* name)
        :
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
        {}

        virtual ~addcommsTypeToOpt() = default;

        virtual void readData(Foam::Istream& is)
        {
            UPstream::defaultCommsType =
                UPstream::commsTypeNames.read(is);
        }

        virtual void writeData(Foam::Ostream& os) const
        {
            os << UPstream::commsTypeNames[UPstream::defaultCommsType];
        }
    };

    addcommsTypeToOpt addcommsTypeToOpt_("commsType");
}
//! \endcond

int Foam::UPstream::maxCommsSize
(
    Foam::debug::optimisationSwitch("maxCommsSize", 0)
);
registerOptSwitch
(
    "maxCommsSize",
    int,
    Foam::UPstream::maxCommsSize
);


const int Foam::UPstream::mpiBufferSize
(
    Foam::debug::optimisationSwitch("mpiBufferSize", 0)
);


// ************************************************************************* //
