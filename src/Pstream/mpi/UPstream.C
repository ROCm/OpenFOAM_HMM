/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"
#include "SubList.H"
#include "UPstreamWrapping.H"
#include "int.H"
#include "collatedFileOperation.H"

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <csignal>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The min value and default for MPI buffers length
constexpr int minBufLen = 20000000;

// Track if we have attached MPI buffers
static bool ourBuffers = false;

// Track if we initialized MPI
static bool ourMpi = false;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static void attachOurBuffers()
{
    if (ourBuffers)
    {
        return;  // Already attached
    }
    ourBuffers = true;

    // Use UPstream::mpiBufferSize (optimisationSwitch),
    // but allow override with MPI_BUFFER_SIZE env variable (int value)

#ifndef SGIMPI
    int len = 0;

    const std::string str(Foam::getEnv("MPI_BUFFER_SIZE"));
    if (str.empty() || !Foam::read(str, len) || len <= 0)
    {
        len = Foam::UPstream::mpiBufferSize;
    }

    if (len < minBufLen)
    {
        len = minBufLen;
    }

    if (Foam::UPstream::debug)
    {
        Foam::Pout<< "UPstream::init : buffer-size " << len << '\n';
    }

    char* buf = new char[len];

    if (MPI_SUCCESS != MPI_Buffer_attach(buf, len))
    {
        delete[] buf;
        Foam::Pout<< "UPstream::init : could not attach buffer\n";
    }
#endif
}


static void detachOurBuffers()
{
    if (!ourBuffers)
    {
        return;  // Nothing to detach
    }
    ourBuffers = false;

    // Some MPI notes suggest that the return code is MPI_SUCCESS when
    // no buffer is attached.
    // Be extra careful and require a non-zero size as well.

#ifndef SGIMPI
    int len = 0;
    char* buf = nullptr;

    if (MPI_SUCCESS == MPI_Buffer_detach(&buf, &len) && len)
    {
        delete[] buf;
    }
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::initNull()
{
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init\n"
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        if (debug)
        {
            Pout<< "UPstream::initNull : was already initialized\n";
        }
    }
    else
    {
        // Not already initialized

        MPI_Init_thread
        (
            nullptr,    // argc
            nullptr,    // argv
            MPI_THREAD_SINGLE,
            &flag       // provided_thread_support
        );

        ourMpi = true;
    }

    // Could also attach buffers etc.

    return true;
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    int numprocs = 0, myRank = 0;
    int provided_thread_support = 0;
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init" << endl
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        // Already initialized.
        // Warn if we've called twice, but skip if initialized externally

        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already initialized - cannot perform MPI_Init" << nl
                << "This could indicate an application programming error!"
                << endl;

            return true;
        }
        else if (debug)
        {
            Pout<< "UPstream::init : was already initialized\n";
        }
    }
    else
    {
        MPI_Init_thread
        (
            &argc,
            &argv,
            (
                needsThread
              ? MPI_THREAD_MULTIPLE
              : MPI_THREAD_SINGLE
            ),
            &provided_thread_support
        );

        ourMpi = true;
    }

    // Check argument list for local world
    label worldIndex = -1;
    word world;
    for (int argi = 1; argi < argc; ++argi)
    {
        if (strcmp(argv[argi], "-world") == 0)
        {
            worldIndex = argi++;
            if (argi >= argc)
            {
                FatalErrorInFunction
                    << "Missing world name to argument \"world\""
                    << Foam::abort(FatalError);
            }
            world = argv[argi];
            break;
        }
    }

    // Filter 'world' option
    if (worldIndex != -1)
    {
        for (label i = worldIndex+2; i < argc; i++)
        {
            argv[i-2] = argv[i];
        }
        argc -= 2;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (debug)
    {
        Pout<< "UPstream::init :"
            << " thread-support : wanted:" << needsThread
            << " obtained:"
            <<  (
                    provided_thread_support == MPI_THREAD_MULTIPLE
                  ? "MPI_THREAD_MULTIPLE"
                  : "MPI_THREAD_SINGLE"
                )
            << " procs:" << numprocs
            << " rank:" << myRank
            << " world:" << world << endl;
    }

    if (worldIndex == -1 && numprocs <= 1)
    {
        FatalErrorInFunction
            << "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }

    // Initialise parallel structure
    setParRun(numprocs, provided_thread_support == MPI_THREAD_MULTIPLE);

    if (worldIndex != -1)
    {
        wordList worlds(numprocs);
        worlds[Pstream::myProcNo()] = world;
        Pstream::gatherList(worlds);
        Pstream::broadcast(worlds);

        // Compact
        if (Pstream::master())
        {
            DynamicList<word> allWorlds(numprocs);
            for (const word& world : worlds)
            {
                allWorlds.appendUniq(world);
            }
            allWorlds_ = std::move(allWorlds);

            worldIDs_.setSize(numprocs);
            forAll(worlds, proci)
            {
                const word& world = worlds[proci];
                worldIDs_[proci] = allWorlds_.find(world);
            }
        }
        Pstream::broadcasts(UPstream::worldComm, allWorlds_, worldIDs_);

        DynamicList<label> subRanks;
        forAll(worlds, proci)
        {
            if (worlds[proci] == worlds[Pstream::myProcNo()])
            {
                subRanks.append(proci);
            }
        }

        // Allocate new communicator 1 with parent 0 (= mpi_world)
        const label subComm = allocateCommunicator(0, subRanks, true);

        // Override worldComm
        UPstream::worldComm = subComm;
        // For testing: warn use of non-worldComm
        UPstream::warnComm = UPstream::worldComm;

        if (debug)
        {
            // Check
            int subNProcs, subRank;
            MPI_Comm_size
            (
                PstreamGlobals::MPICommunicators_[subComm],
                &subNProcs
            );
            MPI_Comm_rank
            (
                PstreamGlobals::MPICommunicators_[subComm],
                &subRank
            );

            Pout<< "UPstream::init : in world:" << world
                << " using local communicator:" << subComm
                << " with procs:" << subNProcs
                << " and rank:" << subRank
                << endl;
        }

        // Override Pout prefix (move to setParRun?)
        Pout.prefix() = '[' + world + '/' +  name(myProcNo(subComm)) + "] ";
        Perr.prefix() = '[' + world + '/' +  name(myProcNo(subComm)) + "] ";
    }
    else
    {
        // All processors use world 0
        worldIDs_.setSize(numprocs, 0);
    }

    attachOurBuffers();

    return true;
}


void Foam::UPstream::shutdown(int errNo)
{
    if (debug)
    {
        Pout<< "UPstream::shutdown\n";
    }

    int flag = 0;

    MPI_Initialized(&flag);
    if (!flag)
    {
        // No MPI initialized - we are done
        return;
    }

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized elsewhere?
        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already finalized (by a connected program?)\n";
        }
        else if (debug)
        {
            Pout<< "UPstream::shutdown : was already finalized\n";
        }
    }
    else
    {
        detachOurBuffers();
    }


    // Warn about any outstanding requests
    {
        label nOutstanding = 0;

        forAll(PstreamGlobals::outstandingRequests_, requestID)
        {
            if (!PstreamGlobals::freedRequests_.found(requestID))
            {
                ++nOutstanding;
            }
        }

        PstreamGlobals::outstandingRequests_.clear();

        if (nOutstanding)
        {
            WarningInFunction
                << "There were still " << nOutstanding
                << " outstanding MPI_Requests." << nl
                << "Which means your code exited before doing a "
                << " UPstream::waitRequests()." << nl
                << "This should not happen for a normal code exit."
                << nl;
        }
    }

    // Clean mpi communicators
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freePstreamCommunicator(communicator);
        }
    }

    if (!flag)
    {
        // MPI not already finalized

        if (!ourMpi)
        {
            WarningInFunction
                << "Finalizing MPI, but was initialized elsewhere\n";
        }

        if (errNo == 0)
        {
            MPI_Finalize();
        }
        else
        {
            // Abort only locally or world?
            MPI_Abort(MPI_COMM_WORLD, errNo);
        }
    }
}


void Foam::UPstream::exit(int errNo)
{
    UPstream::shutdown(errNo);
    std::exit(errNo);
}


void Foam::UPstream::abort()
{
    MPI_Abort(MPI_COMM_WORLD, 1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::allocatePstreamCommunicator
(
    const label parentIndex,
    const label index
)
{
    if (index == PstreamGlobals::MPIGroups_.size())
    {
        // Extend storage with dummy values
        MPI_Group newGroup = MPI_GROUP_NULL;
        PstreamGlobals::MPIGroups_.append(newGroup);
        MPI_Comm newComm = MPI_COMM_NULL;
        PstreamGlobals::MPICommunicators_.append(newComm);
    }
    else if (index > PstreamGlobals::MPIGroups_.size())
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::exit(FatalError);
    }


    if (parentIndex == -1)
    {
        // Allocate world communicator

        if (index != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "world communicator should always be index "
                << UPstream::worldComm << Foam::exit(FatalError);
        }

        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_WORLD;
        MPI_Comm_group(MPI_COMM_WORLD, &PstreamGlobals::MPIGroups_[index]);
        MPI_Comm_rank
        (
            PstreamGlobals::MPICommunicators_[index],
           &myProcNo_[index]
        );

        // Set the number of processes to the actual number
        int numProcs;
        MPI_Comm_size(PstreamGlobals::MPICommunicators_[index], &numProcs);

        //procIDs_[index] = identity(numProcs);
        procIDs_[index].setSize(numProcs);
        forAll(procIDs_[index], i)
        {
            procIDs_[index][i] = i;
        }
    }
    else
    {
        // Create new group
        MPI_Group_incl
        (
            PstreamGlobals::MPIGroups_[parentIndex],
            procIDs_[index].size(),
            procIDs_[index].cdata(),
           &PstreamGlobals::MPIGroups_[index]
        );

        #if defined(MSMPI_VER)
        // ms-mpi (10.0 and others?) does not have MPI_Comm_create_group
        MPI_Comm_create
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            PstreamGlobals::MPIGroups_[index],
            &PstreamGlobals::MPICommunicators_[index]
        );
        #else
        // Create new communicator for this group
        MPI_Comm_create_group
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            PstreamGlobals::MPIGroups_[index],
            Pstream::msgType(),
           &PstreamGlobals::MPICommunicators_[index]
        );
        #endif

        if (PstreamGlobals::MPICommunicators_[index] == MPI_COMM_NULL)
        {
            myProcNo_[index] = -1;
        }
        else
        {
            if
            (
                MPI_Comm_rank
                (
                    PstreamGlobals::MPICommunicators_[index],
                   &myProcNo_[index]
                )
            )
            {
                FatalErrorInFunction
                    << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << procIDs_[index]
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::UPstream::freePstreamCommunicator(const label communicator)
{
    if (communicator != 0)
    {
        if (PstreamGlobals::MPICommunicators_[communicator] != MPI_COMM_NULL)
        {
            // Free communicator. Sets communicator to MPI_COMM_NULL
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[communicator]);
        }
        if (PstreamGlobals::MPIGroups_[communicator] != MPI_GROUP_NULL)
        {
            // Free greoup. Sets group to MPI_GROUP_NULL
            MPI_Group_free(&PstreamGlobals::MPIGroups_[communicator]);
        }
    }
}


Foam::label Foam::UPstream::nRequests()
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label i)
{
    if (i < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.setSize(i);
    }
}


void Foam::UPstream::waitRequests(const label start)
{
    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << PstreamGlobals::outstandingRequests_.size()-start
            << " outstanding requests starting at " << start << endl;
    }

    if (PstreamGlobals::outstandingRequests_.size())
    {
        SubList<MPI_Request> waitRequests
        (
            PstreamGlobals::outstandingRequests_,
            PstreamGlobals::outstandingRequests_.size() - start,
            start
        );

        profilingPstream::beginTiming();

        if
        (
            MPI_Waitall
            (
                waitRequests.size(),
                waitRequests.data(),
                MPI_STATUSES_IGNORE
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Waitall returned with error" << Foam::endl;
        }

        profilingPstream::addWaitTime();

        resetRequests(start);
    }

    if (debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequest : starting wait for request:" << i
            << endl;
    }

    if (i < 0 || i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    profilingPstream::beginTiming();

    if
    (
        MPI_Wait
        (
           &PstreamGlobals::outstandingRequests_[i],
            MPI_STATUS_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error" << Foam::endl;
    }

    profilingPstream::addWaitTime();
    // Push index onto free cache
    PstreamGlobals::freedRequests_.append(i);

    if (debug)
    {
        Pout<< "UPstream::waitRequest : finished wait for request:" << i
            << endl;
    }
}


bool Foam::UPstream::finishedRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::finishedRequest : checking request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    int flag;
    MPI_Test
    (
       &PstreamGlobals::outstandingRequests_[i],
       &flag,
        MPI_STATUS_IGNORE
    );

    if (debug)
    {
        Pout<< "UPstream::finishedRequest : finished request:" << i
            << endl;
    }

    return flag != 0;
}


int Foam::UPstream::allocateTag(const char* s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        Pout<< "UPstream::allocateTag "
            << s << " : tag:" << tag << endl;
    }

    return tag;
}


int Foam::UPstream::allocateTag(const std::string& s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        Pout<< "UPstream::allocateTag "
            << s.c_str() << " : tag:" << tag << endl;
    }

    return tag;
}


void Foam::UPstream::freeTag(const char* s, const int tag)
{
    if (debug)
    {
        Pout<< "UPstream::freeTag "
            << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


void Foam::UPstream::freeTag(const std::string& s, const int tag)
{
    if (debug)
    {
        Pout<< "UPstream::freeTag "
            << s.c_str() << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


// ************************************************************************* //
