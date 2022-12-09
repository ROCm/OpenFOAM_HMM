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
#include "int.H"
#include "SubList.H"
#include "UPstreamWrapping.H"
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
        // During startup, so worldComm == globalComm

        wordList worlds(numprocs);
        worlds[Pstream::myProcNo(UPstream::globalComm)] = world;
        Pstream::gatherList(worlds, UPstream::msgType(), UPstream::globalComm);

        // Compact
        if (Pstream::master(UPstream::globalComm))
        {
            DynamicList<word> worldNames(numprocs);
            worldIDs_.resize_nocopy(numprocs);

            forAll(worlds, proci)
            {
                const word& world = worlds[proci];

                worldIDs_[proci] = worldNames.find(world);

                if (worldIDs_[proci] == -1)
                {
                    worldIDs_[proci] = worldNames.size();
                    worldNames.push_back(world);
                }
            }

            allWorlds_.transfer(worldNames);
        }
        Pstream::broadcasts(UPstream::globalComm, allWorlds_, worldIDs_);

        const label myWorldId =
            worldIDs_[Pstream::myProcNo(UPstream::globalComm)];

        DynamicList<label> subRanks;
        forAll(worldIDs_, proci)
        {
            if (worldIDs_[proci] == myWorldId)
            {
                subRanks.push_back(proci);
            }
        }

        // Allocate new communicator with globalComm as its parent
        const label subComm =
            UPstream::allocateCommunicator
            (
                UPstream::globalComm,  // parent
                subRanks,
                true
            );


        // Override worldComm
        UPstream::worldComm = subComm;
        // For testing: warn use of non-worldComm
        UPstream::warnComm = UPstream::worldComm;

        if (debug)
        {
            // Check
            int subNumProcs, subRank;
            MPI_Comm_size
            (
                PstreamGlobals::MPICommunicators_[subComm],
                &subNumProcs
            );
            MPI_Comm_rank
            (
                PstreamGlobals::MPICommunicators_[subComm],
                &subRank
            );

            Pout<< "UPstream::init : in world:" << world
                << " using local communicator:" << subComm
                << " rank " << subRank
                << " of " << subNumProcs
                << endl;
        }

        // Override Pout prefix (move to setParRun?)
        Pout.prefix() = '[' + world + '/' +  name(myProcNo(subComm)) + "] ";
        Perr.prefix() = Pout.prefix();
    }
    else
    {
        // All processors use world 0
        worldIDs_.resize_nocopy(numprocs);
        worldIDs_ = 0;
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
        if (myProcNo_[communicator] >= 0)
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
        MPI_Comm newComm = MPI_COMM_NULL;
        MPI_Group newGroup = MPI_GROUP_NULL;
        PstreamGlobals::MPIGroups_.push_back(newGroup);
        PstreamGlobals::MPICommunicators_.push_back(newComm);
    }
    else if (index > PstreamGlobals::MPIGroups_.size())
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::exit(FatalError);
    }


    if (parentIndex == -1)
    {
        // Global communicator. Same as world communicator for single-world

        if (index != UPstream::globalComm)
        {
            FatalErrorInFunction
                << "world communicator should always be index "
                << UPstream::globalComm
                << Foam::exit(FatalError);
        }

        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_WORLD;
        MPI_Comm_group(MPI_COMM_WORLD, &PstreamGlobals::MPIGroups_[index]);
        MPI_Comm_rank(MPI_COMM_WORLD, &myProcNo_[index]);

        // Set the number of ranks to the actual number
        int numProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        //procIDs_[index] = identity(numProcs);
        procIDs_[index].resize_nocopy(numProcs);
        forAll(procIDs_[index], i)
        {
            procIDs_[index][i] = i;
        }
    }
    else if (parentIndex == -2)
    {
        // Self communicator

        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_SELF;
        MPI_Comm_group(MPI_COMM_SELF, &PstreamGlobals::MPIGroups_[index]);
        MPI_Comm_rank(MPI_COMM_SELF, &myProcNo_[index]);

        // Number of ranks is always 1 (self communicator)

        #ifdef FULLDEBUG
        int numProcs;
        MPI_Comm_size(MPI_COMM_SELF, &numProcs);

        if (numProcs != 1)
        {
            // Already finalized - this is an error
            FatalErrorInFunction
                << "MPI_COMM_SELF had " << numProcs << " != 1 ranks!\n"
                << Foam::abort(FatalError);
        }
        #endif

        procIDs_[index].resize_nocopy(1);
        procIDs_[index] = 0;
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
            UPstream::msgType(),
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
    // Skip placeholders and pre-defined (not allocated) communicators

    if (UPstream::debug)
    {
        Pout<< "freePstreamCommunicator: " << communicator
            << " from " << PstreamGlobals::MPICommunicators_.size() << endl;
    }

    // Not touching the first two communicators (SELF, WORLD)
    if (communicator > 1)
    {
        if (MPI_COMM_NULL != PstreamGlobals::MPICommunicators_[communicator])
        {
            // Free communicator. Sets communicator to MPI_COMM_NULL
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[communicator]);
        }

        if (MPI_GROUP_NULL != PstreamGlobals::MPIGroups_[communicator])
        {
            // Free group. Sets group to MPI_GROUP_NULL
            MPI_Group_free(&PstreamGlobals::MPIGroups_[communicator]);
        }
    }
}


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
    if (!UPstream::parRun())
    {
        return;  // No-op for non-parallel
    }

    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << PstreamGlobals::outstandingRequests_.size()-start
            << " outstanding requests starting at " << start << endl;
    }

    // TBD: check for
    // (start < 0 || start > PstreamGlobals::outstandingRequests_.size())

    if (PstreamGlobals::outstandingRequests_.size())
    {
        SubList<MPI_Request> waitRequests
        (
            PstreamGlobals::outstandingRequests_,
            PstreamGlobals::outstandingRequests_.size() - start,
            start
        );

        profilingPstream::beginTiming();

        // On success: sets each request to MPI_REQUEST_NULL
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

        PstreamGlobals::outstandingRequests_.resize(start);
    }

    if (debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequest(const label i)
{
    if (!UPstream::parRun() || i < 0)
    {
        return;  // No-op for non-parallel, or placeholder indices
    }

    if (debug)
    {
        Pout<< "UPstream::waitRequest : starting wait for request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "You asked for request=" << i
            << " from " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding requests!" << nl
            << "Mixing use of blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    profilingPstream::beginTiming();

    // On success: sets request to MPI_REQUEST_NULL
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
    PstreamGlobals::freedRequests_.push_back(i);

    if (debug)
    {
        Pout<< "UPstream::waitRequest : finished wait for request:" << i
            << endl;
    }
}


bool Foam::UPstream::finishedRequest(const label i)
{
    if (!UPstream::parRun() || i < 0)
    {
        return true;  // No-op for non-parallel, or placeholder indices
    }

    if (debug)
    {
        Pout<< "UPstream::finishedRequest : checking request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "You asked for request=" << i
            << " from " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding requests!" << nl
            << "Mixing use of blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    // On success: sets request to MPI_REQUEST_NULL
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


int Foam::UPstream::allocateTag(const char* const msg)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.back();
        PstreamGlobals::freedTags_.pop_back();
    }
    else
    {
        tag = ++PstreamGlobals::nTags_;
    }

    if (debug)
    {
        Pout<< "UPstream::allocateTag";
        if (msg) Pout<< ' ' << msg;
        Pout<< " : tag:" << tag << endl;
    }

    return tag;
}


void Foam::UPstream::freeTag(const int tag, const char* const msg)
{
    if (debug)
    {
        Pout<< "UPstream::freeTag ";
        if (msg) Pout<< ' ' << msg;
        Pout<< " : tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.push_back(tag);
}


// ************************************************************************* //
