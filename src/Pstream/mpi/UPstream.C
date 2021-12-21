/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
#include "allReduce.H"
#include "int.H"
#include "collatedFileOperation.H"

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <csignal>

#if defined(WM_SP)
    #define MPI_SCALAR MPI_FLOAT
    #define MPI_SOLVESCALAR MPI_FLOAT
#elif defined(WM_SPDP)
    #define MPI_SCALAR MPI_FLOAT
    #define MPI_SOLVESCALAR MPI_DOUBLE
#elif defined(WM_DP)
    #define MPI_SCALAR MPI_DOUBLE
    #define MPI_SOLVESCALAR MPI_DOUBLE
#endif

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
        Pstream::scatterList(worlds);

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
        Pstream::scatter(allWorlds_);
        Pstream::scatter(worldIDs_);

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


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 1, MPI_SCALAR, MPI_SUM, bop, tag, communicator);
}


void Foam::reduce
(
    scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 1, MPI_SCALAR, MPI_MIN, bop, tag, communicator);
}


void Foam::reduce
(
    vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 2, MPI_SCALAR, MPI_SUM, bop, tag, communicator);
}


void Foam::sumReduce
(
    scalar& Value,
    label& Count,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** sumReduce:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    vector2D twoScalars(Value, scalar(Count));
    reduce(twoScalars, sumOp<vector2D>(), tag, communicator);

    Value = twoScalars.x();
    Count = twoScalars.y();
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
    iallReduce<scalar>(&Value, 1, MPI_SCALAR, MPI_SUM, communicator, requestID);
}


void Foam::reduce
(
    scalar values[],
    const int size,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
    iallReduce<scalar>
    (
        values,
        size,
        MPI_SCALAR,
        MPI_SUM,
        communicator,
        requestID
    );
}


#if defined(WM_SPDP)
void Foam::reduce
(
    solveScalar& Value,
    const sumOp<solveScalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 1, MPI_SOLVESCALAR, MPI_SUM, bop, tag, communicator);
}


void Foam::reduce
(
    solveScalar& Value,
    const minOp<solveScalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 1, MPI_SOLVESCALAR, MPI_MIN, bop, tag, communicator);
}


void Foam::reduce
(
    Vector2D<solveScalar>& Value,
    const sumOp<Vector2D<solveScalar>>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 2, MPI_SOLVESCALAR, MPI_SUM, bop, tag, communicator);
}


void Foam::sumReduce
(
    solveScalar& Value,
    label& Count,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    Vector2D<solveScalar> twoScalars(Value, solveScalar(Count));
    reduce(twoScalars, sumOp<Vector2D<solveScalar>>(), tag, communicator);

    Value = twoScalars.x();
    Count = twoScalars.y();
}


void Foam::reduce
(
    solveScalar& Value,
    const sumOp<solveScalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
    iallReduce<solveScalar>
    (
        &Value,
        1,
        MPI_SOLVESCALAR,
        MPI_SUM,
        communicator,
        requestID
    );
}


void Foam::reduce
(
    solveScalar values[],
    const int size,
    const sumOp<solveScalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
    iallReduce<solveScalar>
    (
        values,
        size,
        MPI_SOLVESCALAR,
        MPI_SUM,
        communicator,
        requestID
    );
}
#endif


void Foam::UPstream::allToAll
(
    const labelUList& sendData,
    labelUList& recvData,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** allToAll :"
            << " np:" << np
            << " sendData:" << sendData.size()
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if (sendData.size() != np || recvData.size() != np)
    {
        FatalErrorInFunction
            << "Size of sendData " << sendData.size()
            << " or size of recvData " << recvData.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        recvData.deepCopy(sendData);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Alltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                const_cast<label*>(sendData.cdata()),
                sizeof(label),
                MPI_BYTE,
                recvData.data(),
                sizeof(label),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoall failed for " << sendData
                << " on communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addAllToAllTime();
    }
}


void Foam::UPstream::allToAll
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,

    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Alltoallv :"
            << " sendSizes:" << sendSizes
            << " sendOffsets:" << sendOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        sendSizes.size() != np
     || sendOffsets.size() != np
     || recvSizes.size() != np
     || recvOffsets.size() != np
    )
    {
        FatalErrorInFunction
            << "Size of sendSize " << sendSizes.size()
            << ", sendOffsets " << sendOffsets.size()
            << ", recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        if (recvSizes[0] != sendSizes[0])
        {
            FatalErrorInFunction
                << "Bytes to send " << sendSizes[0]
                << " does not equal bytes to receive " << recvSizes[0]
                << Foam::abort(FatalError);
        }
        std::memmove(recvData, &sendData[sendOffsets[0]], recvSizes[0]);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Alltoallv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoallv failed for sendSizes " << sendSizes
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addAllToAllTime();
    }
}


void Foam::UPstream::mpiGather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Gather :"
            << " np:" << np
            << " recvSize:" << recvSize
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvSize);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Gather
            (
                const_cast<char*>(sendData),
                sendSize,
                MPI_BYTE,
                recvData,
                recvSize,
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gather failed for sendSize " << sendSize
                << " recvSize " << recvSize
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addGatherTime();
    }
}


void Foam::UPstream::mpiScatter
(
    const char* sendData,
    int sendSize,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Scatter :"
            << " np:" << np
            << " recvSize:" << recvSize
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvSize);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Scatter
            (
                const_cast<char*>(sendData),
                sendSize,
                MPI_BYTE,
                recvData,
                recvSize,
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatter failed for sendSize " << sendSize
                << " recvSize " << recvSize
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addScatterTime();
    }
}


void Foam::UPstream::gather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Gatherv :"
            << " np:" << np
            << " recvSizes:" << recvSizes
            << " recvOffsets:" << recvOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        UPstream::master(communicator)
     && (recvSizes.size() != np || recvOffsets.size() < np)
    )
    {
        // Note: allow recvOffsets to be e.g. 1 larger than np so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Size of recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        // recvSizes[0] may be invalid - use sendSize instead
        std::memmove(recvData, sendData, sendSize);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Gatherv
            (
                const_cast<char*>(sendData),
                sendSize,
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gatherv failed for sendSize " << sendSize
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addGatherTime();
    }
}


void Foam::UPstream::scatter
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Scatterv :"
            << " np:" << np
            << " sendSizes:" << sendSizes
            << " sendOffsets:" << sendOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        UPstream::master(communicator)
     && (sendSizes.size() != np || sendOffsets.size() != np)
    )
    {
        FatalErrorInFunction
            << "Size of sendSizes " << sendSizes.size()
            << " or sendOffsets " << sendOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvSize);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Scatterv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                MPI_BYTE,
                recvData,
                recvSize,
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatterv failed for sendSizes " << sendSizes
                << " sendOffsets " << sendOffsets
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addScatterTime();
    }
}


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
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


int Foam::UPstream::allocateTag(const word& s)
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
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


void Foam::UPstream::freeTag(const char* s, const int tag)
{
    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


void Foam::UPstream::freeTag(const word& s, const int tag)
{
    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


// ************************************************************************* //
