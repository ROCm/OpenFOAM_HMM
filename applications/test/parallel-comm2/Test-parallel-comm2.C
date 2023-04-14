/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

Application
    Test-parallel-comm2

Description
    Basic communicator tests

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "Pair.H"
#include "Tuple2.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"
#include "SHA1.H"
#include <mpi.h>

using namespace Foam;

void rankInfo(const label comm)
{
    const int ranki = UPstream::myProcNo(comm);

    Pout<< "comm:" << comm
        << "(parent:" << UPstream::parent(comm) << ')'
        << " rank:" << ranki
        << "(sub:" << UPstream::is_subrank(comm)
        << ") nProcs:" << UPstream::nProcs(comm);
      // << " baseProcNo:" << UPstream::baseProcNo(comm, ranki);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();
    argList::addBoolOption("verbose", "Set debug level");
    argList::addBoolOption("comm-split", "Test simple comm split");
    argList::addBoolOption("host-comm", "Test DIY host-comm split");

    // Capture manually. We need values before proper startup
    int nVerbose = 0;
    for (int argi = 1; argi < argc; ++argi)
    {
        if (strcmp(argv[argi], "-verbose") == 0)
        {
            ++nVerbose;
        }
    }

    UPstream::debug = nVerbose;

    #include "setRootCase.H"

    Info<< nl
        << "parallel:" << UPstream::parRun()
        << "nProcs = " << UPstream::nProcs()
        << " with " << UPstream::nComms() << " predefined comm(s)."
        << " proc:" << UPstream::myProcNo() << nl;

    Info<< nl;

    //- Process IDs within a given communicator
    Info<< "procIDs: "
        << flatOutput(UPstream::procID(UPstream::commWorld())) << endl;

    rankInfo(UPstream::commWorld());
    Pout<< endl;

    const int myProci = UPstream::myProcNo(UPstream::commWorld());
    int localRanki = myProci;

    labelList subRanks;
    UPstream::communicator newComm;

    if (!args.found("comm-split") && !args.found("host-comm"))
    {
        #if 1
        // With first ranks
        subRanks = identity(UPstream::nProcs(UPstream::commWorld()) / 2);

        newComm.reset(UPstream::commWorld(), subRanks);
        localRanki = UPstream::myProcNo(newComm);

        Pout.prefix() =
        (
            '[' + Foam::name(myProci) + " a:" + Foam::name(localRanki) + "] "
        );

        Info<< "procIDs: "
            << flatOutput(UPstream::procID(newComm)) << endl;

        rankInfo(newComm);
        Pout<< endl;
        #endif

        #if 1
        // With every other rank
        subRanks = identity(UPstream::nProcs(UPstream::commWorld()));

        for (label& val : subRanks)
        {
            if (val % 2) val = -1;
        }

        newComm.reset(UPstream::commWorld(), subRanks);
        localRanki = UPstream::myProcNo(newComm);

        Pout.prefix() =
        (
            '[' + Foam::name(myProci) + " b:" + Foam::name(localRanki) + "] "
        );

        Info<< "procIDs: "
            << flatOutput(UPstream::procID(newComm)) << endl;

        rankInfo(newComm);
        Pout<< endl;
        #endif
    }

    if (Pstream::parRun() && args.found("comm-split"))
    {
        int world_nprocs = 0;
        int world_rank = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &world_nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        int host_nprocs = 0;
        int host_rank = -1;
        MPI_Comm hostComm;
        MPI_Comm_split_type
        (
            MPI_COMM_WORLD,
            MPI_COMM_TYPE_SHARED,  // OMPI_COMM_TYPE_NODE
            0, MPI_INFO_NULL, &hostComm
        );

        MPI_Comm_size(hostComm, &host_nprocs);
        MPI_Comm_rank(hostComm, &host_rank);

        int leader_nprocs = 0;
        int leader_rank = -1;
        MPI_Comm hostMasterComm;

        if (false)
        {
            // Easy enough to use MPI_Comm_split, but slightly annoying
            // that it returns MPI_COMM_NULL for unused ranks...
            MPI_Comm hostMasterComm;
            MPI_Comm_split
            (
                MPI_COMM_WORLD,
                (host_rank == 0) ? 0 : MPI_UNDEFINED,
                0, &hostMasterComm
            );

            if (hostMasterComm != MPI_COMM_NULL)
            {
                MPI_Comm_size(hostMasterComm, &leader_nprocs);
                MPI_Comm_rank(hostMasterComm, &leader_rank);
            }
        }
        else
        {
            boolList isHostLeader(world_nprocs, false);
            isHostLeader[world_rank] = (host_rank == 0);

            MPI_Allgather
            (
                // recv is also send
                MPI_IN_PLACE, 1, MPI_C_BOOL,
                isHostLeader.data(), 1, MPI_C_BOOL,
                MPI_COMM_WORLD
            );

            Pout<< "leaders: " << isHostLeader << endl;

            DynamicList<int> subRanks(isHostLeader.size());
            forAll(isHostLeader, proci)
            {
                if (isHostLeader[proci])
                {
                    subRanks.push_back(proci);
                }
            }
            // Starting from parent
            MPI_Group parent_group;
            MPI_Comm_group(MPI_COMM_WORLD, &parent_group);

            MPI_Group active_group;
            MPI_Group_incl
            (
                parent_group,
                subRanks.size(),
                subRanks.cdata(),
                &active_group
            );

            // Create new communicator for this group
            MPI_Comm_create_group
            (
                MPI_COMM_WORLD,
                active_group,
                UPstream::msgType(),
                &hostMasterComm
            );

            // Groups not needed after this...
            MPI_Group_free(&parent_group);
            MPI_Group_free(&active_group);

            MPI_Comm_size(hostMasterComm, &leader_nprocs);
            MPI_Comm_rank(hostMasterComm, &leader_rank);
        }

        Pout<< nl << "[MPI_Comm_split_type]" << nl
            << "Host comm with " << host_rank << " / " << host_nprocs
            << " on " << hostName()
            << " master:" << (host_rank == 0)
            << " leader rank: " << leader_rank
            << " / " << leader_nprocs
            << " host leader:" << (leader_rank == 0)
            << " sub-rank:" << (leader_rank > 0)
            << nl;

        if (hostMasterComm != MPI_COMM_NULL)
        {
            MPI_Comm_free(&hostMasterComm);
        }
        MPI_Comm_free(&hostComm);
    }
    if (Pstream::parRun() && args.found("host-comm"))
    {
        // Host communicator, based on the current world communicator
        // Use hostname
        // Lowest rank per hostname is the IO rank

        label numprocs = UPstream::nProcs(UPstream::commGlobal());

        // Option 1: using hostnames
        // - pro: trivial coding
        // - con: unequal lengths, more allocations and 'hops'
        stringList hosts(numprocs);
        hosts[Pstream::myProcNo(UPstream::commGlobal())] = hostName();
        Pstream::gatherList(hosts, UPstream::msgType(), UPstream::commGlobal());


        // Option 2: using SHA1 of hostnames
        // - con: uglier coding (but only needed locally!)
        // - pro: fixed digest length enables direct MPI calls
        //        can avoid Pstream::gatherList() during setup...

        SHA1Digest myHostDigest(SHA1(hostName()).digest());

        List<SHA1Digest> digests;
        if (UPstream::master(UPstream::commGlobal()))
        {
            digests.resize(numprocs);
        }

        UPstream::mpiGather
        (
            reinterpret_cast<const char*>(myHostDigest.cdata_bytes()),
            SHA1Digest::max_size(),     // Num send per proc
            digests.data_bytes(),       // Recv
            SHA1Digest::max_size(),     // Num recv per proc
            UPstream::commGlobal()
        );

        // MPI_Gather
        // (
        //     myHostDigest.cdata_bytes(), // Send
        //     SHA1Digest::max_size(),     // Num send per proc
        //     MPI_BYTE,
        //     digests.data_bytes(),       // Recv
        //     SHA1Digest::max_size(),     // Num recv per proc
        //     MPI_BYTE,
        //     0,  //  root
        //     MPI_COMM_WORLD
        // );

        Info<< "digests: " << digests << nl;

        labelList hostIDs_(numprocs);

        // Compact
        if (UPstream::master(UPstream::commGlobal()))
        {
            DynamicList<word> hostNames(numprocs);

            forAll(hosts, proci)
            {
                const word& host = hosts[proci];

                hostIDs_[proci] = hostNames.find(host);

                if (hostIDs_[proci] < 0)
                {
                    // First appearance of host (encode as leader)
                    hostIDs_[proci] = -(hostNames.size() + 1);
                    hostNames.push_back(host);
                }
            }

            DynamicList<SHA1Digest> hostDigest(numprocs);

            forAll(digests, proci)
            {
                const SHA1Digest& dig = digests[proci];

                hostIDs_[proci] = hostDigest.find(dig);

                if (hostIDs_[proci] < 0)
                {
                    // First appearance of host (encode as leader)
                    hostIDs_[proci] = -(hostDigest.size() + 1);
                    hostDigest.push_back(dig);
                }
            }
        }

        Info<< "hosts =  " << hosts << endl;
        Info<< "hostIDs_ =  " << hostIDs_ << endl;

        UPstream::broadcast
        (
            hostIDs_.data_bytes(),
            hostIDs_.size_bytes(),
            UPstream::commGlobal(),
            UPstream::masterNo()
        );

        DynamicList<label> subRanks(numprocs);

        // Ranks for world to hostMaster
        forAll(hostIDs_, proci)
        {
            // Is host leader?
            if (hostIDs_[proci] < 0)
            {
                subRanks.push_back(proci);

                // Flip back to generic host id
                hostIDs_[proci] = -(hostIDs_[proci] + 1);
            }
        }

        // From world to hostMaster
        const label hostMasterComm =
            UPstream::allocateCommunicator
            (
                UPstream::commGlobal(),
                subRanks,
                true
            );


        const label myHostId =
            hostIDs_[Pstream::myProcNo(UPstream::commGlobal())];

        // Ranks for within a host
        subRanks.clear();
        forAll(hostIDs_, proci)
        {
            if (hostIDs_[proci] == myHostId)
            {
                subRanks.push_back(proci);
            }
        }

        // The intra-host ranks
        const label hostComm =
            UPstream::allocateCommunicator
            (
                UPstream::commGlobal(),
                subRanks,
                true
            );

        Pout<< nl << "[manual split]" << nl
            << nl << "Host comm with "
            << UPstream::myProcNo(hostComm)
            << " / " << UPstream::nProcs(hostComm)
            << " on " << hostName()
            << " master:" << UPstream::master(hostComm)
            << " leader rank: " << UPstream::myProcNo(hostMasterComm)
            << " / " << UPstream::nProcs(hostMasterComm)
            << " host leader:" << UPstream::master(hostMasterComm)
            << " sub-rank:" << UPstream::is_subrank(hostMasterComm)
            << nl;

        UPstream::freeCommunicator(hostMasterComm, true);
        UPstream::freeCommunicator(hostComm, true);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
