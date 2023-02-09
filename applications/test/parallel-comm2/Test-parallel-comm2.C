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
#include <mpi.h>

using namespace Foam;

void rankInfo(const label comm)
{
    const int ranki = UPstream::myProcNo(comm);

    Pout<< "comm:" << comm
        << "(parent:" << UPstream::parent(comm) << ')'
        << " rank:" << ranki
        << "(sub:" << UPstream::is_subrank(comm)
        << ") nProcs:" << UPstream::nProcs(comm)
        << " baseProcNo:" << UPstream::baseProcNo(comm, ranki);
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
        << "nProcs = " << UPstream::nProcs()
        << " with " << UPstream::nComms() << " predefined comm(s)" << nl;

    Info<< nl;

    //- Process IDs within a given communicator
    Info<< "procIDs: "
        << flatOutput(UPstream::procID(UPstream::worldComm)) << endl;

    rankInfo(UPstream::worldComm);
    Pout<< endl;

    const int myProci = UPstream::myProcNo(UPstream::worldComm);
    int localRanki = myProci;

    labelList subRanks;
    UPstream::communicator newComm;

    #if 1
    // With first ranks
    subRanks = identity(UPstream::nProcs(UPstream::worldComm) / 2);

    newComm.reset(UPstream::worldComm, subRanks);
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
    subRanks = identity(UPstream::nProcs(UPstream::worldComm));

    for (label& val : subRanks)
    {
        if (val % 2) val = -1;
    }

    newComm.reset(UPstream::worldComm, subRanks);
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

    if (Pstream::parRun() && args.found("comm-split"))
    {
        MPI_Comm hostComm;
        MPI_Comm_split_type
        (
            MPI_COMM_WORLD,
            MPI_COMM_TYPE_SHARED,  // OMPI_COMM_TYPE_NODE
            0, MPI_INFO_NULL, &hostComm
        );

        int host_nprocs = 0;
        int host_rank = 0;
        MPI_Comm_size(hostComm, &host_nprocs);
        MPI_Comm_rank(hostComm, &host_rank);

        Pout<< nl << "Host comm with "
            << host_rank << " / " << host_nprocs
            << " (using MPI_Comm_split_type)" << endl;

        MPI_Comm_free(&hostComm);
    }
    if (Pstream::parRun() && args.found("host-comm"))
    {
        // Host communicator, based on the current worldComm
        // Use hostname
        // Lowest rank per hostname is the IO rank

        label numprocs = UPstream::nProcs(UPstream::globalComm);

        stringList hosts(numprocs);
        hosts[Pstream::myProcNo(UPstream::globalComm)] = hostName();

        labelList hostIDs_;

        // Compact
        if (Pstream::master(UPstream::globalComm))
        {
            DynamicList<word> hostNames(numprocs);
            hostIDs_.resize_nocopy(numprocs);

            forAll(hosts, proci)
            {
                const word& host = hosts[proci];

                hostIDs_[proci] = hostNames.find(host);

                if (hostIDs_[proci] == -1)
                {
                    hostIDs_[proci] = hostNames.size();
                    hostNames.push_back(host);
                }
            }
        }

        Pstream::broadcasts(UPstream::globalComm, hostIDs_);

        const label myHostId =
            hostIDs_[Pstream::myProcNo(UPstream::globalComm)];

        DynamicList<label> subRanks;
        forAll(hostIDs_, proci)
        {
            if (hostIDs_[proci] == myHostId)
            {
                subRanks.push_back(proci);
            }
        }

        // Allocate new communicator with globalComm as its parent
        const label hostComm =
            UPstream::allocateCommunicator
            (
                UPstream::globalComm,  // parent
                subRanks,
                true
            );

        Pout<< nl << "Host comm with "
            << UPstream::myProcNo(hostComm)
            << " / " << UPstream::nProcs(hostComm)
            << nl;

        UPstream::freeCommunicator(hostComm, true);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
