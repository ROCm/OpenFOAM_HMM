/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Test-parallel-comm0

Description
    Very basic checks on standard communicators

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "Pair.H"
#include "Tuple2.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"

using namespace Foam;


void printInfo(const label comm)
{
    Info<< "comm:" << comm
        << " nprocs:" << UPstream::nProcs(comm)
        << " all:" << UPstream::allProcs(comm)
        << " sub:" << UPstream::subProcs(comm) << nl;


    if (UPstream::selfComm == comm)
    {
        Pout<< "self all:" << UPstream::allProcs(comm)
            << " sub:" << UPstream::subProcs(comm) << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();
    argList::addBoolOption("verbose", "Set debug level");

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

    Info<< "worldComm : ";
    printInfo(UPstream::worldComm);

    Info<< "selfComm  : ";
    printInfo(UPstream::selfComm);

    Info<< nl;

    // Reductions (using MPI intrinsics)
    {
        label val = Pstream::myProcNo(UPstream::worldComm);

        label worldVal = returnReduce
        (
            val,
            sumOp<label>(),
            Pstream::msgType(),
            UPstream::worldComm
        );

        label selfVal = returnReduce
        (
            val,
            sumOp<label>(),
            Pstream::msgType(),
            UPstream::selfComm
        );

        Pout<< "value " << val
            << " (world) reduced " << worldVal
            << " (self) reduced " << selfVal << nl;
    }

    // Reductions (not using MPI intrinsics)
    {
        Pair<label> val
        (
            Pstream::myProcNo(UPstream::worldComm),
            Pstream::myProcNo(UPstream::worldComm)
        );

        Pair<label> worldVal = val;

        Pstream::combineReduce
        (
            worldVal,
            minFirstEqOp<label>(),
            Pstream::msgType(),
            UPstream::worldComm
        );

        Pair<label> selfVal = val;

        Pstream::combineReduce
        (
            worldVal,
            minFirstEqOp<label>(),
            Pstream::msgType(),
            UPstream::selfComm
        );

        Pout<< "value " << val
            << " (world) reduced " << worldVal
            << " (self) reduced " << selfVal << nl;
    }

    Pout<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
