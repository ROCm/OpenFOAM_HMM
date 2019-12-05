/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Test-parallel-external-init

Description
    Simulate starting MPI outside of OpenFOAM

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Pstream.H"

#include <mpi.h>
#include <iostream>

using namespace Foam;


bool startMPI()
{
    int nprocs = 0, rank = 0;

    MPI_Init(nullptr, nullptr);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (nprocs && rank == 0)
    {
        std::cout<< nl << "Using MPI with " << nprocs << " procs" << nl << nl;
    }

    return true;
}


bool stopMPI()
{
    Info<< nl << "Stopping MPI" << nl << nl;

    MPI_Finalize();

    return true;
}


string message()
{
    return
    (
        "rank " + name(Pstream::myProcNo())
      + " / " + name(Pstream::nProcs()) + "\n"
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    UPstream::debug = 1;

    startMPI();

    #include "setRootCase.H"

    Pout<< message().c_str();

    stopMPI();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
