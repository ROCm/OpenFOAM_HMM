/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "hostCollatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "bitSet.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(hostCollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostCollatedFileOperation,
        word
    );

    // Register initialisation routine. Signals need for threaded mpi and
    // handles command line arguments
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        hostCollatedFileOperationInitialise,
        word,
        hostCollated
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fileOperations::hostCollatedFileOperation::subRanks
(
    const label n
)
{
    labelList mainIOranks(fileOperation::getGlobalIORanks());

    if (mainIOranks.empty())
    {
        mainIOranks = fileOperation::getGlobalHostIORanks();
    }

    labelRange subRange = fileOperation::subRanks(mainIOranks);

    return identity(subRange);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::fileOperations::hostCollatedFileOperation::init(bool verbose)
{
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (verbose)
    {
        this->printBanner(ioRanks_.size());
    }
}


Foam::fileOperations::hostCollatedFileOperation::hostCollatedFileOperation
(
    bool verbose
)
:
    collatedFileOperation
    (
        UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            subRanks(UPstream::nProcs())
        ),
        (UPstream::parRun() ? labelList() : getGlobalIORanks()),
        false  // verbose
    ),
    managedComm_(getManagedComm(comm_))  // Possibly locally allocated
{
    init(verbose);
}


void Foam::fileOperations::hostCollatedFileOperation::storeComm() const
{
    // From externally -> locally managed
    managedComm_ = getManagedComm(comm_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::hostCollatedFileOperation::~hostCollatedFileOperation()
{
    UPstream::freeCommunicator(managedComm_);
}


// ************************************************************************* //
