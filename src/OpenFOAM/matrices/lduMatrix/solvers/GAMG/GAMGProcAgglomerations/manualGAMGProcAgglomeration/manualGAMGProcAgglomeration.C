/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "manualGAMGProcAgglomeration.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manualGAMGProcAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        manualGAMGProcAgglomeration,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manualGAMGProcAgglomeration::manualGAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    procAgglomMaps_(controlDict.lookup("processorAgglomeration"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::manualGAMGProcAgglomeration::
~manualGAMGProcAgglomeration()
{
    forAllReverse(comms_, i)
    {
        if (comms_[i] != -1)
        {
            UPstream::freeCommunicator(comms_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::manualGAMGProcAgglomeration::agglomerate()
{
    if (debug)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        forAll(procAgglomMaps_, i)
        {
            const label fineLevelIndex = procAgglomMaps_[i].first();

            if (fineLevelIndex >= agglom_.size())
            {
                WarningIn("manualGAMGProcAgglomeration::agglomerate()")
                    << "Ignoring specification for level " << fineLevelIndex
                    << " since outside agglomeration." << endl;

                continue;
            }


            if (agglom_.hasMeshLevel(fineLevelIndex))
            {
                const labelList& procAgglomMap = procAgglomMaps_[i].second();
                // Get the fine mesh
                const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex);


                if (procAgglomMap.size() != Pstream::nProcs(levelMesh.comm()))
                {
                    FatalErrorIn("manualGAMGProcAgglomeration::agglomerate()")
                        << "At level " << fineLevelIndex
                        << "The fine-to-coarse agglomeration map size "
                        << procAgglomMap.size() << " differs from the"
                        << " number of processors "
                        << Pstream::nProcs(levelMesh.comm())
                        << exit(FatalError);
                }

                // Master processor
                labelList masterProcs;
                // Local processors that agglomerate. agglomProcIDs[0] is in
                // masterProc.
                List<int> agglomProcIDs;
                GAMGAgglomeration::calculateRegionMaster
                (
                    levelMesh.comm(),
                    procAgglomMap,
                    masterProcs,
                    agglomProcIDs
                );

                // Allocate a communicator for the processor-agglomerated matrix
                comms_.append
                (
                    UPstream::allocateCommunicator
                    (
                        levelMesh.comm(),
                        masterProcs
                    )
                );

                // Use procesor agglomeration maps to do the actual collecting.
                if (Pstream::myProcNo(levelMesh.comm()) != -1)
                {
                    GAMGProcAgglomeration::agglomerate
                    (
                        fineLevelIndex,
                        procAgglomMap,
                        masterProcs,
                        agglomProcIDs,
                        comms_.last()
                    );
                }
            }
        }

        // Print a bit
        if (debug)
        {
            Pout<< nl << "Agglomerated mesh overview" << endl;
            printStats(Pout, agglom_);
        }
    }

    return true;
}


// ************************************************************************* //
