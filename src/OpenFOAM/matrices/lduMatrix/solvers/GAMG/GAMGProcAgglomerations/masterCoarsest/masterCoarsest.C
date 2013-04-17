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

#include "masterCoarsest.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(masterCoarsest, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        masterCoarsest,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterCoarsest::masterCoarsest
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::masterCoarsest::agglomerate()
{
    if (debug)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        // Agglomerate one but last level (since also agglomerating
        // restrictAddressing)
        label fineLevelIndex = agglom_.size()-1;

        // Get the coarse mesh
        const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex-1);
        label levelComm = levelMesh.comm();

        // Processor restriction map: per processor the coarse processor
        labelList procAgglomMap(UPstream::nProcs(levelComm));
        //{
        //    label half = (procAgglomMap.size()+1)/2;
        //    for (label i = 0; i < half; i++)
        //    {
        //        procAgglomMap[i] = 0;
        //    }
        //    for (label i = half; i < procAgglomMap.size(); i++)
        //    {
        //        procAgglomMap[i] = 1;
        //    }
        //}
        procAgglomMap = 0;

        // Master processor
        labelList masterProcs;
        // Local processors that agglomerate. agglomProcIDs[0] is in
        // masterProc.
        List<int> agglomProcIDs;
        GAMGAgglomeration::calculateRegionMaster
        (
            levelComm,
            procAgglomMap,
            masterProcs,
            agglomProcIDs
        );

        // Allocate a communicator for the processor-agglomerated matrix
        label procAgglomComm = UPstream::allocateCommunicator
        (
            levelComm,
            masterProcs
        );


        // Above should really be determined by a procesor-agglomeration
        // engine. Use procesor agglomeration maps to do the actual
        // collecting.

        if (Pstream::myProcNo(levelComm) != -1)
        {
            GAMGProcAgglomeration::agglomerate
            (
                fineLevelIndex,
                procAgglomMap,
                masterProcs,
                agglomProcIDs,
                procAgglomComm
            );
        }
    }

    // Print a bit
    if (debug)
    {
        Pout<< nl << "Agglomerated mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    return true;
}


// ************************************************************************* //
