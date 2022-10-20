/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "masterCoarsestGAMGProcAgglomeration.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(masterCoarsestGAMGProcAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        masterCoarsestGAMGProcAgglomeration,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterCoarsestGAMGProcAgglomeration::masterCoarsestGAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    nProcessorsPerMaster_
    (
        controlDict.getOrDefault<label>
        (
            "nProcessorsPerMaster",
            0,
            keyType::LITERAL
        )
    )
{
    const auto* ePtr = controlDict.findEntry("nMasters", keyType::LITERAL);
    if (ePtr)
    {
        if (nProcessorsPerMaster_ > 0)
        {
            FatalIOErrorInFunction(controlDict)
                << "Cannot specify both \"nMasters\" and"
                << " \"nProcessorsPerMaster\"" << exit(FatalIOError);
        }

        const label nMasters(readLabel(ePtr->stream()));

        if (nMasters <= 0)
        {
            FatalIOErrorInFunction(controlDict)
                << "Illegal value \"nMasters\" "
                << nMasters << exit(FatalIOError);
        }

        nProcessorsPerMaster_ =
            (Pstream::nProcs(agglom.mesh().comm())+nMasters-1)
          / nMasters;
    }

    if (nProcessorsPerMaster_ < 0)
    {
        FatalIOErrorInFunction(controlDict)
            << "Illegal value \"nProcessorsPerMaster\" "
            << nProcessorsPerMaster_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterCoarsestGAMGProcAgglomeration::
~masterCoarsestGAMGProcAgglomeration()
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

bool Foam::masterCoarsestGAMGProcAgglomeration::agglomerate()
{
    if (debug & 2)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        // Agglomerate one but last level (since also agglomerating
        // restrictAddressing)
        label fineLevelIndex = agglom_.size()-1;

        if (agglom_.hasMeshLevel(fineLevelIndex))
        {
            // Get the fine mesh
            const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex);
            label levelComm = levelMesh.comm();
            label nProcs = UPstream::nProcs(levelComm);

            if (nProcs > 1)
            {
                // Processor restriction map: per processor the coarse processor
                labelList procAgglomMap(nProcs);

                if (nProcessorsPerMaster_ > 0)
                {
                    forAll(procAgglomMap, fineProci)
                    {
                        procAgglomMap[fineProci] =
                        (
                            fineProci
                          / nProcessorsPerMaster_
                        );
                    }
                }
                else
                {
                    procAgglomMap = Zero;
                }

                // Master processor
                labelList masterProcs;
                // Local processors that agglomerate. agglomProcIDs[0] is in
                // masterProc.
                List<label> agglomProcIDs;
                GAMGAgglomeration::calculateRegionMaster
                (
                    levelComm,
                    procAgglomMap,
                    masterProcs,
                    agglomProcIDs
                );

                if (debug)
                {
                    if (masterProcs.size())
                    {
                        labelListList masterToProcs
                        (
                            invertOneToMany
                            (
                                masterProcs.size(),
                                procAgglomMap
                            )
                        );
                        Info<< typeName << " : agglomerating" << nl
                            << "\tmaster\tnProcs\tprocIDs" << endl;
                        for (const auto& p : masterToProcs)
                        {
                            Info<< '\t' << p[0]
                                << '\t' << p.size()
                                << '\t'
                                << flatOutput(SubList<label>(p, p.size()-1, 1))
                                << endl;
                        }
                    }
                }


                // Allocate a communicator for the processor-agglomerated matrix
                comms_.append
                (
                    UPstream::allocateCommunicator
                    (
                        levelComm,
                        masterProcs
                    )
                );

                // Use processor agglomeration maps to do the actual collecting.
                if (Pstream::myProcNo(levelComm) != -1)
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
    }

    // Print a bit
    if (debug & 2)
    {
        Pout<< nl << "Agglomerated mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    return true;
}


// ************************************************************************* //
