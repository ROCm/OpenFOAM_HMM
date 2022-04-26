/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

InClass
    Foam::processorTopology

Description
    Define the processor-processor connection table by walking a list of
    patches and detecting the processor ones.

    This has been split off as a separate include from processorTopology
    to minimize dependencies.

Warning
    Does not currently correctly support multiple processor
    patches connecting two processors.

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorTopologyNew_H
#define Foam_processorTopologyNew_H

#include "processorTopology.H"
#include "bitSet.H"
#include "commSchedule.H"
#include "DynamicList.H"
#include "Pstream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class ProcPatch, class PatchListType>
Foam::processorTopology Foam::processorTopology::New
(
    const PatchListType& patches,
    const label comm
)
{
    const label myProci = Pstream::myProcNo(comm);
    const label nProcs = Pstream::nProcs(comm);

    processorTopology procTopo;

    auto& procToProcs = procTopo.procNeighbours_;
    auto& procToPatch = procTopo.procPatchMap_;
    auto& schedule = procTopo.patchSchedule_;

    procToProcs.resize(nProcs);
    schedule.resize(2*patches.size());


    if (Pstream::parRun())
    {
        // Fill my 'slot' with my neighbours
        auto& procSlot = procToProcs[myProci];

        bitSet isNeighbour(procToProcs.size());

        forAll(patches, patchi)
        {
            const auto* cpp = isA<ProcPatch>(patches[patchi]);
            if (cpp)
            {
                const label nbrProci = cpp->neighbProcNo();

                isNeighbour.set(nbrProci);
            }
        }

        // The neighbouring procs in sorted (ascending) order
        procSlot = isNeighbour.sortedToc();

        const label maxNbrProci = procSlot.empty() ? -1 : procSlot.last();

        // Note could also use Map<label> if desired
        procToPatch.resize_nocopy(maxNbrProci + 1);
        procToPatch = -1;

        forAll(patches, patchi)
        {
            const auto* cpp = isA<ProcPatch>(patches[patchi]);
            if (cpp)
            {
                const label nbrProci = cpp->neighbProcNo();

                // Reverse map
                procToPatch[nbrProci] = patchi;
            }
        }

        // Synchronize on all processors
        Pstream::allGatherList(procToProcs, UPstream::msgType(), comm);
    }


    // Define the schedule

    label patchEvali = 0;

    // 1. All non-processor patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(patches, patchi)
    {
        if (!isA<ProcPatch>(patches[patchi]))
        {
            schedule[patchEvali++].setInitEvaluate(patchi);
            schedule[patchEvali++].setEvaluate(patchi);
        }
    }


    // 2. All processor patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    if
    (
        Pstream::parRun()
     && Pstream::defaultCommsType == Pstream::commsTypes::scheduled
    )
    {
        // Determine the schedule for all processor patches.
        // Insert processor pair once to determine the schedule.
        // Each processor pair stands for both send and receive.

        label nComms = 0;
        for (const labelList& nbrProcs : procToProcs)
        {
            nComms += nbrProcs.size();
        }
        DynamicList<labelPair> comms(nComms/2);

        forAll(procToProcs, proci)
        {
            for (const label nbrProci : procToProcs[proci])
            {
                if (proci < nbrProci)
                {
                    comms.append(labelPair(proci, nbrProci));
                }
            }
        }

        // Determine a schedule.
        labelList mySchedule
        (
            commSchedule
            (
                nProcs,
                comms
            ).procSchedule()[myProci]
        );

        for (const label scheduleIndex : mySchedule)
        {
            // Get the other processor
            label nbrProci = comms[scheduleIndex][0];
            if (nbrProci == myProci)
            {
                nbrProci = comms[scheduleIndex][1];
            }
            const label patchi = procToPatch[nbrProci];

            if (myProci > nbrProci)
            {
                schedule[patchEvali++].setInitEvaluate(patchi);
                schedule[patchEvali++].setEvaluate(patchi);
            }
            else
            {
                schedule[patchEvali++].setEvaluate(patchi);
                schedule[patchEvali++].setInitEvaluate(patchi);
            }
        }
    }
    else
    {
        // Non-blocking schedule for processor patches

        // initEvaluate
        forAll(patches, patchi)
        {
            if (isA<ProcPatch>(patches[patchi]))
            {
                schedule[patchEvali++].setInitEvaluate(patchi);
            }
        }

        // evaluate
        forAll(patches, patchi)
        {
            if (isA<ProcPatch>(patches[patchi]))
            {
                schedule[patchEvali++].setEvaluate(patchi);
            }
        }
    }

    return procTopo;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
