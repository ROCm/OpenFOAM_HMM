/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "dataSchedule.H"
#include "SortableList.H"
#include "commSchedule.H"
#include "Pstream.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dataSchedule, 0);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from separate addressing
Foam::dataSchedule::dataSchedule
(
    const labelList& sendProcs,
    const labelList& receiveProcs
)
:
    sendOrder_(Pstream::nProcs()),
    receiveOrder_(Pstream::nProcs())
{
    // Determine the schedule
    // ~~~~~~~~~~~~~~~~~~~~~~

    // Get all pairs of processors. Leave out local comms.
    List<labelPair> allComms;
    {
        HashSet<labelPair, labelPair::Hash<> > commsSet(Pstream::nProcs());

        forAll(receiveProcs, i)
        {
            if (sendProcs[i] != receiveProcs[i])
            {
                commsSet.insert(labelPair(sendProcs[i], receiveProcs[i]));
            }
        }
        allComms = commsSet.toc();
    }


    // Determine my schedule.
    labelList mySchedule
    (
        commSchedule
        (
            Pstream::nProcs(),
            allComms
        ).procSchedule()[Pstream::myProcNo()]
    );

    // Processors involved in my schedule
    List<labelPair>::operator=
    (
        IndirectList<labelPair>(allComms, mySchedule)
    );

    if (debug)
    {
        Pout<< "I need to:" << endl;
        forAll(*this, i)
        {
            const labelPair& twoProcs = operator[](i);
            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (recvProc == Pstream::myProcNo())
            {
                Pout<< "    receive from " << sendProc << endl;
            }
            else
            {
                Pout<< "    send to " << recvProc << endl;
            }
        }
    }


    // Determine the addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // Per processor the indices we have to send/receive.
    List<DynamicList<label> > dynSendOrder(Pstream::nProcs());
    List<DynamicList<label> > dynReceiveOrder(Pstream::nProcs());
    // Estimate size
    forAll(dynSendOrder, procI)
    {
        dynSendOrder[procI].setSize(sendProcs.size()/Pstream::nProcs());
        dynReceiveOrder[procI].setSize(sendProcs.size()/Pstream::nProcs());
    }

    forAll(sendProcs, sampleI)
    {
        // Note that also need to include local communication (both receiveProc
        // and sendProc on local processor)

        if (Pstream::myProcNo() == sendProcs[sampleI])
        {
            // I am the sender
            dynSendOrder[receiveProcs[sampleI]].append(sampleI);
        }
        if (Pstream::myProcNo() == receiveProcs[sampleI])
        {
            // I am the receiver
            dynReceiveOrder[sendProcs[sampleI]].append(sampleI);
        }
    }

    forAll(dynSendOrder, procI)
    {
        sendOrder_[procI].transfer(dynSendOrder[procI].shrink());
        receiveOrder_[procI].transfer(dynReceiveOrder[procI].shrink());
    }
}


// ************************************************************************* //
