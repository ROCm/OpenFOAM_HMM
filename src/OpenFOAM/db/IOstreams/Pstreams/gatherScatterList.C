/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a list with element procID the data from processor
    procID. Before calling every processor should insert its value into
    Values[UPstream::myProcNo()].
    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use scatter to distribute the data.

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
void Pstream::gatherList
(
    const List<UPstream::commsStruct>& comms,
    List<T>& Values
)
{
    if (UPstream::parRun())
    {
        if (Values.size() != UPstream::nProcs())
        {
            FatalErrorIn
            (
                "UPstream::gatherList(const List<UPstream::commsStruct>&"
                ", List<T>)"
            )   << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << UPstream::nProcs()
                << Foam::abort(FatalError);
        }

        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo()];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];
            const labelList& belowLeaves = comms[belowID].allBelow();

            if (contiguous<T>())
            {
                List<T> receivedValues(belowLeaves.size() + 1);

                UIPstream::read
                (
                    UPstream::scheduled,
                    belowID,
                    reinterpret_cast<char*>(receivedValues.begin()),
                    receivedValues.byteSize()
                );

                Values[belowID] = receivedValues[0];

                forAll(belowLeaves, leafI)
                {
                    Values[belowLeaves[leafI]] = receivedValues[leafI + 1];
                }
            }
            else
            {
                IPstream fromBelow(UPstream::scheduled, belowID);
                fromBelow >> Values[belowID];

                if (debug & 2)
                {
                    Pout<< " received through "
                        << belowID << " data from:" << belowID
                        << " data:" << Values[belowID] << endl;
                }

                // Receive from all other processors below belowID
                forAll(belowLeaves, leafI)
                {
                    label leafID = belowLeaves[leafI];
                    fromBelow >> Values[leafID];

                    if (debug & 2)
                    {
                        Pout<< " received through "
                            << belowID << " data from:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }

        // Send up from Values:
        // - my own value first
        // - all belowLeaves next
        if (myComm.above() != -1)
        {
            const labelList& belowLeaves = myComm.allBelow();

            if (debug & 2)
            {
                Pout<< " sending to " << myComm.above()
                    << " data from me:" << UPstream::myProcNo()
                    << " data:" << Values[UPstream::myProcNo()] << endl;
            }

            if (contiguous<T>())
            {
                List<T> sendingValues(belowLeaves.size() + 1);
                sendingValues[0] = Values[UPstream::myProcNo()];

                forAll(belowLeaves, leafI)
                {
                    sendingValues[leafI + 1] = Values[belowLeaves[leafI]];
                }

                OPstream::write
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(sendingValues.begin()),
                    sendingValues.byteSize()
                );
            }
            else
            {
                OPstream toAbove(UPstream::scheduled, myComm.above());
                toAbove << Values[UPstream::myProcNo()];

                forAll(belowLeaves, leafI)
                {
                    label leafID = belowLeaves[leafI];

                    if (debug & 2)
                    {
                        Pout<< " sending to "
                            << myComm.above() << " data from:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                    toAbove << Values[leafID];
                }
            }
        }
    }
}


template <class T>
void Pstream::gatherList(List<T>& Values)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        gatherList(UPstream::linearCommunication(), Values);
    }
    else
    {
        gatherList(UPstream::treeCommunication(), Values);
    }
}


template <class T>
void Pstream::scatterList
(
    const List<UPstream::commsStruct>& comms,
    List<T>& Values
)
{
    if (UPstream::parRun())
    {
        if (Values.size() != UPstream::nProcs())
        {
            FatalErrorIn
            (
                "UPstream::scatterList(const List<UPstream::commsStruct>&"
                ", List<T>)"
            )   << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << UPstream::nProcs()
                << Foam::abort(FatalError);
        }

        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            const labelList& notBelowLeaves = myComm.allNotBelow();

            if (contiguous<T>())
            {
                List<T> receivedValues(notBelowLeaves.size());

                UIPstream::read
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(receivedValues.begin()),
                    receivedValues.byteSize()
                );

                forAll(notBelowLeaves, leafI)
                {
                    Values[notBelowLeaves[leafI]] = receivedValues[leafI];
                }
            }
            else
            {
                IPstream fromAbove(UPstream::scheduled, myComm.above());

                forAll(notBelowLeaves, leafI)
                {
                    label leafID = notBelowLeaves[leafI];
                    fromAbove >> Values[leafID];

                    if (debug)
                    {
                        Pout<< " received through "
                            << myComm.above() << " data for:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];
            const labelList& notBelowLeaves = comms[belowID].allNotBelow();

            if (contiguous<T>())
            {
                List<T> sendingValues(notBelowLeaves.size());

                forAll(notBelowLeaves, leafI)
                {
                    sendingValues[leafI] = Values[notBelowLeaves[leafI]];
                }

                OPstream::write
                (
                    UPstream::scheduled,
                    belowID,
                    reinterpret_cast<const char*>(sendingValues.begin()),
                    sendingValues.byteSize()
                );
            }
            else
            {
                OPstream toBelow(UPstream::scheduled, belowID);

                // Send data destined for all other processors below belowID
                forAll(notBelowLeaves, leafI)
                {
                    label leafID = notBelowLeaves[leafI];
                    toBelow << Values[leafID];

                    if (debug)
                    {
                        Pout<< " sent through "
                            << belowID << " data for:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }
    }
}


template <class T>
void Pstream::scatterList(List<T>& Values)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        scatterList(UPstream::linearCommunication(), Values);
    }
    else
    {
        scatterList(UPstream::treeCommunication(), Values);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
