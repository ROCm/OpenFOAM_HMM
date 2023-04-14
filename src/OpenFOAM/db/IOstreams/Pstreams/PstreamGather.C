/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a single value constructed from the values
    on individual processors using a user-specified operator.

\*---------------------------------------------------------------------------*/

#include "OPstream.H"
#include "IPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BinaryOp>
void Foam::Pstream::gather
(
    const List<UPstream::commsStruct>& comms,
    T& value,
    const BinaryOp& bop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        // My communication order
        const commsStruct& myComm = comms[UPstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        for (const label belowID : myComm.below())
        {
            T received;

            if (is_contiguous<T>::value)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    reinterpret_cast<char*>(&received),
                    sizeof(T),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream fromBelow
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    0,
                    tag,
                    comm
                );
                fromBelow >> received;
            }

            value = bop(value, received);
        }

        // Send up value
        if (myComm.above() != -1)
        {
            if (is_contiguous<T>::value)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(&value),
                    sizeof(T),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toAbove
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    0,
                    tag,
                    comm
                );
                toAbove << value;
            }
        }
    }
}


template<class T>
void Foam::Pstream::scatter
(
    const List<UPstream::commsStruct>& comms,
    T& value,
    const int tag,
    const label comm
)
{
    #ifndef Foam_Pstream_scatter_nobroadcast
    Pstream::broadcast(value, comm);
    #else
    if (UPstream::is_parallel(comm))
    {
        // My communication order
        const commsStruct& myComm = comms[UPstream::myProcNo(comm)];

        // Receive from up
        if (myComm.above() != -1)
        {
            if (is_contiguous<T>::value)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(&value),
                    sizeof(T),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream fromAbove
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    0,
                    tag,
                    comm
                );
                fromAbove >> value;
            }
        }

        // Send to my downstairs neighbours. Note reverse order (compared to
        // receiving). This is to make sure to send to the critical path
        // (only when using a tree schedule!) first.
        forAllReverse(myComm.below(), belowI)
        {
            const label belowID = myComm.below()[belowI];

            if (is_contiguous<T>::value)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    reinterpret_cast<const char*>(&value),
                    sizeof(T),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toBelow
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    0,
                    tag,
                    comm
                );
                toBelow << value;
            }
        }
    }
    #endif
}


template<class T, class BinaryOp>
void Foam::Pstream::gather
(
    T& value,
    const BinaryOp& bop,
    const int tag,
    const label comm
)
{
    Pstream::gather(UPstream::whichCommunication(comm), value, bop, tag, comm);
}


template<class T>
void Foam::Pstream::scatter(T& value, const int tag, const label comm)
{
    #ifndef Foam_Pstream_scatter_nobroadcast
    Pstream::broadcast(value, comm);
    #else
    Pstream::scatter(UPstream::whichCommunication(comm), value, tag, comm);
    #endif
}


// ************************************************************************* //
