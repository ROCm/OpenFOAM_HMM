/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Variant of gather, scatter.
    Normal gather uses:
    - default construct and read (>>) from Istream
    - binary operator and assignment operator to combine values

    combineGather uses:
    - construct from Istream
    - modify operator which modifies its lhs

\*---------------------------------------------------------------------------*/

#include "OPstream.H"
#include "IPstream.H"
#include "IOstreams.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class CombineOp>
void Foam::Pstream::combineGather
(
    const List<UPstream::commsStruct>& comms,
    T& value,
    const CombineOp& cop,
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
            if (is_contiguous<T>::value)
            {
                T received;

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    reinterpret_cast<char*>(&received),
                    sizeof(T),
                    tag,
                    comm
                );

                if (debug & 2)
                {
                    Pout<< " received from "
                        << belowID << " data:" << received << endl;
                }

                cop(value, received);
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
                T received(fromBelow);

                if (debug & 2)
                {
                    Pout<< " received from "
                        << belowID << " data:" << received << endl;
                }

                cop(value, received);
            }
        }

        // Send up value
        if (myComm.above() != -1)
        {
            if (debug & 2)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << value << endl;
            }

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
void Foam::Pstream::combineScatter
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
        const UPstream::commsStruct& myComm = comms[UPstream::myProcNo(comm)];

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
                value = T(fromAbove);
            }
        }

        // Send to my downstairs neighbours
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


template<class T, class CombineOp>
void Foam::Pstream::combineGather
(
    T& value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    Pstream::combineGather
    (
        UPstream::whichCommunication(comm),
        value,
        cop,
        tag,
        comm
    );
}


template<class T>
void Foam::Pstream::combineScatter
(
    T& value,
    const int tag,
    const label comm
)
{
    #ifndef Foam_Pstream_scatter_nobroadcast
    Pstream::broadcast(value, comm);
    #else
    Pstream::combineScatter
    (
        UPstream::whichCommunication(comm),
        value,
        tag,
        comm
    );
    #endif
}


template<class T, class CombineOp>
void Foam::Pstream::combineReduce
(
    const List<UPstream::commsStruct>& comms,
    T& value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    Pstream::combineGather(comms, value, cop, tag, comm);
    Pstream::broadcast(value, comm);
}


template<class T, class CombineOp>
void Foam::Pstream::combineReduce
(
    T& value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        const auto& comms = UPstream::whichCommunication(comm);

        Pstream::combineGather(comms, value, cop, tag, comm);
        Pstream::broadcast(value, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class CombineOp>
void Foam::Pstream::listCombineGather
(
    const List<UPstream::commsStruct>& comms,
    List<T>& values,
    const CombineOp& cop,
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
            if (is_contiguous<T>::value)
            {
                List<T> received(values.size());

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    received.data_bytes(),
                    received.size_bytes(),
                    tag,
                    comm
                );

                if (debug & 2)
                {
                    Pout<< " received from "
                        << belowID << " data:" << received << endl;
                }

                forAll(values, i)
                {
                    cop(values[i], received[i]);
                }
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
                List<T> received(fromBelow);

                if (debug & 2)
                {
                    Pout<< " received from "
                        << belowID << " data:" << received << endl;
                }

                forAll(values, i)
                {
                    cop(values[i], received[i]);
                }
            }
        }

        // Send up values
        if (myComm.above() != -1)
        {
            if (debug & 2)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << values << endl;
            }

            if (is_contiguous<T>::value)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    values.cdata_bytes(),
                    values.size_bytes(),
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
                toAbove << values;
            }
        }
    }
}


template<class T>
void Foam::Pstream::listCombineScatter
(
    const List<UPstream::commsStruct>& comms,
    List<T>& values,
    const int tag,
    const label comm
)
{
    #ifndef Foam_Pstream_scatter_nobroadcast
    Pstream::broadcast(values, comm);
    #else
    if (UPstream::is_parallel(comm))
    {
        // My communication order
        const UPstream::commsStruct& myComm = comms[UPstream::myProcNo(comm)];

        // Receive from up
        if (myComm.above() != -1)
        {
            if (is_contiguous<T>::value)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    values.data_bytes(),
                    values.size_bytes(),
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
                fromAbove >> values;
            }
        }

        // Send to my downstairs neighbours
        forAllReverse(myComm.below(), belowI)
        {
            const label belowID = myComm.below()[belowI];

            if (is_contiguous<T>::value)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    values.cdata_bytes(),
                    values.size_bytes(),
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
                toBelow << values;
            }
        }
    }
    #endif
}


template<class T, class CombineOp>
void Foam::Pstream::listCombineGather
(
    List<T>& values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    Pstream::listCombineGather
    (
        UPstream::whichCommunication(comm),
        values,
        cop,
        tag,
        comm
    );
}


template<class T>
void Foam::Pstream::listCombineScatter
(
    List<T>& values,
    const int tag,
    const label comm
)
{
    #ifndef Foam_Pstream_scatter_nobroadcast
    Pstream::broadcast(values, comm);
    #else
    Pstream::listCombineScatter
    (
        UPstream::whichCommunication(comm),
        values,
        tag,
        comm
    );
    #endif
}


template<class T, class CombineOp>
void Foam::Pstream::listCombineReduce
(
    List<T>& values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        const auto& comms = UPstream::whichCommunication(comm);

        Pstream::listCombineGather(comms, values, cop, tag, comm);
        Pstream::broadcast(values, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Container, class CombineOp>
void Foam::Pstream::mapCombineGather
(
    const List<UPstream::commsStruct>& comms,
    Container& values,
    const CombineOp& cop,
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
            // Map/HashTable: non-contiguous

            IPstream fromBelow
            (
                UPstream::commsTypes::scheduled,
                belowID,
                0,
                tag,
                comm
            );
            Container received(fromBelow);

            if (debug & 2)
            {
                Pout<< " received from "
                    << belowID << " data:" << received << endl;
            }

            for
            (
                auto recvIter = received.cbegin();
                recvIter != received.cend();
                ++recvIter
            )
            {
                auto masterIter = values.find(recvIter.key());

                if (masterIter.good())
                {
                    // Combine with existing
                    cop(masterIter.val(), recvIter.val());
                }
                else
                {
                    // Insert new key/value
                    values.insert(recvIter.key(), recvIter.val());
                }
            }
        }

        // Send up values
        if (myComm.above() != -1)
        {
            if (debug & 2)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << values << endl;
            }

            OPstream toAbove
            (
                UPstream::commsTypes::scheduled,
                myComm.above(),
                0,
                tag,
                comm
            );
            toAbove << values;
        }
    }
}


template<class Container>
void Foam::Pstream::mapCombineScatter
(
    const List<UPstream::commsStruct>& comms,
    Container& values,
    const int tag,
    const label comm
)
{
    #ifndef Foam_Pstream_scatter_nobroadcast
    Pstream::broadcast(values, comm);
    #else
    if (UPstream::is_parallel(comm))
    {
        // My communication order
        const UPstream::commsStruct& myComm = comms[UPstream::myProcNo(comm)];

        // Receive from up
        if (myComm.above() != -1)
        {
            IPstream fromAbove
            (
                UPstream::commsTypes::scheduled,
                myComm.above(),
                0,
                tag,
                comm
            );
            fromAbove >> values;

            if (debug & 2)
            {
                Pout<< " received from "
                    << myComm.above() << " data:" << values << endl;
            }
        }

        // Send to my downstairs neighbours
        forAllReverse(myComm.below(), belowI)
        {
            const label belowID = myComm.below()[belowI];

            if (debug & 2)
            {
                Pout<< " sending to " << belowID << " data:" << values << endl;
            }

            OPstream toBelow
            (
                UPstream::commsTypes::scheduled,
                belowID,
                0,
                tag,
                comm
            );
            toBelow << values;
        }
    }
    #endif
}


template<class Container, class CombineOp>
void Foam::Pstream::mapCombineGather
(
    Container& values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    Pstream::mapCombineGather
    (
        UPstream::whichCommunication(comm),
        values,
        cop,
        tag,
        comm
    );
}


template<class Container>
void Foam::Pstream::mapCombineScatter
(
    Container& values,
    const int tag,
    const label comm
)
{
    #ifndef Foam_Pstream_scatter_nobroadcast
    Pstream::broadcast(values, comm);
    #else
    Pstream::mapCombineScatter
    (
        UPstream::whichCommunication(comm),
        values,
        tag,
        comm
    );
    #endif
}


template<class Container, class CombineOp>
void Foam::Pstream::mapCombineReduce
(
    Container& values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        const auto& comms = UPstream::whichCommunication(comm);

        Pstream::mapCombineGather(comms, values, cop, tag, comm);
        Pstream::broadcast(values, comm);
    }
}


// ************************************************************************* //
