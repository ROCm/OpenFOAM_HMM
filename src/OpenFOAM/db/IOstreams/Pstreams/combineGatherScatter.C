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
    Variant of gather, scatter.
    Normal gather uses:
    - construct null and read (>>) from Istream
    - binary operator and assignment operator to combine values

    combineGather uses:
    - construct from Istream
    - modify operator which modifies its lhs

\*---------------------------------------------------------------------------*/

#include "OPstream.H"
#include "IPstream.H"
#include "IOstreams.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T, class CombineOp>
void Pstream::combineGather
(
    const List<UPstream::commsStruct>& comms,
    T& Value,
    const CombineOp& cop
)
{
    if (UPstream::parRun())
    {
        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo()];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (contiguous<T>())
            {
                T value;
                UIPstream::read
                (
                    UPstream::scheduled,
                    belowID,
                    reinterpret_cast<char*>(&value),
                    sizeof(T)
                );

                if (debug & 2)
                {
                    Pout<< " received from "
                        << belowID << " data:" << value << endl;
                }

                cop(Value, value);
            }
            else
            {
                IPstream fromBelow(UPstream::scheduled, belowID);
                T value(fromBelow);

                if (debug & 2)
                {
                    Pout<< " received from "
                        << belowID << " data:" << value << endl;
                }

                cop(Value, value);
            }
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            if (debug & 2)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << Value << endl;
            }

            if (contiguous<T>())
            {
                UOPstream::write
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T)
                );
            }
            else
            {
                OPstream toAbove(UPstream::scheduled, myComm.above());
                toAbove << Value;
            }
        }
    }
}


template <class T, class CombineOp>
void Pstream::combineGather(T& Value, const CombineOp& cop)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        combineGather(UPstream::linearCommunication(), Value, cop);
    }
    else
    {
        combineGather(UPstream::treeCommunication(), Value, cop);
    }
}


template <class T>
void Pstream::combineScatter
(
    const List<UPstream::commsStruct>& comms,
    T& Value
)
{
    if (UPstream::parRun())
    {
        // Get my communication order
        const UPstream::commsStruct& myComm = comms[UPstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            if (contiguous<T>())
            {
                UIPstream::read
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(&Value),
                    sizeof(T)
                );
            }
            else
            {
                IPstream fromAbove(UPstream::scheduled, myComm.above());
                Value = T(fromAbove);
            }

            if (debug & 2)
            {
                Pout<< " received from "
                    << myComm.above() << " data:" << Value << endl;
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (debug & 2)
            {
                Pout<< " sending to " << belowID << " data:" << Value << endl;
            }

            if (contiguous<T>())
            {
                UOPstream::write
                (
                    UPstream::scheduled,
                    belowID,
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T)
                );
            }
            else
            {
                OPstream toBelow(UPstream::scheduled, belowID);
                toBelow << Value;
            }
        }
    }
}


template <class T>
void Pstream::combineScatter(T& Value)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        combineScatter(UPstream::linearCommunication(), Value);
    }
    else
    {
        combineScatter(UPstream::treeCommunication(), Value);
    }
}


// Same thing but for whole list at a time
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


template <class T, class CombineOp>
void Pstream::listCombineGather
(
    const List<UPstream::commsStruct>& comms,
    List<T>& Values,
    const CombineOp& cop
)
{
    if (UPstream::parRun())
    {
        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo()];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (contiguous<T>())
            {
                List<T> receivedValues(Values.size());

                UIPstream::read
                (
                    UPstream::scheduled,
                    belowID,
                    reinterpret_cast<char*>(receivedValues.begin()),
                    receivedValues.byteSize()
                );

                if (debug & 2)
                {
                    Pout<< " received from "
                        << belowID << " data:" << receivedValues << endl;
                }

                forAll(Values, i)
                {
                    cop(Values[i], receivedValues[i]);
                }
            }
            else
            {
                IPstream fromBelow(UPstream::scheduled, belowID);
                List<T> receivedValues(fromBelow);

                if (debug & 2)
                {
                    Pout<< " received from "
                        << belowID << " data:" << receivedValues << endl;
                }

                forAll(Values, i)
                {
                    cop(Values[i], receivedValues[i]);
                }
            }
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            if (debug & 2)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << Values << endl;
            }

            if (contiguous<T>())
            {
                UOPstream::write
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(Values.begin()),
                    Values.byteSize()
                );
            }
            else
            {
                OPstream toAbove(UPstream::scheduled, myComm.above());
                toAbove << Values;
            }
        }
    }
}


template <class T, class CombineOp>
void Pstream::listCombineGather(List<T>& Values, const CombineOp& cop)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        listCombineGather(UPstream::linearCommunication(), Values, cop);
    }
    else
    {
        listCombineGather(UPstream::treeCommunication(), Values, cop);
    }
}


template <class T>
void Pstream::listCombineScatter
(
    const List<UPstream::commsStruct>& comms,
    List<T>& Values
)
{
    if (UPstream::parRun())
    {
        // Get my communication order
        const UPstream::commsStruct& myComm = comms[UPstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            if (contiguous<T>())
            {
                UIPstream::read
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(Values.begin()),
                    Values.byteSize()
                );
            }
            else
            {
                IPstream fromAbove(UPstream::scheduled, myComm.above());
                fromAbove >> Values;
            }

            if (debug & 2)
            {
                Pout<< " received from "
                    << myComm.above() << " data:" << Values << endl;
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (debug & 2)
            {
                Pout<< " sending to " << belowID << " data:" << Values << endl;
            }

            if (contiguous<T>())
            {
                UOPstream::write
                (
                    UPstream::scheduled,
                    belowID,
                    reinterpret_cast<const char*>(Values.begin()),
                    Values.byteSize()
                );
            }
            else
            {
                OPstream toBelow(UPstream::scheduled, belowID);
                toBelow << Values;
            }
        }
    }
}


template <class T>
void Pstream::listCombineScatter(List<T>& Values)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        listCombineScatter(UPstream::linearCommunication(), Values);
    }
    else
    {
        listCombineScatter(UPstream::treeCommunication(), Values);
    }
}




// Same thing but for sparse list (map)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


template <class Container, class CombineOp>
void Pstream::mapCombineGather
(
    const List<UPstream::commsStruct>& comms,
    Container& Values,
    const CombineOp& cop
)
{
    if (UPstream::parRun())
    {
        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo()];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            IPstream fromBelow(UPstream::scheduled, belowID);
            Container receivedValues(fromBelow);

            if (debug & 2)
            {
                Pout<< " received from "
                    << belowID << " data:" << receivedValues << endl;
            }

            for
            (
                typename Container::const_iterator slaveIter =
                    receivedValues.begin();
                slaveIter != receivedValues.end();
                ++slaveIter
            )
            {
                typename Container::iterator
                    masterIter = Values.find(slaveIter.key());

                if (masterIter != Values.end())
                {
                    cop(masterIter(), slaveIter());
                }
                else
                {
                    Values.insert(slaveIter.key(), slaveIter());
                }
            }
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            if (debug & 2)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << Values << endl;
            }

            OPstream toAbove(UPstream::scheduled, myComm.above());
            toAbove << Values;
        }
    }
}


template <class Container, class CombineOp>
void Pstream::mapCombineGather(Container& Values, const CombineOp& cop)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        mapCombineGather(UPstream::linearCommunication(), Values, cop);
    }
    else
    {
        mapCombineGather(UPstream::treeCommunication(), Values, cop);
    }
}


template <class Container>
void Pstream::mapCombineScatter
(
    const List<UPstream::commsStruct>& comms,
    Container& Values
)
{
    if (UPstream::parRun())
    {
        // Get my communication order
        const UPstream::commsStruct& myComm = comms[UPstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            IPstream fromAbove(UPstream::scheduled, myComm.above());
            fromAbove >> Values;

            if (debug & 2)
            {
                Pout<< " received from "
                    << myComm.above() << " data:" << Values << endl;
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (debug & 2)
            {
                Pout<< " sending to " << belowID << " data:" << Values << endl;
            }

            OPstream toBelow(UPstream::scheduled, belowID);
            toBelow << Values;
        }
    }
}


template <class Container>
void Pstream::mapCombineScatter(Container& Values)
{
    if (UPstream::nProcs() < UPstream::nProcsSimpleSum)
    {
        mapCombineScatter(UPstream::linearCommunication(), Values);
    }
    else
    {
        mapCombineScatter(UPstream::treeCommunication(), Values);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
