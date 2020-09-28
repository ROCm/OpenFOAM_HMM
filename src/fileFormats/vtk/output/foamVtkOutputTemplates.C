/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "Pstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void Foam::vtk::write
(
    vtk::formatter& fmt,
    const Type& val,
    const label n
)
{
    const direction nCmpt = pTraits<Type>::nComponents;

    for (label i=0; i < n; ++i)
    {
        for (direction cmpt=0; cmpt < nCmpt; ++cmpt)
        {
            fmt.write(component(val, cmpt));
        }
    }
}


template<class Type>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<Type>& values
)
{
    for (const Type& val : values)
    {
        vtk::write(fmt, val);
    }
}


template<class Type, unsigned N>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const FixedList<Type, N>& values
)
{
    for (const Type& val : values)
    {
        vtk::write(fmt, val);
    }
}


template<class Type>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const labelUList& addressing
)
{
    for (const label idx : addressing)
    {
        vtk::write(fmt, values[idx]);
    }
}


template<class Type>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const bitSet& selected
)
{
    for (const label idx : selected)
    {
        vtk::write(fmt, values[idx]);
    }
}


template<class Type>
void Foam::vtk::writeLists
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const UList<Type>& indirect,
    const labelUList& addressing
)
{
    vtk::writeList(fmt, values);
    vtk::writeList(fmt, indirect, addressing);
}


template<class Type>
void Foam::vtk::writeListParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values
)
{
    if (Pstream::master())
    {
        vtk::writeList(fmt, values);

        List<Type> recv;

        // Receive and write
        for (const int slave : Pstream::subProcs())
        {
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);

            fromSlave >> recv;

            vtk::writeList(fmt, recv);
        }
    }
    else
    {
        // Send to master
        OPstream toMaster
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo()
        );

        toMaster << values;
    }
}


template<class Type>
void Foam::vtk::writeListParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const labelUList& addressing
)
{
    if (Pstream::master())
    {
        vtk::writeList(fmt, values, addressing);

        List<Type> recv;

        // Receive and write
        for (const int slave : Pstream::subProcs())
        {
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);

            fromSlave >> recv;

            vtk::writeList(fmt, recv);
        }
    }
    else
    {
        // Send to master
        OPstream toMaster
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo()
        );

        toMaster << List<Type>(values, addressing);
    }
}


template<class Type>
void Foam::vtk::writeListParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const bitSet& selected
)
{
    if (Pstream::master())
    {
        vtk::writeList(fmt, values, selected);

        List<Type> recv;

        // Receive and write
        for (const int slave : Pstream::subProcs())
        {
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);

            fromSlave >> recv;

            vtk::writeList(fmt, recv);
        }
    }
    else
    {
        // Send to master
        OPstream toMaster
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo()
        );

        toMaster << subset(selected, values);
    }
}


template<class Type>
void Foam::vtk::writeListsParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values1,
    const UList<Type>& values2
)
{
    if (Pstream::master())
    {
        vtk::writeList(fmt, values1);
        vtk::writeList(fmt, values2);

        List<Type> recv1, recv2;

        // Receive and write
        for (const int slave : Pstream::subProcs())
        {
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);

            fromSlave >> recv1 >> recv2;

            vtk::writeList(fmt, recv1);
            vtk::writeList(fmt, recv2);
        }
    }
    else
    {
        // Send to master
        OPstream toMaster
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo()
        );

        toMaster << values1 << values2;
    }
}


template<class Type>
void Foam::vtk::writeListsParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values1,
    const UList<Type>& values2,
    const labelUList& addressing
)
{
    if (Pstream::master())
    {
        vtk::writeList(fmt, values1);
        vtk::writeList(fmt, values2, addressing);

        List<Type> recv1, recv2;

        // Receive and write
        for (const int slave : Pstream::subProcs())
        {
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);

            fromSlave >> recv1 >> recv2;

            vtk::writeList(fmt, recv1);
            vtk::writeList(fmt, recv2);
        }
    }
    else
    {
        // Send to master
        OPstream toMaster
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo()
        );

        toMaster << values1 << List<Type>(values2, addressing);
    }
}


// ************************************************************************* //
