/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "OPstream.H"
#include "IPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
void Foam::Pstream::genericBroadcast(T& value, const label comm)
{
    // Generic: use stream interface
    if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        if (UPstream::master(comm))
        {
            OPBstream toAll(UPstream::masterNo(), comm);
            toAll << value;
        }
        else
        {
            IPBstream fromMaster(UPstream::masterNo(), comm);
            fromMaster >> value;
        }
    }
}


template<class T>
void Foam::Pstream::broadcast(T& value, const label comm)
{
    if (!is_contiguous<T>::value)
    {
        Pstream::genericBroadcast(value, comm);
    }
    else if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        UPstream::broadcast
        (
            reinterpret_cast<char*>(&value),
            sizeof(T),
            comm,
            UPstream::masterNo()
        );
    }
}


template<class T>
void Foam::Pstream::broadcast(List<T>& values, const label comm)
{
    if (!is_contiguous<T>::value)
    {
        Pstream::genericBroadcast(values, comm);
    }
    else if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        // Broadcast the size
        label len(values.size());
        UPstream::broadcast
        (
            reinterpret_cast<char*>(&len),
            sizeof(label),
            comm,
            UPstream::masterNo()
        );
        values.resize_nocopy(len);  // A no-op on master

        if (len)
        {
            UPstream::broadcast
            (
                values.data_bytes(),
                values.size_bytes(),
                comm,
                UPstream::masterNo()
            );
        }
    }
}


template<class T, int SizeMin>
void Foam::Pstream::broadcast(DynamicList<T, SizeMin>& values, const label comm)
{
    if (!is_contiguous<T>::value)
    {
        Pstream::genericBroadcast(values, comm);
    }
    else if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        // Broadcast the size
        label len(values.size());
        UPstream::broadcast
        (
            reinterpret_cast<char*>(&len),
            sizeof(label),
            comm,
            UPstream::masterNo()
        );
        values.resize_nocopy(len);  // A no-op on master

        if (len)
        {
            UPstream::broadcast
            (
                values.data_bytes(),
                values.size_bytes(),
                comm,
                UPstream::masterNo()
            );
        }
    }
}


// ************************************************************************* //
