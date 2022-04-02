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

template<class Type>
void Foam::Pstream::broadcast(Type& value, const label comm)
{
    if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        if (is_contiguous<Type>::value)
        {
            // Note: contains parallel guard internally as well
            UPstream::broadcast
            (
                reinterpret_cast<char*>(&value),
                sizeof(Type),
                comm,
                UPstream::masterNo()
            );
        }
        else
        {
            if (UPstream::master(comm))
            {
                OPBstream os(UPstream::masterNo(), comm);
                os << value;
            }
            else
            {
                IPBstream is(UPstream::masterNo(), comm);
                is >> value;
            }
        }
    }
}


template<class Type, class... Args>
void Foam::Pstream::broadcasts(const label comm, Type& arg1, Args&&... args)
{
    if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        if (UPstream::master(comm))
        {
            OPBstream os(UPstream::masterNo(), comm);
            Detail::outputLoop(os, arg1, std::forward<Args>(args)...);
        }
        else
        {
            IPBstream is(UPstream::masterNo(), comm);
            Detail::inputLoop(is, arg1, std::forward<Args>(args)...);
        }
    }
}


// ************************************************************************* //
