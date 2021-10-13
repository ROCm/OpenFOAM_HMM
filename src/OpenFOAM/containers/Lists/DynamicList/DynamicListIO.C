/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "List.H"
#include "Istream.H"
#include "token.H"
#include "SLList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::DynamicList<T, SizeMin>::DynamicList(Istream& is)
:
    List<T>(),
    capacity_(0)
{
    this->readList(is);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::Istream& Foam::DynamicList<T, SizeMin>::readList
(
    Istream& is
)
{
    DynamicList<T, SizeMin>& list = *this;

    // Needs rewrite (2021-10)
    // Use entire storage - ie, resize(capacity())
    (void) list.expandStorage();

    static_cast<List<T>&>(list).readList(is);
    list.capacity_ = list.size();

    return is;
}


// ************************************************************************* //
