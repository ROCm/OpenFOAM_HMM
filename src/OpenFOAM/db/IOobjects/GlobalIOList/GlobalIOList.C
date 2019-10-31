/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "GlobalIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    readHeaderOk(IOstream::BINARY, typeName);
}


template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList(const IOobject& io, const label size)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    if (!readHeaderOk(IOstream::BINARY, typeName))
    {
        List<Type>::setSize(size);
    }
}


template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList
(
    const IOobject& io,
    const UList<Type>& content
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    if (!readHeaderOk(IOstream::BINARY, typeName))
    {
        List<Type>::operator=(content);
    }
}


template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList
(
    const IOobject& io,
    List<Type>&& content
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    List<Type>::transfer(content);

    readHeaderOk(IOstream::BINARY, typeName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::GlobalIOList<Type>::readData(Istream& is)
{
    is >> *this;
    return is.good();
}


template<class Type>
bool Foam::GlobalIOList<Type>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::GlobalIOList<Type>::operator=(const GlobalIOList<Type>& rhs)
{
    List<Type>::operator=(rhs);
}


// ************************************************************************* //
