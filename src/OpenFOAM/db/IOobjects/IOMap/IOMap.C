/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "IOMap.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
bool Foam::IOMap<T>::readContents()
{
    if
    (
        (
            readOpt() == IOobject::MUST_READ
         || readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        readStream(typeName) >> *this;
        close();

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io)
:
    regIOobject(io)
{
    readContents();
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const label size)
:
    regIOobject(io)
{
    if (!readContents())
    {
        Map<T>::resize(size);
    }
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const Map<T>& content)
:
    regIOobject(io)
{
    if (!readContents())
    {
        Map<T>::operator=(content);
    }
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, Map<T>&& content)
:
    regIOobject(io)
{
    Map<T>::transfer(content);

    readContents();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::IOMap<T>::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::IOMap<T>::operator=(const IOMap<T>& rhs)
{
    Map<T>::operator=(rhs);
}


// ************************************************************************* //
