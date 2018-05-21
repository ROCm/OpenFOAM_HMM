/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashPtrTable<T, Key, Hash>::HashPtrTable()
:
    parent_type()
{}


template<class T, class Key, class Hash>
inline Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(const label size)
:
    parent_type(size)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::insert
(
    const Key& key,
    autoPtr<T>& aptr
)
{
    if (parent_type::insert(key, aptr.get()))
    {
        aptr.release(); // Now owned by HashPtrTable
        return true;
    }

    return false;
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::insert
(
    const Key& key,
    autoPtr<T>&& aptr
)
{
    if (parent_type::insert(key, aptr.get()))
    {
        aptr.release(); // Now owned by HashPtrTable
        return true;
    }

    return false;
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    T* ptr
)
{
    return this->parent_type::set(key, ptr);
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    autoPtr<T>& aptr
)
{
    return this->set(key, aptr.release());
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    autoPtr<T>&& aptr
)
{
    return this->set(key, aptr.release());
}


// ************************************************************************* //
