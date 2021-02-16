/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
inline Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(const label size)
:
    parent_type(size)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline const T* Foam::HashPtrTable<T, Key, Hash>::get(const Key& key) const
{
    const const_iterator iter(this->cfind(key));
    if (iter.good())
    {
        return iter.val();
    }
    return nullptr;
}


template<class T, class Key, class Hash>
template<class... Args>
inline bool Foam::HashPtrTable<T, Key, Hash>::emplace
(
    const Key& key,
    Args&&... args
)
{
    // Use insertion semantics
    return
    (
        !parent_type::found(key)
     && this->parent_type::set(key, new T(std::forward<Args>(args)...))
    );
}


template<class T, class Key, class Hash>
template<class... Args>
inline bool Foam::HashPtrTable<T, Key, Hash>::emplace_set
(
    const Key& key,
    Args&&... args
)
{
    return this->set(key, new T(std::forward<Args>(args)...));
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::insert
(
    const Key& key,
    autoPtr<T>& ptr
)
{
    if (parent_type::insert(key, ptr.get()))
    {
        ptr.release(); // Now owned by HashPtrTable
        return true;
    }

    return false;
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::insert
(
    const Key& key,
    autoPtr<T>&& ptr
)
{
    if (parent_type::insert(key, ptr.get()))
    {
        ptr.release(); // Now owned by HashPtrTable
        return true;
    }

    return false;
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::insert
(
    const Key& key,
    std::unique_ptr<T>&& ptr
)
{
    if (parent_type::insert(key, ptr.get()))
    {
        ptr.release(); // Now owned by HashPtrTable
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
    const T* old = this->get(key);

    const bool ok = this->parent_type::set(key, ptr);

    if (ok && old != ptr)
    {
        delete const_cast<T*>(old);
    }

    return ok;
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    autoPtr<T>& ptr
)
{
    return this->set(key, ptr.release());
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    autoPtr<T>&& ptr
)
{
    return this->set(key, ptr.release());
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    std::unique_ptr<T>&& ptr
)
{
    return this->set(key, ptr.release());
}


// ************************************************************************* //
