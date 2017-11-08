/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "error.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label
Foam::HashTable<T, Key, Hash>::hashKeyIndex(const Key& key) const
{
    // capacity is always a power of two - this is the modulus
    return Hash()(key) & (capacity_ - 1);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::capacity() const
{
    return capacity_;
}


template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::size() const
{
    return size_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::empty() const
{
    return !size_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::found(const Key& key) const
{
    if (size_)
    {
        return Iterator<true>(this, key).found();
    }

    return false;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::find
(
    const Key& key
)
{
    if (size_)
    {
        return iterator(Iterator<false>(this, key));
    }

    return iterator();
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::find
(
    const Key& key
) const
{
    return this->cfind(key);
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::cfind
(
    const Key& key
) const
{
    if (size_)
    {
        return const_iterator(Iterator<true>(this, key));
    }

    return const_iterator();
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::insert
(
    const Key& key,
    const T& obj
)
{
    return this->set(key, obj, false);  // No overwrite
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& obj
)
{
    return this->set(key, obj, true);   // Overwrite
}


template<class T, class Key, class Hash>
inline Foam::Xfer<Foam::HashTable<T, Key, Hash>>
Foam::HashTable<T, Key, Hash>::xfer()
{
    return xferMove(*this);
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::lookup
(
    const Key& key,
    const T& deflt
) const
{
    const const_iterator iter(this->cfind(key));
    return iter.found() ? iter.object() : deflt;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key)
{
    const iterator iter(this->find(key));

    if (!iter.found())
    {
        FatalErrorInFunction
            << key << " not found in table.  Valid entries: "
            << toc()
            << exit(FatalError);
    }

    return iter.object();
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key) const
{
    const const_iterator iter(this->cfind(key));

    if (!iter.found())
    {
        FatalErrorInFunction
            << key << " not found in table.  Valid entries: "
            << toc()
            << exit(FatalError);
    }

    return iter.object();
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator()(const Key& key)
{
    const iterator iter(this->find(key));

    if (iter.found())
    {
        return iter.object();
    }

    this->insert(key, mapped_type());
    return find(key).object();
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator()
(
    const Key& key,
    const T& deflt
)
{
    const iterator iter(this->find(key));

    if (iter.found())
    {
        return iter.object();
    }

    this->insert(key, deflt);
    return find(key).object();
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::operator()
(
    const Key& key,
    const T& deflt
) const
{
    return this->lookup(key, deflt);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline void Foam::Swap
(
    HashTable<T, Key, Hash>& a,
    HashTable<T, Key, Hash>& b
)
{
    a.swap(b);
}


// ************************************************************************* //
