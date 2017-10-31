/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#ifndef HashSet_C
#define HashSet_C

#include "HashSet.H"
#include "FixedList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Key, class Hash>
template<class InputIter>
inline Foam::label Foam::HashSet<Key, Hash>::insertMultiple
(
    const InputIter begIter,
    const InputIter endIter
)
{
    label changed = 0;
    for (InputIter iter = begIter; iter != endIter; ++iter)
    {
        if (insert(*iter))
        {
            ++changed;
        }
    }
    return changed;
}


template<class Key, class Hash>
template<class InputIter>
inline Foam::label Foam::HashSet<Key, Hash>::assignMultiple
(
    const InputIter begIter,
    const InputIter endIter,
    const label sz
)
{
    if (!this->capacity())
    {
        // Could be zero-sized from a previous transfer()?
        this->resize(sz);
    }
    else
    {
        this->clear();
    }

    return insertMultiple(begIter, endIter);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Key, class Hash>
Foam::HashSet<Key, Hash>::HashSet(const UList<Key>& lst)
:
    parent_type(2*lst.size())
{
    for (const auto& k : lst)
    {
        this->insert(k);
    }
}


template<class Key, class Hash>
template<unsigned Size>
Foam::HashSet<Key, Hash>::HashSet(const FixedList<Key, Size>& lst)
:
    parent_type(2*lst.size())
{
    for (const auto& k : lst)
    {
        this->insert(k);
    }
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>::HashSet(std::initializer_list<Key> lst)
:
    parent_type(2*lst.size())
{
    for (const auto& k : lst)
    {
        this->insert(k);
    }
}


template<class Key, class Hash>
template<class AnyType, class AnyHash>
Foam::HashSet<Key, Hash>::HashSet
(
    const HashTable<AnyType, Key, AnyHash>& tbl
)
:
    parent_type(tbl.capacity())
{
    using other_iter =
        typename HashTable<AnyType, Key, AnyHash>::const_iterator;

    for (other_iter iter = tbl.cbegin(); iter != tbl.cend(); ++iter)
    {
        this->insert(iter.key());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Key, class Hash>
Foam::label Foam::HashSet<Key, Hash>::insert(const UList<Key>& lst)
{
    return insertMultiple(lst.begin(), lst.end());
}


template<class Key, class Hash>
template<unsigned Size>
Foam::label Foam::HashSet<Key, Hash>::insert(const FixedList<Key, Size>& lst)
{
    return insertMultiple(lst.begin(), lst.end());
}


template<class Key, class Hash>
Foam::label Foam::HashSet<Key, Hash>::insert(std::initializer_list<Key> lst)
{
    return insertMultiple(lst.begin(), lst.end());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Key, class Hash>
inline bool Foam::HashSet<Key, Hash>::operator()(const Key& key) const
{
    return this->found(key);
}


template<class Key, class Hash>
inline bool Foam::HashSet<Key, Hash>::operator[](const Key& key) const
{
    return this->found(key);
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator=(const UList<Key>& lst)
{
    assignMultiple(lst.begin(), lst.end(), 2*lst.size());
}


template<class Key, class Hash>
template<unsigned Size>
void Foam::HashSet<Key, Hash>::operator=(const FixedList<Key, Size>& lst)
{
    assignMultiple(lst.begin(), lst.end(), 2*lst.size());
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator=(std::initializer_list<Key> lst)
{
    assignMultiple(lst.begin(), lst.end(), 2*lst.size());
}



template<class Key, class Hash>
bool Foam::HashSet<Key, Hash>::operator==(const HashSet<Key, Hash>& rhs) const
{
    // Sizes (number of keys) must match
    if (this->size() != rhs.size())
    {
        return false;
    }

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        if (!this->found(iter.key()))
        {
            return false;
        }
    }

    return true;
}


template<class Key, class Hash>
bool Foam::HashSet<Key, Hash>::operator!=(const HashSet<Key, Hash>& rhs) const
{
    return !operator==(rhs);
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator|=(const HashSet<Key, Hash>& rhs)
{
    // Add rhs elements into lhs
    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        this->insert(iter.key());
    }
}


template<class Key, class Hash>
inline void Foam::HashSet<Key, Hash>::operator&=(const HashSet<Key, Hash>& rhs)
{
    this->parent_type::retain(rhs);
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator^=(const HashSet<Key, Hash>& rhs)
{
    // Add missed rhs elements, remove duplicate elements
    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        if (this->found(iter.key()))
        {
            this->erase(iter.key());
        }
        else
        {
            this->insert(iter.key());
        }
    }
}


template<class Key, class Hash>
inline void Foam::HashSet<Key, Hash>::operator-=(const HashSet<Key, Hash>& rhs)
{
    this->parent_type::erase(rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Key, class Hash>
Foam::Ostream& Foam::operator<<(Ostream& os, const HashSet<Key, Hash>& tbl)
{
    return tbl.writeList(os, 10);  // 10=consistent with UList
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Key, class Hash>
inline void Foam::Swap
(
    HashSet<Key, Hash>& a,
    HashSet<Key, Hash>& b
)
{
    a.swap(b);
}


/* * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * */

template<class Key, class Hash>
Foam::HashSet<Key, Hash>
Foam::operator|
(
    const HashSet<Key, Hash>& hash1,
    const HashSet<Key, Hash>& hash2
)
{
    HashSet<Key, Hash> out(hash1);
    out |= hash2;
    return out;
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>
Foam::operator&
(
    const HashSet<Key, Hash>& hash1,
    const HashSet<Key, Hash>& hash2
)
{
    HashSet<Key, Hash> out(hash1);
    out &= hash2;
    return out;
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>
Foam::operator^
(
    const HashSet<Key, Hash>& hash1,
    const HashSet<Key, Hash>& hash2
)
{
    HashSet<Key, Hash> out(hash1);
    out ^= hash2;
    return out;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Key, class Hash>
inline typename Foam::HashSet<Key, Hash>::iterator
Foam::HashSet<Key, Hash>::begin()
{
    return HashTableCore::iterator_begin<iterator>
    (
        static_cast<parent_type&>(*this)
    );
}


template<class Key, class Hash>
inline typename Foam::HashSet<Key, Hash>::const_iterator
Foam::HashSet<Key, Hash>::begin() const
{
    return HashTableCore::iterator_cbegin<const_iterator>
    (
        static_cast<const parent_type&>(*this)
    );
}


template<class Key, class Hash>
inline typename Foam::HashSet<Key, Hash>::const_iterator
Foam::HashSet<Key, Hash>::cbegin() const
{
    return HashTableCore::iterator_cbegin<const_iterator>
    (
        static_cast<const parent_type&>(*this)
    );
}


template<class Key, class Hash>
inline const typename Foam::HashSet<Key, Hash>::iterator&
Foam::HashSet<Key, Hash>::end()
{
    return HashTableCore::iterator_end<iterator>();
}


template<class Key, class Hash>
inline const typename Foam::HashSet<Key, Hash>::const_iterator&
Foam::HashSet<Key, Hash>::end() const
{
    return HashTableCore::iterator_cend<const_iterator>();
}


template<class Key, class Hash>
inline const typename Foam::HashSet<Key, Hash>::const_iterator&
Foam::HashSet<Key, Hash>::cend() const
{
    return HashTableCore::iterator_cend<const_iterator>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
