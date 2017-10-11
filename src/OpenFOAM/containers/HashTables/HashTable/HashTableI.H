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

// * * * * * * * * * * * * * Private Member Classes * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::hashedEntry::hashedEntry
(
    const Key& key,
    const T& obj,
    hashedEntry* next
)
:
    key_(key),
    obj_(obj),
    next_(next)
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label
Foam::HashTable<T, Key, Hash>::hashKeyIndex(const Key& key) const
{
    // size is power of two - this is the modulus
    return Hash()(key) & (tableSize_ - 1);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::capacity() const
{
    return tableSize_;
}


template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::size() const
{
    return nElmts_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::empty() const
{
    return !nElmts_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::insert
(
    const Key& key,
    const T& obj
)
{
    return this->set(key, obj, true);
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& obj
)
{
    return this->set(key, obj, false);
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
    const_iterator iter = this->find(key);
    return iter.found() ? iter.object() : deflt;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key)
{
    iterator iter = this->find(key);

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
    const_iterator iter = this->find(key);

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
    iterator iter = this->find(key);

    if (iter.found())
    {
        return iter.object();
    }

    this->insert(key, T());
    return find(key).object();
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator()
(
    const Key& key,
    const T& deflt
)
{
    iterator iter = this->find(key);

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


// * * * * * * * * * * * * * * * iterator base * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator_base::iterator_base()
:
    entryPtr_(nullptr),
    hashTable_(nullptr),
    hashIndex_(0)
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator_base::iterator_base
(
    const table_type* hashTbl,
    const entry_type* elmt,
    const label hashIndex
)
:
    entryPtr_(const_cast<entry_type*>(elmt)),
    hashTable_(const_cast<table_type*>(hashTbl)),
    hashIndex_(hashIndex)
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator_base::iterator_base
(
    const table_type* hashTbl
)
:
    entryPtr_(nullptr),
    hashTable_(const_cast<table_type*>(hashTbl)),
    hashIndex_(0)
{
    if (hashTable_ && hashTable_->nElmts_)
    {
        // find first non-nullptr table entry
        while
        (
            !(entryPtr_ = hashTable_->table_[hashIndex_])
         && ++hashIndex_ < hashTable_->tableSize_
        )
        {}

        if (hashIndex_ >= hashTable_->tableSize_)
        {
            // make into an end iterator
            entryPtr_ = nullptr;
            hashIndex_ = 0;
        }
    }
}


template<class T, class Key, class Hash>
inline void
Foam::HashTable<T, Key, Hash>::iterator_base::increment()
{
    // A negative index is a special value from erase
    if (hashIndex_ < 0)
    {
        // the markPos='-curPos-1', but we wish to continue at 'curPos-1'
        // thus use '-(markPos+1) -1'
        hashIndex_ = -(hashIndex_+1) - 1;
    }
    else if (entryPtr_)
    {
        if (entryPtr_->next_)
        {
            // Move to next element on the SLList
            entryPtr_ = entryPtr_->next_;
            return;
        }
    }
    // else
    // {
    //     // if we reach here (entryPtr_ is nullptr) it is already at the end()
    //     // we should probably stop
    // }


    // Step to the next table entry
    while
    (
        ++hashIndex_ < hashTable_->tableSize_
     && !(entryPtr_ = hashTable_->table_[hashIndex_])
    )
    {}

    if (hashIndex_ >= hashTable_->tableSize_)
    {
        // make into an end iterator
        entryPtr_ = nullptr;
        hashIndex_ = 0;
    }
}


template<class T, class Key, class Hash>
inline bool
Foam::HashTable<T, Key, Hash>::iterator_base::found() const
{
    return entryPtr_;
}


template<class T, class Key, class Hash>
inline const Key& Foam::HashTable<T, Key, Hash>::iterator_base::key() const
{
    return entryPtr_->key_;
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::iterator_base::element() const
{
    return entryPtr_->obj_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::iterator_base::operator==
(
    const iterator_base& iter
) const
{
    return entryPtr_ == iter.entryPtr_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::iterator_base::operator!=
(
    const iterator_base& iter
) const
{
    return entryPtr_ != iter.entryPtr_;
}


// * * * * * * * * * * * * * * key iterator base * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class WrappedIterator>
inline Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::key_iterator_base
(
    const WrappedIterator& iter
)
:
    WrappedIterator(iter)
{}


template<class T, class Key, class Hash>
template<class WrappedIterator>
inline const Key&
Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::operator*() const
{
    return this->key();
}


template<class T, class Key, class Hash>
template<class WrappedIterator>
inline const Key&
Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::operator()() const
{
    return this->key();
}


template<class T, class Key, class Hash>
template<class WrappedIterator>
inline Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>&
Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::operator++()
{
    this->increment();
    return *this;
}


template<class T, class Key, class Hash>
template<class WrappedIterator>
inline Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
Foam::HashTable<T, Key, Hash>::key_iterator_base<WrappedIterator>
::operator++(int)
{
    key_iterator_base old = *this;
    this->increment();
    return old;
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator::iterator()
:
    iterator_base()
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator::iterator
(
    table_type* hashTbl
)
:
    iterator_base(hashTbl)
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator::iterator
(
    table_type* hashTbl,
    entry_type* elmt,
    const label hashIndex
)
:
    iterator_base(hashTbl, elmt, hashIndex)
{}


template<class T, class Key, class Hash>
inline T&
Foam::HashTable<T, Key, Hash>::iterator::object() const
{
    return this->element();
}


template<class T, class Key, class Hash>
inline T&
Foam::HashTable<T, Key, Hash>::iterator::operator*() const
{
    return this->object();
}


template<class T, class Key, class Hash>
inline T&
Foam::HashTable<T, Key, Hash>::iterator::operator()() const
{
    return this->object();
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator&
Foam::HashTable<T, Key, Hash>::iterator::operator++()
{
    this->increment();
    return *this;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::iterator::operator++(int)
{
    iterator old = *this;
    this->increment();
    return old;
}


// * * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::const_iterator::const_iterator()
:
    iterator_base()
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::const_iterator::const_iterator
(
    const HashTable<T, Key, Hash>::iterator& iter
)
:
    iterator_base(iter)
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::const_iterator::const_iterator
(
    table_type* hashTbl
)
:
    iterator_base(hashTbl)
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::const_iterator::const_iterator
(
    table_type* hashTbl,
    entry_type* elmt,
    const label hashIndex
)
:
    iterator_base(hashTbl, elmt, hashIndex)
{}


template<class T, class Key, class Hash>
inline const T&
Foam::HashTable<T, Key, Hash>::const_iterator::object() const
{
    return this->element();
}


template<class T, class Key, class Hash>
inline const T&
Foam::HashTable<T, Key, Hash>::const_iterator::operator*() const
{
    return this->object();
}


template<class T, class Key, class Hash>
inline const T&
Foam::HashTable<T, Key, Hash>::const_iterator::operator()() const
{
    return this->object();
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator&
Foam::HashTable<T, Key, Hash>::const_iterator::operator++()
{
    this->increment();
    return *this;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::const_iterator::operator++(int)
{
    const_iterator old = *this;
    this->increment();
    return old;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::begin()
{
    return iterator(this);
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::begin() const
{
    return const_iterator(this);
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::cbegin() const
{
    return const_iterator(this);
}


template<class T, class Key, class Hash>
inline const typename Foam::HashTable<T, Key, Hash>::iterator&
Foam::HashTable<T, Key, Hash>::end()
{
    return iterator_end<iterator>();
}


template<class T, class Key, class Hash>
inline const typename Foam::HashTable<T, Key, Hash>::const_iterator&
Foam::HashTable<T, Key, Hash>::end() const
{
    return iterator_end<const_iterator>();
}


template<class T, class Key, class Hash>
inline const typename Foam::HashTable<T, Key, Hash>::const_iterator&
Foam::HashTable<T, Key, Hash>::cend() const
{
    return iterator_cend<const_iterator>();
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
