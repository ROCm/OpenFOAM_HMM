/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#ifndef HashTbl_C
#define HashTbl_C

#include "HashTbl.H"
#include "List.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::label Foam::HashTbl<T, Key, Hash>::canonicalSize(const label size)
{
    if (size < 1)
    {
        return 0;
    }

    // enforce power of two
    unsigned int goodSize = size;

    if (goodSize & (goodSize - 1))
    {
        // brute-force is fast enough
        goodSize = 1;
        while (goodSize < unsigned(size))
        {
            goodSize <<= 1;
        }
    }

    return goodSize;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTbl<T, Key, Hash>::HashTbl(const label size)
:
    HashTblName(),
    nElmts_(0),
    tableSize_(canonicalSize(size)),
    table_(NULL)
{
    if (tableSize_)
    {
        table_ = new hashedEntry*[tableSize_];

        for (label hashIdx = 0; hashIdx < tableSize_; hashIdx++)
        {
            table_[hashIdx] = 0;
        }
    }
}


template<class T, class Key, class Hash>
Foam::HashTbl<T, Key, Hash>::HashTbl(const HashTbl<T, Key, Hash>& ht)
:
    HashTblName(),
    nElmts_(0),
    tableSize_(ht.tableSize_),
    table_(NULL)
{
    if (tableSize_)
    {
        table_ = new hashedEntry*[tableSize_];

        for (label hashIdx = 0; hashIdx < tableSize_; hashIdx++)
        {
            table_[hashIdx] = 0;
        }

        for (const_iterator iter = ht.cbegin(); iter != ht.cend(); ++iter)
        {
            insert(iter.key(), *iter);
        }
    }
}

template<class T, class Key, class Hash>
Foam::HashTbl<T, Key, Hash>::HashTbl
(
    const Xfer<HashTbl<T, Key, Hash> >& ht
)
:
    HashTblName(),
    nElmts_(0),
    tableSize_(0),
    table_(NULL)
{
    transfer(ht());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTbl<T, Key, Hash>::~HashTbl()
{
    if (table_)
    {
        clear();
        delete[] table_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
bool Foam::HashTbl<T, Key, Hash>::found(const Key& key) const
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);

        for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
        {
            if (key == ep->key_)
            {
                return true;
            }
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "HashTbl<T, Key, Hash>::found(const Key& key) : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return false;
}


template<class T, class Key, class Hash>
typename Foam::HashTbl<T, Key, Hash>::iterator
Foam::HashTbl<T, Key, Hash>::find
(
    const Key& key
)
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);

        for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
        {
            if (key == ep->key_)
            {
                return iterator(this, ep, hashIdx);
            }
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "HashTbl<T, Key, Hash>::find(const Key& key) : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return end();
}


template<class T, class Key, class Hash>
typename Foam::HashTbl<T, Key, Hash>::const_iterator
Foam::HashTbl<T, Key, Hash>::find
(
    const Key& key
) const
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);

        for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
        {
            if (key == ep->key_)
            {
                return const_iterator(this, ep, hashIdx);
            }
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "HashTbl<T, Key, Hash>::find(const Key& key) const : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return cend();
}


template<class T, class Key, class Hash>
Foam::List<Key> Foam::HashTbl<T, Key, Hash>::toc() const
{
    List<Key> tofc(nElmts_);
    label i = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        tofc[i++] = iter.key();
    }

    return tofc;
}


template<class T, class Key, class Hash>
Foam::List<Key> Foam::HashTbl<T, Key, Hash>::sortedToc() const
{
    List<Key> sortedList = this->toc();
    sort(sortedList);

    return sortedList;
}


template<class T, class Key, class Hash>
bool Foam::HashTbl<T, Key, Hash>::set
(
    const Key& key,
    const T& newEntry,
    const bool protect
)
{
    if (!tableSize_)
    {
        resize(2);
    }

    const label hashIdx = hashKeyIndex(key);

    hashedEntry* existing = 0;
    hashedEntry* prev = 0;

    for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
    {
        if (key == ep->key_)
        {
            existing = ep;
            break;
        }
        prev = ep;
    }

    // not found, insert it at the head
    if (!existing)
    {
        table_[hashIdx] = new hashedEntry(key, table_[hashIdx], newEntry);
        nElmts_++;

        if (double(nElmts_)/tableSize_ > 0.8)
        {
#           ifdef FULLDEBUG
            if (debug)
            {
                Info<< "HashTbl<T, Key, Hash>::set"
                    "(const Key& key, T newEntry) : "
                    "Doubling table size\n";
            }
#           endif

            resize(2*tableSize_);
        }
    }
    else if (protect)
    {
        // found - but protected from overwriting
        // this corresponds to the STL 'insert' convention
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTbl<T, Key, Hash>::set"
                "(const Key& key, T newEntry, true) : "
                "Cannot insert " << key << " already in hash table\n";
        }
#       endif
        return false;
    }
    else
    {
        // found - overwrite existing entry
        // this corresponds to the Perl convention
        hashedEntry* ep = new hashedEntry(key, existing->next_, newEntry);

        // replace existing element - within list or insert at the head
        if (prev)
        {
            prev->next_ = ep;
        }
        else
        {
            table_[hashIdx] = ep;
        }

        delete existing;
    }

    return true;
}


template<class T, class Key, class Hash>
bool Foam::HashTbl<T, Key, Hash>::iteratorBase::erase()
{
    // note: entryPtr_ is NULL for end(), so this catches that too
    if (entryPtr_)
    {
        // Search element before entryPtr_
        hashedEntry* prev = 0;

        for
        (
            hashedEntry* ep = hashTable_->table_[hashIndex_];
            ep;
            ep = ep->next_
        )
        {
            if (ep == entryPtr_)
            {
                break;
            }
            prev = ep;
        }

        if (prev)
        {
            // has an element before entryPtr - reposition to there
            prev->next_ = entryPtr_->next_;
            delete entryPtr_;
            entryPtr_ = prev;
        }
        else
        {
            // entryPtr was first element on SLList
            hashTable_->table_[hashIndex_] = entryPtr_->next_;
            delete entryPtr_;

            // assign any non-NULL pointer value so it doesn't look
            // like end()/cend()
            entryPtr_ = reinterpret_cast<hashedEntry*>(this);

            // Mark with special hashIndex value to signal it has been rewound.
            // The next increment will bring it back to the present location.
            //
            // For the current position 'X', mark it as '-(X+1)', which is
            // written as '-X-1' to avoid overflow.
            // The extra '-1' is needed to avoid ambiguity for position '0'.
            // To retrieve the previous position 'X-1' we would later
            // use '-(-X-1) - 2'
            hashIndex_ = -hashIndex_ - 1;
        }

        hashTable_->nElmts_--;

        return true;
    }
    else
    {
        return false;
    }
}



// NOTE:
// We use (const iterator&) here, but manipulate its contents anyhow.
// The parameter should be (iterator&), but then the compiler doesn't find
// it correctly and tries to call as (iterator) instead.
//
template<class T, class Key, class Hash>
bool Foam::HashTbl<T, Key, Hash>::erase(const iterator& cit)
{
    // adjust iterator after erase
    return const_cast<iterator&>(cit).erase();
}


template<class T, class Key, class Hash>
bool Foam::HashTbl<T, Key, Hash>::erase(const Key& key)
{
    return erase(find(key));
}


template<class T, class Key, class Hash>
Foam::label Foam::HashTbl<T, Key, Hash>::erase(const UList<Key>& keys)
{
    const label nTotal = nElmts_;
    label count = 0;

    // Remove listed keys from this table - terminates early if possible
    for (label keyI = 0; count < nTotal && keyI < keys.size(); ++keyI)
    {
        if (erase(keys[keyI]))
        {
            count++;
        }
    }

    return count;
}


template<class T, class Key, class Hash>
template<class AnyType, class AnyHash>
Foam::label Foam::HashTbl<T, Key, Hash>::erase
(
    const HashTbl<AnyType, Key, AnyHash>& rhs
)
{
    label count = 0;

    // Remove rhs keys from this table - terminates early if possible
    // Could optimize depending on which hash is smaller ...
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (rhs.found(iter.key()) && erase(iter))
        {
            count++;
        }
    }

    return count;
}


template<class T, class Key, class Hash>
void Foam::HashTbl<T, Key, Hash>::resize(const label sz)
{
    label newSize = canonicalSize(sz);

    if (newSize == tableSize_)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTbl<T, Key, Hash>::resize(const label) : "
                << "new table size == old table size\n";
        }
#       endif

        return;
    }

    HashTbl<T, Key, Hash>* tmpTable = new HashTbl<T, Key, Hash>(newSize);

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        tmpTable->insert(iter.key(), *iter);
    }

    label oldSize = tableSize_;
    tableSize_ = tmpTable->tableSize_;
    tmpTable->tableSize_ = oldSize;

    hashedEntry** oldTable = table_;
    table_ = tmpTable->table_;
    tmpTable->table_ = oldTable;

    delete tmpTable;
}


template<class T, class Key, class Hash>
void Foam::HashTbl<T, Key, Hash>::clear()
{
    if (nElmts_)
    {
        for (label hashIdx = 0; hashIdx < tableSize_; hashIdx++)
        {
            if (table_[hashIdx])
            {
                hashedEntry* ep = table_[hashIdx];
                while (hashedEntry* next = ep->next_)
                {
                    delete ep;
                    ep = next;
                }
                delete ep;
                table_[hashIdx] = 0;
            }
        }
        nElmts_ = 0;
    }
}


template<class T, class Key, class Hash>
void Foam::HashTbl<T, Key, Hash>::clearStorage()
{
    clear();
    resize(0);
}


template<class T, class Key, class Hash>
void Foam::HashTbl<T, Key, Hash>::shrink()
{
    const label newSize = canonicalSize(nElmts_);

    if (newSize < tableSize_)
    {
        // avoid having the table disappear on us
        resize(newSize ? newSize : 2);
    }
}


template<class T, class Key, class Hash>
void Foam::HashTbl<T, Key, Hash>::transfer(HashTbl<T, Key, Hash>& ht)
{
    // as per the Destructor
    if (table_)
    {
        clear();
        delete[] table_;
    }

    tableSize_ = ht.tableSize_;
    ht.tableSize_ = 0;

    table_ = ht.table_;
    ht.table_ = NULL;

    nElmts_ = ht.nElmts_;
    ht.nElmts_ = 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void Foam::HashTbl<T, Key, Hash>::operator=
(
    const HashTbl<T, Key, Hash>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "HashTbl<T, Key, Hash>::operator="
            "(const HashTbl<T, Key, Hash>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    // could be zero-sized from a previous transfer()
    if (!tableSize_)
    {
        resize(rhs.tableSize_);
    }
    else
    {
        clear();
    }

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        insert(iter.key(), *iter);
    }
}


template<class T, class Key, class Hash>
bool Foam::HashTbl<T, Key, Hash>::operator==
(
    const HashTbl<T, Key, Hash>& rhs
) const
{
    // sizes (number of keys) must match
    if (size() != rhs.size())
    {
        return false;
    }

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        const_iterator fnd = find(iter.key());

        if (fnd == cend() || fnd() != iter())
        {
            return false;
        }
    }

    return true;
}


template<class T, class Key, class Hash>
bool Foam::HashTbl<T, Key, Hash>::operator!=
(
    const HashTbl<T, Key, Hash>& rhs
) const
{
    return !(operator==(rhs));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "HashTblIO.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
