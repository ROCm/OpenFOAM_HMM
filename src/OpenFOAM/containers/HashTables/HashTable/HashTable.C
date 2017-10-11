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

#ifndef HashTable_C
#define HashTable_C

#include "HashTable.H"
#include "List.H"
#include "FixedList.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class InputIter>
Foam::label Foam::HashTable<T, Key, Hash>::eraseMultiple
(
    const InputIter begIter,
    const InputIter endIter
)
{
    const label nTotal = this->size();
    label changed = 0;

    for
    (
        InputIter iter = begIter;
        changed < nTotal && iter != endIter; // terminate early
        ++iter
    )
    {
        if (this->erase(*iter))
        {
            ++changed;
        }
    }
    return changed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(const label size)
:
    HashTableCore(),
    nElmts_(0),
    tableSize_(HashTableCore::canonicalSize(size)),
    table_(nullptr)
{
    if (tableSize_)
    {
        table_ = new hashedEntry*[tableSize_];

        for (label hashIdx = 0; hashIdx < tableSize_; ++hashIdx)
        {
            table_[hashIdx] = nullptr;
        }
    }
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(const HashTable<T, Key, Hash>& ht)
:
    HashTable<T, Key, Hash>(ht.tableSize_)
{
    for (const_iterator iter = ht.cbegin(); iter != ht.cend(); ++iter)
    {
        insert(iter.key(), iter.object());
    }
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(HashTable<T, Key, Hash>&& ht)
:
    HashTableCore(),
    nElmts_(0),
    tableSize_(0),
    table_(nullptr)
{
    transfer(ht);
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable
(
    const Xfer<HashTable<T, Key, Hash>>& ht
)
:
    HashTableCore(),
    nElmts_(0),
    tableSize_(0),
    table_(nullptr)
{
    transfer(ht());
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable
(
    std::initializer_list<std::pair<Key, T>> lst
)
:
    HashTable<T, Key, Hash>(2*lst.size())
{
    for (const auto& pair : lst)
    {
        insert(pair.first, pair.second);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::~HashTable()
{
    if (table_)
    {
        clear();
        delete[] table_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::found(const Key& key) const
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

    #ifdef FULLDEBUG
    if (debug)
    {
        InfoInFunction << "Entry " << key << " not found in hash table\n";
    }
    #endif

    return false;
}


template<class T, class Key, class Hash>
typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::find
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

    #ifdef FULLDEBUG
    if (debug)
    {
        InfoInFunction << "Entry " << key << " not found in hash table\n";
    }
    #endif

    return iterator();
}


template<class T, class Key, class Hash>
typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::find
(
    const Key& key
) const
{
    return this->cfind(key);
}


template<class T, class Key, class Hash>
typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::cfind
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

    #ifdef FULLDEBUG
    if (debug)
    {
        InfoInFunction << "Entry " << key << " not found in hash table\n";
    }
    #endif

    return const_iterator();
}


template<class T, class Key, class Hash>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::toc() const
{
    List<Key> keyLst(nElmts_);
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        keyLst[count++] = iter.key();
    }

    return keyLst;
}


template<class T, class Key, class Hash>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::sortedToc() const
{
    List<Key> keyLst = this->toc();
    Foam::sort(keyLst);

    return keyLst;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::tocKeys
(
    const UnaryPredicate& pred,
    const bool invert
) const
{
    List<Key> keyLst(nElmts_);
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.key()) ? !invert : invert))
        {
            keyLst[count++] = iter.key();
        }
    }

    keyLst.setSize(count);
    Foam::sort(keyLst);

    return keyLst;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::tocValues
(
    const UnaryPredicate& pred,
    const bool invert
) const
{
    List<Key> keyLst(nElmts_);
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.object()) ? !invert : invert))
        {
            keyLst[count++] = iter.key();
        }
    }

    keyLst.setSize(count);
    Foam::sort(keyLst);

    return keyLst;
}


template<class T, class Key, class Hash>
template<class BinaryPredicate>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::tocEntries
(
    const BinaryPredicate& pred,
    const bool invert
) const
{
    List<Key> keyLst(nElmts_);
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.key(), iter.object()) ? !invert : invert))
        {
            keyLst[count++] = iter.key();
        }
    }

    keyLst.setSize(count);
    Foam::sort(keyLst);

    return keyLst;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::countKeys
(
    const UnaryPredicate& pred,
    const bool invert
) const
{
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.key()) ? !invert : invert))
        {
            ++count;
        }
    }

    return count;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::countValues
(
    const UnaryPredicate& pred,
    const bool invert
) const
{
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.object()) ? !invert : invert))
        {
            ++count;
        }
    }

    return count;
}


template<class T, class Key, class Hash>
template<class BinaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::countEntries
(
    const BinaryPredicate& pred,
    const bool invert
) const
{
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.key(), iter.object()) ? !invert : invert))
        {
            ++count;
        }
    }

    return count;
}


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& obj,
    const bool protect
)
{
    if (!tableSize_)
    {
        resize(2);
    }

    const label hashIdx = hashKeyIndex(key);

    hashedEntry* existing = nullptr;
    hashedEntry* prev = nullptr;

    for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
    {
        if (key == ep->key_)
        {
            existing = ep;
            break;
        }
        prev = ep;
    }

    if (!existing)
    {
        // Not found, insert it at the head
        table_[hashIdx] = new hashedEntry(key, obj, table_[hashIdx]);
        nElmts_++;

        if (double(nElmts_)/tableSize_ > 0.8 && tableSize_ < maxTableSize)
        {
            #ifdef FULLDEBUG
            if (debug)
            {
                InfoInFunction << "Doubling table size\n";
            }
            #endif

            resize(2*tableSize_);
        }
    }
    else if (protect)
    {
        // Found - but protected from overwriting
        // this corresponds to the STL 'insert' convention
        #ifdef FULLDEBUG
        if (debug)
        {
            InfoInFunction
                << "Cannot insert " << key << " already in hash table\n";
        }
        #endif
        return false;
    }
    else
    {
        // Found - overwrite existing entry
        // this corresponds to the Perl convention
        hashedEntry* ep = new hashedEntry(key, obj, existing->next_);

        // Replace existing element - within list or insert at the head
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
bool Foam::HashTable<T, Key, Hash>::iterator_base::erase()
{
    // Note: entryPtr_ is nullptr for end(), so this catches that too
    if (entryPtr_)
    {
        // Search element before entryPtr_
        entry_type* prev = nullptr;

        for
        (
            entry_type* ep = hashTable_->table_[hashIndex_];
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
            // Has an element before entryPtr - reposition to there
            prev->next_ = entryPtr_->next_;
            delete entryPtr_;
            entryPtr_ = prev;
        }
        else
        {
            // entryPtr was first element on SLList
            hashTable_->table_[hashIndex_] = entryPtr_->next_;
            delete entryPtr_;

            // Assign any non-nullptr value so it doesn't look like end()
            entryPtr_ = reinterpret_cast<hashedEntry*>(this);

            // Mark with special hashIndex value to signal it has been rewound.
            // The next increment will bring it back to the present location.
            //
            // From the current position 'curPos', we wish to continue at
            // prevPos='curPos-1', which we mark as markPos='-curPos-1'.
            // The negative lets us notice it is special, the extra '-1'
            // is needed to avoid ambiguity for position '0'.
            // To retrieve prevPos, we would later use '-(markPos+1) - 1'
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


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::erase(const iterator& iter)
{
    // NOTE: We use (const iterator&) here, but manipulate its contents anyhow.
    // The parameter should be (iterator&), but then the compiler doesn't find
    // it correctly and tries to call as (iterator) instead.
    //
    // Adjust iterator after erase
    return const_cast<iterator&>(iter).erase();
}


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::erase(const Key& key)
{
    auto iter = find(key);
    return erase(iter);
}


template<class T, class Key, class Hash>
Foam::label Foam::HashTable<T, Key, Hash>::erase(const UList<Key>& keys)
{
    return eraseMultiple(keys.begin(), keys.end());
}


template<class T, class Key, class Hash>
template<unsigned Size>
Foam::label Foam::HashTable<T, Key, Hash>::erase
(
    const FixedList<Key, Size>& keys
)
{
    return eraseMultiple(keys.begin(), keys.end());
}


template<class T, class Key, class Hash>
Foam::label Foam::HashTable<T, Key, Hash>::erase
(
    std::initializer_list<Key> keys
)
{
    return eraseMultiple(keys.begin(), keys.end());
}


template<class T, class Key, class Hash>
template<class AnyType, class AnyHash>
Foam::label Foam::HashTable<T, Key, Hash>::erase
(
    const HashTable<AnyType, Key, AnyHash>& other
)
{
    const label nTotal = this->size();
    label changed = 0;

    using other_iter =
        typename HashTable<AnyType, Key, AnyHash>::const_iterator;

    if (other.size() <= nTotal)
    {
        // The other is smaller/same-size, use its keys for removal

        for
        (
            other_iter iter = other.begin();
            changed < nTotal && iter != other.end(); // terminate early
            ++iter
        )
        {
            if (erase(iter.key()))
            {
                ++changed;
            }
        }
    }
    else
    {
        // We are smaller: remove if found in the other hash
        for
        (
            iterator iter = begin();
            changed < nTotal && iter != end(); // terminate early
            ++iter
        )
        {
            if (other.found(iter.key()) && erase(iter))
            {
                ++changed;
            }
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
template<class AnyType, class AnyHash>
Foam::label Foam::HashTable<T, Key, Hash>::retain
(
    const HashTable<AnyType, Key, AnyHash>& other
)
{
    const label nTotal = this->size();
    label changed = 0;

    if (other.empty())
    {
        // Trivial case
        changed = nTotal;
        this->clear();
    }
    else
    {
        // Inverted logic: remove if *not* found in the other hash

        for (iterator iter = begin(); iter != end(); ++iter)
        {
            if (!other.found(iter.key()) && erase(iter))
            {
                ++changed;
            }
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::resize(const label sz)
{
    const label newSize = HashTableCore::canonicalSize(sz);

    if (newSize == tableSize_)
    {
        #ifdef FULLDEBUG
        if (debug)
        {
            InfoInFunction << "New table size == old table size\n";
        }
        #endif

        return;
    }

    HashTable<T, Key, Hash>* tmpTable = new HashTable<T, Key, Hash>(newSize);

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        tmpTable->insert(iter.key(), iter.object());
    }

    const label oldSize = tableSize_;
    tableSize_ = tmpTable->tableSize_;
    tmpTable->tableSize_ = oldSize;

    hashedEntry** oldTable = table_;
    table_ = tmpTable->table_;
    tmpTable->table_ = oldTable;

    delete tmpTable;
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::clear()
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
                table_[hashIdx] = nullptr;
            }
        }
        nElmts_ = 0;
    }
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::clearStorage()
{
    clear();
    resize(0);
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::swap(HashTable<T, Key, Hash>& ht)
{
    Foam::Swap(table_,     ht.table_);
    Foam::Swap(tableSize_, ht.tableSize_);
    Foam::Swap(nElmts_,    ht.nElmts_);
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::transfer(HashTable<T, Key, Hash>& ht)
{
    // As per the Destructor
    if (table_)
    {
        clear();
        delete[] table_;
    }

    tableSize_ = ht.tableSize_;
    ht.tableSize_ = 0;

    table_ = ht.table_;
    ht.table_ = nullptr;

    nElmts_ = ht.nElmts_;
    ht.nElmts_ = 0;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::filterKeys
(
    const UnaryPredicate& pred,
    const bool pruning
)
{
    label changed = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        // Matches? either prune (pruning) or keep (!pruning)
        if
        (
            (pred(iter.key()) ? pruning : !pruning)
         && erase(iter)
        )
        {
            ++changed;
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::filterValues
(
    const UnaryPredicate& pred,
    const bool pruning
)
{
    label changed = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        // Matches? either prune (pruning) or keep (!pruning)
        if
        (
            (pred(iter.object()) ? pruning : !pruning)
         && erase(iter)
        )
        {
            ++changed;
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
template<class BinaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::filterEntries
(
    const BinaryPredicate& pred,
    const bool pruning
)
{
    label changed = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        // Matches? either prune (pruning) or keep (!pruning)
        if
        (
            (pred(iter.key(), iter.object()) ? pruning : !pruning)
         && erase(iter)
        )
        {
            ++changed;
        }
    }

    return changed;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::operator=
(
    const HashTable<T, Key, Hash>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    // Could be zero-sized from a previous transfer()
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
        insert(iter.key(), iter.object());
    }
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::operator=
(
    std::initializer_list<std::pair<Key, T>> lst
)
{
    // Could be zero-sized from a previous transfer()
    if (!tableSize_)
    {
        resize(2*lst.size());
    }
    else
    {
        clear();
    }

    for (const auto& pair : lst)
    {
        insert(pair.first, pair.second);
    }
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::operator=
(
    HashTable<T, Key, Hash>&& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    transfer(rhs);
}


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::operator==
(
    const HashTable<T, Key, Hash>& rhs
) const
{
    // Sizes (number of keys) must match
    if (size() != rhs.size())
    {
        return false;
    }

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        const_iterator other = this->cfind(iter.key());

        if (!other.found() || other.object() != iter.object())
        {
            return false;
        }
    }

    return true;
}


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::operator!=
(
    const HashTable<T, Key, Hash>& rhs
) const
{
    return !operator==(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "HashTableIO.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
