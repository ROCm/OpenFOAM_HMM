/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

Class
    Foam::HashSet

Description
    A HashTable with keys but without contents that is similar to
    \c std::unordered_set.

    The entries are considered \a unordered since their placement
    depends on the method used to generate the hash key index, the
    table capacity, insertion order etc. When the key order is
    important, use the sortedToc() method to obtain a list of sorted
    keys and use that for further access.

Note
    The HashSet iterator dereferences to the key, so the following
    range-for works as expected:
    \code
        labelHashSet someLabels{10, 20, 30, 40, ...};
        for (const label i : someLabels)
        {
            Info<< "val:" << i << nl;
        }
    \endcode

Typedef
    Foam::wordHashSet

Description
    A HashSet with (the default) word keys.

Typedef
    Foam::labelHashSet

Description
    A HashSet with label keys and label hasher.

\*---------------------------------------------------------------------------*/

#ifndef HashSet_H
#define HashSet_H

#include "HashTable.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class Key, class Hash> class HashSet;

template<class Key, class Hash>
Ostream& operator<<(Ostream& os, const HashSet<Key, Hash>& tbl);


/*---------------------------------------------------------------------------*\
                           Class HashSet Declaration
\*---------------------------------------------------------------------------*/

template<class Key=word, class Hash=string::hash>
class HashSet
:
    public HashTable<zero::null, Key, Hash>
{
    // Private Member Functions

        //- Assign from the input iterator range
        template<class InputIter>
        inline label assignMany
        (
            const label nItems,
            InputIter first,
            InputIter last
        );


public:

    //- The template instance used for this HashSet
    typedef HashSet<Key, Hash> this_type;

    //- The template instance used for the parent HashTable
    typedef HashTable<zero::null, Key, Hash> parent_type;

    //- An iterator, returning reference to the key
    using iterator = typename parent_type::key_iterator;

    //- A const_iterator, returning reference to the key
    using const_iterator = typename parent_type::const_key_iterator;


    // Constructors

        //- Construct null with default (128) table capacity
        HashSet()
        :
            parent_type()
        {}

        //- Construct given initial table capacity
        explicit HashSet(const label size)
        :
            parent_type(size)
        {}

        //- Construct from Istream with default table capacity
        HashSet(Istream& is)
        :
            parent_type(is)
        {}

        //- Construct from FixedList of Key
        template<unsigned Size>
        explicit HashSet(const FixedList<Key, Size>& list);

        //- Construct from UList of Key
        explicit HashSet(const UList<Key>& list);

        //- Construct from an indirect list
        explicit HashSet(const UIndirectList<Key>& list);

        //- Construct from an initializer list of Key
        HashSet(std::initializer_list<Key> list);

        //- Copy construct
        HashSet(const this_type& hs)
        :
            parent_type(hs)
        {}

        //- Move construct
        HashSet(this_type&& hs)
        :
            parent_type(std::move(hs))
        {}

        //- Construct from the keys of another HashTable,
        //- the type of values held is arbitrary.
        template<class AnyType, class AnyHash>
        explicit HashSet(const HashTable<AnyType, Key, AnyHash>& tbl);


    // Member Functions

    // Edit

        //- Insert a new entry, not overwriting existing entries.
        //  \return True if the entry inserted, which means that it did
        //  not previously exist in the set.
        bool insert(const Key& key)
        {
            return this->parent_type::insert(key, zero::null());
        }

        //- Same as insert (no value to overwrite)
        bool set(const Key& key)
        {
            return insert(key);
        }

        //- Unset the specified key - same as erase
        //  \return True if the entry existed and was removed
        bool unset(const Key& key)
        {
            return this->parent_type::erase(key);
        }


    // Convenience

        //- Insert keys from the input iterator range
        //  \return The number of new elements inserted
        template<class InputIter>
        inline label insert(InputIter first, InputIter last);

        //- Insert keys from a initializer list of Key
        //  \return The number of new elements inserted
        inline label insert(std::initializer_list<Key> list);

        //- Insert keys from the list of Key
        //  \return The number of new elements inserted
        template<unsigned Size>
        inline label insert(const FixedList<Key, Size>& list);

        //- Insert keys from the list of Key
        //  \return The number of new elements inserted
        inline label insert(const UList<Key>& list);

        //- Insert keys from the list of Key
        //  \return The number of new elements inserted
        inline label insert(const UIndirectList<Key>& list);

        //- Same as insert (no value to overwrite)
        template<class InputIter>
        inline label set(InputIter first, InputIter last)
        {
            return insert(first, last);
        }

        //- Same as insert (no value to overwrite)
        inline label set(std::initializer_list<Key> list)
        {
            return insert(list);
        }

        //- Same as insert (no value to overwrite)
        template<unsigned Size>
        inline label set(const FixedList<Key, Size>& list)
        {
            return insert(list);
        }

        //- Same as insert (no value to overwrite)
        inline label set(const UList<Key>& list)
        {
            return insert(list);
        }

        //- Same as insert (no value to overwrite)
        inline label set(const UIndirectList<Key>& list)
        {
            return insert(list);
        }

        //- Same as insert (no value to overwrite)
        //  \note Method name compatibility with bitSet
        template<class InputIter>
        inline label setMany(InputIter first, InputIter last)
        {
            return insert(first, last);
        }

        //- Unset the keys listed in the input iterator range
        //  \return The number of items removed
        template<class InputIter>
        inline label unset(InputIter first, InputIter last);

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        inline label unset(std::initializer_list<Key> list);

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        template<unsigned Size>
        inline label unset(const FixedList<Key, Size>& list);

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        inline label unset(const UList<Key>& list);

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        inline label unset(const UIndirectList<Key>& list);


    // STL iterators

        iterator begin();
        const_iterator begin() const;
        const_iterator cbegin() const;

        const iterator& end();
        const const_iterator& end() const;
        const const_iterator& cend() const;


    // Writing

        //- Write the unordered keys as a list, with line-breaks if list length
        //- exceeds shortListLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortListLen=0) const
        {
            return this->writeKeys(os, shortListLen);
        }


    // Member Operators

        //- Return true if the entry exists, same as found()
        inline bool operator()(const Key& key) const;

        //- Return true if the entry exists, same as found().
        inline bool operator[](const Key& key) const;

        using parent_type::operator=;

        //- Copy assignment
        void operator=(const this_type& rhs)
        {
            parent_type::operator=(rhs);
        }

        //- Move assignment
        void operator=(this_type&& rhs)
        {
            parent_type::operator=(std::move(rhs));
        }


    // Comparison

        //- Sets are equal if all keys are equal,
        //- independent of order or underlying storage size.
        bool operator==(const this_type& rhs) const;

        //- The opposite of the equality operation.
        bool operator!=(const this_type& rhs) const;


    // Assignment

        //- Assignment from a UList of keys
        void operator=(const UList<Key>& rhs);

        //- Assignment from a FixedList of keys
        template<unsigned Size>
        void operator=(const FixedList<Key, Size>& rhs);

        //- Assignment from an initializer list of keys
        void operator=(std::initializer_list<Key> rhs);


    // Logical and set operations

        //- Add entries to this HashSet
        this_type& operator|=(const this_type& rhs);

        //- Only retain entries found in both HashSets
        inline this_type& operator&=(const this_type& rhs);

        //- Only retain unique entries (xor)
        this_type& operator^=(const this_type& rhs);

        //- Add entries to this HashSet
        inline this_type& operator+=(const this_type& rhs)
        {
            return this->operator|=(rhs);
        }

        //- Remove entries from this HashSet
        inline this_type& operator-=(const this_type& rhs);


    // IOstream Operator

        friend Ostream& operator<< <Key, Hash>
        (
            Ostream& os,
            const HashSet<Key, Hash>& tbl
        );


    // Housekeeping

        //- Not applicable for HashSet
        template<class UnaryPredicate>
        List<Key> tocValues(const UnaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class BinaryPredicate>
        List<Key> tocEntries(const BinaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class UnaryPredicate>
        label countValues(const UnaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class BinaryPredicate>
        label countEntries(const BinaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class UnaryPredicate>
        label filterValues(const UnaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class BinaryPredicate>
        label filterEntries(const BinaryPredicate&, const bool) = delete;

};


// Global Operators

//- Combine entries from HashSets
template<class Key, class Hash>
HashSet<Key,Hash> operator|
(
    const HashSet<Key,Hash>& hash1,
    const HashSet<Key,Hash>& hash2
);


//- Create a HashSet that only contains entries found in both HashSets
template<class Key, class Hash>
HashSet<Key,Hash> operator&
(
    const HashSet<Key,Hash>& hash1,
    const HashSet<Key,Hash>& hash2
);


//- Create a HashSet that only contains unique entries (xor)
template<class Key, class Hash>
HashSet<Key,Hash> operator^
(
    const HashSet<Key,Hash>& hash1,
    const HashSet<Key,Hash>& hash2
);


//- A HashSet with word keys.
typedef HashSet<word> wordHashSet;

//- A HashSet with label keys and label hasher.
typedef HashSet<label, Hash<label>> labelHashSet;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HashSet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
