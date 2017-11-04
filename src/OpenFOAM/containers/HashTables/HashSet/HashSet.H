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
        HashSet<label> someLabels{10, 20, 30, 40, ...};
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
    A HashSet with label keys.

\*---------------------------------------------------------------------------*/

#ifndef HashSet_H
#define HashSet_H

#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
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

        //- Insert values, using begin/end iterators.
        template<class InputIter>
        inline label insertMultiple
        (
            const InputIter begIter,
            const InputIter endIter
        );

        //- Assign using begin/end iterators.
        template<class InputIter>
        inline label assignMultiple
        (
            const InputIter begIter,
            const InputIter endIter,
            const label sz
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

        //- Construct from UList of Key
        explicit HashSet(const UList<Key>& lst);

        //- Construct from FixedList of Key
        template<unsigned Size>
        explicit HashSet(const FixedList<Key, Size>& lst);

        //- Construct from an initializer list of Key
        HashSet(std::initializer_list<Key> lst);

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

        //- Construct by transferring the parameter contents
        HashSet(const Xfer<this_type>& hs)
        :
            parent_type(hs)
        {}

        //- Construct by transferring the parameter contents
        HashSet(const Xfer<parent_type>& hs)
        :
            parent_type(hs)
        {}

        //- Construct from the keys of another HashTable,
        //  the type of values held is arbitrary.
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

        //- Insert keys from the list of Key
        //  \return The number of new elements inserted
        label insert(const UList<Key>& lst);

        //- Insert keys from the list of Key
        //  \return The number of new elements inserted
        template<unsigned Size>
        label insert(const FixedList<Key, Size>& lst);

        //- Insert keys from a initializer list of Key
        //  \return The number of new elements inserted
        label insert(std::initializer_list<Key> lst);

        //- Same as insert (no value to overwrite)
        bool set(const Key& key)
        {
            return insert(key);
        }

        //- Same as insert (no value to overwrite)
        label set(const UList<Key>& lst)
        {
            return insert(lst);
        }

        //- Same as insert (no value to overwrite)
        template<unsigned Size>
        label set(const FixedList<Key, Size>& lst)
        {
            return insert(lst);
        }

        //- Same as insert (no value to overwrite)
        label set(std::initializer_list<Key> lst)
        {
            return insert(lst);
        }

        //- Unset the specified key - same as erase
        //  \return True if the entry existed and was removed
        bool unset(const Key& key)
        {
            return this->parent_type::erase(key);
        }

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        label unset(const UList<Key>& lst)
        {
            return this->parent_type::erase(lst);
        }

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        template<unsigned Size>
        label unset(const FixedList<Key, Size>& lst)
        {
            return this->parent_type::erase(lst);
        }

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        label unset(std::initializer_list<Key> lst)
        {
            return this->parent_type::erase(lst);
        }

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
        void operator=(const UList<Key>& lst);

        //- Assignment from a FixedList of keys
        template<unsigned Size>
        void operator=(const FixedList<Key, Size>& lst);

        //- Assignment from an initializer list of keys
        void operator=(std::initializer_list<Key> lst);


      // Logical operations

        //- Combine entries from HashSets
        void operator|=(const this_type& rhs);

        //- Only retain entries found in both HashSets
        inline void operator&=(const this_type& rhs);

        //- Only retain unique entries (xor)
        void operator^=(const this_type& rhs);

        //- Add entries listed in the given HashSet to this HashSet
        inline void operator+=(const this_type& rhs)
        {
            this->operator|=(rhs);
        }

        //- Remove entries listed in the given HashSet from this HashSet
        inline void operator-=(const this_type& rhs);


    // IOstream Operator

        friend Ostream& operator<< <Key, Hash>
        (
            Ostream& os,
            const HashSet<Key, Hash>& tbl
        );

};


// Global Functions

// Exchange contents of HashSets - see HashSet::swap().
template<class Key, class Hash>
inline void Swap
(
    HashSet<Key, Hash>& a,
    HashSet<Key, Hash>& b
);


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
typedef HashSet<> wordHashSet;

//- A HashSet with label keys.
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
