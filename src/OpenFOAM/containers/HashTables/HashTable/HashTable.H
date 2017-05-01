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

Class
    Foam::HashTable

Description
    An STL-conforming hash table.

Note
    Hashing index collisions are handled via chaining using a singly-linked
    list with the colliding entry being added to the head of the linked
    list. Thus copying the hash table (or indeed even resizing it) will
    often result in a different hash order. Use a sorted table-of-contents
    when the hash order is important.

SourceFiles
    HashTableI.H
    HashTable.C
    HashTableCore.C
    HashTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashTable_H
#define HashTable_H

#include "label.H"
#include "uLabel.H"
#include "word.H"
#include "Xfer.H"
#include "className.H"
#include "nullObject.H"

#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class List;
template<class T> class UList;
template<class T, unsigned Size> class FixedList;
template<class T1, class T2> class Tuple2;
template<class T, class Key, class Hash> class HashTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream& is, HashTable<T, Key, Hash>& L);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream& os, const HashTable<T, Key, Hash>& tbl);


/*---------------------------------------------------------------------------*\
                        Class HashTableCore Declaration
\*---------------------------------------------------------------------------*/

//- Template-invariant bits for HashTable.
//  This also includes a global end-iterator.
//
// The end iterator of all hash-tables has a nullptr to the hash entry.
// Thus avoid separate allocation for each table and use a single one with
// a nullptr.
// The hash-table iterators always have this as its first member data,
// so we can reinterpret_cast from anything else that has a nullptr for its
// first data member. This is now also the case for the NullObject.
class HashTableCore
{
public:

    //- Return a canonical (power-of-two) size
    static label canonicalSize(const label size);

    //- Maximum allowable table size
    static const label maxTableSize;

    //- Construct null
    HashTableCore()
    {}

    //- Define template name and debug
    ClassName("HashTable");


protected:

    static_assert
    (
        sizeof(NullObject) >= sizeof(void*),
        "NullObject is too small to reinterpret_cast as HashTable::iterator"
    );

    //- Reinterpret a NullObject as a hash-table iterator.
    template<class Iterator>
    inline static const Iterator& endIteratorRef()
    {
        return *reinterpret_cast<const Iterator*>(nullObjectPtr);
    }
};


/*---------------------------------------------------------------------------*\
                           Class HashTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class HashTable
:
    public HashTableCore
{
    // Private data type for table entries

        //- Structure to hold a hashed entry, with a SLList for collisions
        struct hashedEntry
        {
            //- The lookup key
            Key key_;

            //- The data object
            T obj_;

            //- Pointer to next hashedEntry in sub-list
            hashedEntry* next_;

            //- Construct from key, next pointer and object
            inline hashedEntry(const Key& key, const T& obj, hashedEntry* next);


        private:
            //- Disallow default bitwise copy construct
            hashedEntry(const hashedEntry&) = delete;

            //- Disallow default bitwise assignment
            void operator=(const hashedEntry&) = delete;
        };


    // Private data: size of table, the table and current number of elements

        //- The current number of elements in table
        label nElmts_;

        //- Number of primary entries allocated in table
        label tableSize_;

        //- The table of primary entries
        hashedEntry** table_;


    // Private Member Functions

        //- Return the hash index of the Key within the current table size.
        //  No checks for zero-sized tables.
        inline label hashKeyIndex(const Key& key) const;

        //- Assign a new hash-entry to a possibly already existing key.
        //  Return true if the new entry was set.
        bool set(const Key& key, const T& newEntry, const bool protect);


protected:

    // Protected Member Functions

        //- Remove using begin/end iterators of listed keys
        template<class InputIter>
        inline label eraseMultiple
        (
            const InputIter begIter,
            const InputIter endIter
        );


public:

    // Forward declaration of iterators

        class iteratorBase;
        class iterator;
        class const_iterator;

        //- Declare friendship with the iteratorBase
        friend class iteratorBase;


    // Constructors

        //- Construct given initial table size
        HashTable(const label size = 128);

        //- Construct from Istream
        HashTable(Istream& is, const label size = 128);

        //- Construct as copy
        HashTable(const HashTable<T, Key, Hash>& ht);

        //- Construct by transferring the parameter contents
        HashTable(const Xfer<HashTable<T, Key, Hash>>& ht);

        //- Construct from an initializer list
        HashTable(std::initializer_list<Tuple2<Key, T>> lst);


    //- Destructor
    ~HashTable();


    // Member Functions

        // Access

            //- The size of the underlying table
            inline label capacity() const;

            //- Return number of elements in table
            inline label size() const;

            //- Return true if the hash table is empty
            inline bool empty() const;

            //- Return true if hashedEntry is found in table
            bool found(const Key& key) const;

            //- Find and return an iterator set at the hashedEntry
            //  If not found iterator = end()
            iterator find(const Key& key);

            //- Find and return an const_iterator set at the hashedEntry
            //  If not found iterator = end()
            const_iterator find(const Key& key) const;

            //- Return the table of contents
            List<Key> toc() const;

            //- Return the table of contents as a sorted list
            List<Key> sortedToc() const;


        // Edit

            //- Insert a new hashedEntry
            //  Return true if the entry inserted, which means that it did
            //  not previously exist in the table.
            inline bool insert(const Key& key, const T& newEntry);

            //- Assign a new hashedEntry, overwriting existing entries.
            //  Returns true.
            inline bool set(const Key& key, const T& newEntry);

            //- Erase a hashedEntry specified by given iterator
            //  This invalidates the iterator until the next operator++
            bool erase(const iterator& iter);

            //- Erase a hashedEntry specified by the given key
            bool erase(const Key& key);

            //- Remove entries given by the listed keys from this HashTable
            //  Return the number of elements removed
            label erase(const UList<Key>& keys);

            //- Remove entries given by the listed keys from this HashTable
            //  Return the number of elements removed
            template<unsigned Size>
            label erase(const FixedList<Key, Size>& keys);

            //- Remove entries given by the listed keys from this HashTable
            //  Return the number of elements removed
            label erase(std::initializer_list<Key> keys);

            //- Remove entries given by the given keys from this HashTable
            //  Return the number of elements removed.
            //  The parameter HashTable needs the same type of key, but the
            //  type of values held and the hashing function are arbitrary.
            template<class AnyType, class AnyHash>
            label erase(const HashTable<AnyType, Key, AnyHash>& other);

            //- Resize the hash table for efficiency
            void resize(const label sz);

            //- Clear all entries from table
            void clear();

            //- Clear the table entries and the table itself.
            //  Equivalent to clear() followed by resize(0)
            void clearStorage();

            //- Shrink the allocated table to approx. twice number of elements
            void shrink();

            //- Transfer the contents of the argument table into this table
            //  and annul the argument table.
            void transfer(HashTable<T, Key, Hash>& ht);

            //- Transfer contents to the Xfer container
            inline Xfer<HashTable<T, Key, Hash>> xfer();


    // Member Operators

        //- Find and return a hashedEntry
        inline T& operator[](const Key& key);

        //- Find and return a hashedEntry
        inline const T& operator[](const Key& key) const;

        //- Return existing entry or create a new entry.
        //  A newly created entry is created as a nameless T() and is thus
        //  value-initialized. For primitives, this will be zero.
        inline T& operator()(const Key& key);

        //- Assignment
        void operator=(const HashTable<T, Key, Hash>& rhs);

        //- Assignment from an initializer list
        void operator=(std::initializer_list<Tuple2<Key, T>> lst);

        //- Equality. Hash tables are equal if the keys and values are equal.
        //  Independent of table storage size and table order.
        bool operator==(const HashTable<T, Key, Hash>& rhs) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const HashTable<T, Key, Hash>& rhs) const;


    // STL type definitions

        //- Type of keys that the HashTable uses.
        typedef Key key_type;

        //- Type of values that the HashTable contains.
        typedef T value_type;

        //- Type that can be used for storing into HashTable::value_type
        //  objects.  This type is usually List::value_type&.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  HashTable::value_type objects.  This type is usually const
        //  HashTable::value_type&.
        typedef const T& const_reference;

        //- The type that can represent the size of a HashTable.
        typedef label size_type;


    // Iterator access

        //- Iterator set to the beginning of the HashTable
        inline iterator begin();

        //- const_iterator set to the beginning of the HashTable
        inline const_iterator begin() const;

        //- const_iterator set to the beginning of the HashTable
        inline const_iterator cbegin() const;

        //- iterator to signal the end for any HashTable
        inline const iterator& end();

        //- const_iterator to signal the end for any HashTable
        inline const const_iterator& end() const;

        //- const_iterator to signal the end for any HashTable
        inline const const_iterator& cend() const;


    // Iterators and helpers

        //- The iterator base for HashTable
        //  Note: data and functions are protected, to allow reuse by iterator
        //  and prevent most external usage.
        //  iterator and const_iterator have the same size, allowing
        //  us to reinterpret_cast between them (if desired)
        class iteratorBase
        {
            using entry_type = hashedEntry;

        public:
            // Public typedefs
            using table_type = HashTable<T, Key, Hash>;
            using key_type   = HashTable<T, Key, Hash>::key_type;

        private:

          // Private Data

            //- Currently selected entry.
            //  MUST be the first member for easy comparison between iterators
            //  and for reinterpret_cast from nullObject
            entry_type* entryPtr_;

            //- Pointer to the hash-table for which this is an iterator
            //  This allows use of the default bitwise copy/assignment
            table_type* hashTable_;

            //- Current hash index within the hash-table data.
            //  A signed value, since erase() uses a negative value to signal
            //  the erasure state.
            label hashIndex_;

        protected:

          // Protected Member Functions

            //- Increment to the next position
            inline void increment();

            //- Erase the entry at the current position
            bool erase();

            //- The referenced object/value element
            inline T& element() const;


        public:

          // Constructors

            //- Construct null (end iterator)
            inline iteratorBase();

            //- Construct from begin of hash-table
            inline explicit iteratorBase(const table_type* hashTbl);

            //- Construct from hash table, element and hash index
            inline iteratorBase
            (
                const table_type* hashTbl,
                const entry_type* elmt,
                const label hashIndex
            );

          // Member functions/operators

            //- True if iterator points to an entry
            //  This can be used directly instead of comparing to end()
            inline bool found() const;

            //- Return the Key corresponding to the iterator
            inline const Key& key() const;

            //- Compare hash-entry element pointers
            inline bool operator==(const iteratorBase& iter) const;
            inline bool operator!=(const iteratorBase& iter) const;
        };


    // STL iterator

        //- An STL-conforming iterator
        class iterator
        :
            public iteratorBase
        {
            friend class HashTable;  // uses iterator::erase() method
            using entry_type = hashedEntry;

        public:

          // Public typedefs
            using table_type = HashTable<T, Key, Hash>;
            using key_type   = HashTable<T, Key, Hash>::key_type;
            using reference  = HashTable<T, Key, Hash>::reference;

          // Constructors

            //- Construct null (end iterator)
            inline iterator();

            //- Construct from begin of hash-table
            inline explicit iterator(table_type* hashTbl);

            //- Construct from hash table, element and hash index
            //  Used by the hash-table find() method.
            inline iterator
            (
                table_type* hashTbl,
                entry_type* elmt,
                const label hashIndex
            );

          // Member functions/operators

            //- Return non-const access to referenced object
            inline reference object() const;

            //- Return non-const access to referenced object
            inline reference operator*() const;
            inline reference operator()() const;

            inline iterator& operator++();
            inline iterator operator++(int);
        };


    // STL const_iterator

        //- An STL-conforming const_iterator
        class const_iterator
        :
            public iteratorBase
        {
            using entry_type = const hashedEntry;

        public:

          // Public typedefs
            using table_type = const HashTable<T, Key, Hash>;
            using key_type   = HashTable<T, Key, Hash>::key_type;
            using reference  = HashTable<T, Key, Hash>::const_reference;

          // Constructors

            //- Construct null (end iterator)
            inline const_iterator();

            //- Construct from begin of hash-table
            inline explicit const_iterator(table_type* hashTbl);

            //- Construct from hash table, element and hash index.
            //  Used by the hash-table find() method.
            inline const_iterator
            (
                table_type* hashTbl,
                entry_type* elmt,
                const label hashIndex
            );

            //- Copy construct from iterator
            inline const_iterator(const iterator& iter);

          // Member functions/operators

            //- Return const access to referenced object
            inline reference object() const;

            //- Return const access to referenced object
            inline reference operator*() const;
            inline reference operator()() const;

            inline const_iterator& operator++();
            inline const_iterator operator++(int);
        };


    // Writing

        //- Print information
        Ostream& printInfo(Ostream& os) const;

        //- Write the unordered keys as a list, with line-breaks if list length
        //  exceeds shortListLen. Using '0' suppresses line-breaks entirely.
        Ostream& writeKeys(Ostream& os, const label shortListLen=0) const;


    // IOstream Operator

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream& is,
            HashTable<T, Key, Hash>& L
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream& os,
            const HashTable<T, Key, Hash>& tbl
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HashTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef NoHashTableC
#ifdef NoRepository
    #include "HashTable.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
