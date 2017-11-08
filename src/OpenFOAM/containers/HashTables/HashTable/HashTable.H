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
    A HashTable similar to \c std::unordered_map.

    The entries are considered \a unordered since their placement
    depends on the method used to generate the hash key index, the
    table capacity, insertion order etc. When the key order is
    important, use the sortedToc() method to obtain a list of sorted
    keys and use that for further access.

    Internally the table uses closed addressing into a flat storage space
    with collisions handled by linked-list chaining.

    The end iterator of all hash-tables has a nullptr to the hash entry.
    Thus avoid separate allocation for each table and use a single one with
    a nullptr. The hash-table iterators always have an entry-pointer as the
    first member data, which allows reinterpret_cast from anything else with
    a nullptr as its first data member.
    The nullObject is such an item (with a nullptr data member).

Note
    For historical reasons, dereferencing the table iterator
    (eg, \a *iter) returns a reference to the stored object value
    value rather than the stored key/object value like std::unordered_map
    does.

    The HashTable iterator:
    \code
        forAllConstIters(table, iter)
        {
            Info<< "val:" << *iter << nl
                << "key:" << iter.key() << nl;
                << "val:" << iter.object() << nl;
        }
    \endcode
    whereas for the \c std::unordered_map iterator:
    \code
        forAllConstIters(stdmap, iter)
        {
            Info<< "key/val:" << *iter << nl
                << "key:" << iter->first << nl
                << "val:" << iter->second << nl;
        }
    \endcode
    This difference is most evident when using range-for syntax.

SourceFiles
    HashTableI.H
    HashTableIterI.H
    HashTable.C
    HashTableIO.C
    HashTableIter.C

\*---------------------------------------------------------------------------*/

#ifndef HashTable_H
#define HashTable_H

#include "word.H"
#include "zero.H"
#include "Xfer.H"
#include "Hash.H"
#include "HashTableCore.H"

#include <initializer_list>
#include <iterator>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class List;
template<class T> class UList;
template<class T, unsigned Size> class FixedList;
template<class T, class Key, class Hash> class HashTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream& is, HashTable<T, Key, Hash>& L);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream& os, const HashTable<T, Key, Hash>& tbl);


/*---------------------------------------------------------------------------*\
                          Class HashTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class HashTable
:
    public HashTableCore
{
    // Private types for table entries

        //- Structure with a (K,V) tuple and a linked-list for collisions
        //  Could store key/object as std::pair, but no particular advantage
        //  unless the iterator dereference type changes.
        struct pair_entry
        {
            //- Type of key
            typedef Key key_type;

            //- Object content type
            typedef T mapped_type;

            //- The lookup key
            key_type key_;

            //- The data object
            mapped_type obj_;

            //- Addressing (next in collision list)
            pair_entry* next_;

            //- Construct from key, object, next pointer
            pair_entry(const Key& key, const T& obj, pair_entry* next)
            :
                key_(key),
                obj_(obj),
                next_(next)
            {}

            //- The key
            const key_type& key() const
            {
                return key_;
            }

            //- The mapped object
            const mapped_type& mapped() const
            {
                return obj_;
            }
            mapped_type& mapped()
            {
                return obj_;
            }

        private:
            //- Disallow default bitwise copy construct / assignment
            pair_entry(const pair_entry&) = delete;
            void operator=(const pair_entry&) = delete;
        };


        //- Structure with a single (K) value and a linked-list for collisions
        struct unary_entry
        {
            //- Type of key
            typedef Key key_type;

            //- Object content type
            typedef zero::null mapped_type;

            //- Content storage type to the entry
            typedef key_type value_type;

            //- The lookup key == content
            key_type key_;

            //- Addressing (next in collision list)
            unary_entry* next_;

            //- Construct from key, (ununsed) object, next pointer
            unary_entry(const Key& key, const T&, unary_entry* next)
            :
                key_(key),
                next_(next)
            {}

            //- The key
            const key_type& key() const
            {
                return key_;
            }

            //- Dummy mapped object
            const mapped_type& mapped() const
            {
                return zeroNullElement;
            }
            mapped_type& mapped()
            {
                return zeroNullElement;
            }


        private:
            //- Disallow default bitwise copy construct / assignment
            unary_entry(const unary_entry&) = delete;
            void operator=(const unary_entry&) = delete;
        };


        //- Hashed node with a linked-list for collisions
        typedef typename std::conditional
        <
            std::is_same<zero::null, typename std::remove_cv<T>::type>::value,
            unary_entry,
            pair_entry
        >::type node_type;


    // Private Data

        //- The number of nodes currently stored in table
        label size_;

        //- Number of nodes allocated in table
        label capacity_;

        //- The table of primary nodes
        node_type** table_;


    // Private Member Functions

        //- Return the hash index of the Key within the current table size.
        //  No checks for zero-sized tables.
        inline label hashKeyIndex(const Key& key) const;

        //- Assign a new hash-entry to a possibly already existing key.
        //  \return True if the new entry was set.
        bool set(const Key& key, const T& obj, const bool overwrite);


public:

    //- The template instance used for this HashTable
    typedef HashTable<T, Key, Hash> this_type;


    // STL type definitions

        //- The second template parameter, type of keys used.
        typedef Key key_type;

        //- The first template parameter, type of objects contained.
        typedef T mapped_type;

        //- Same as mapped_type for OpenFOAM HashTables
        //  Note that this is different than the std::map definition.
        typedef T value_type;

        //- The third template parameter, the hash index method.
        typedef Hash hasher;

        //- The type used for storing into value_type objects.
        //  This type is usually 'value_type*'.
        typedef T* pointer;

        //- The type used for storing into value_type objects.
        //  This type is usually 'value_type&'.
        typedef T& reference;

        //- The type used for reading from constant value_type objects.
        typedef const T* const_pointer;

        //- The type used for reading from constant value_type objects.
        typedef const T& const_reference;

        //- The type to represent the difference between two iterators
        typedef label difference_type;

        //- The type that can represent the size of a HashTable.
        typedef label size_type;


        //- Forward iterator with non-const access
        class iterator;

        //- Forward iterator with const access
        class const_iterator;

protected:

        //- Internally used base for iterator and const_iterator
        template<bool Const> class Iterator;

        //- An iterator with const access to HashTable internals.
        friend class Iterator<true>;

        //- An iterator with non-const access to HashTable internals.
        friend class Iterator<false>;


    // Protected Member Functions

        //- Remove using begin/end iterators of listed keys
        template<class InputIter>
        inline label eraseMultiple
        (
            const InputIter begIter,
            const InputIter endIter
        );


public:

    // Constructors

        //- Construct null with default (128) table capacity
        HashTable();

        //- Construct given initial table capacity
        explicit HashTable(const label size);

        //- Construct from Istream with default table capacity
        HashTable(Istream& is, const label size = 128);

        //- Construct as copy
        HashTable(const this_type& ht);

        //- Move construct
        HashTable(this_type&& ht);

        //- Construct by transferring the parameter contents
        HashTable(const Xfer<this_type>& ht);

        //- Construct from an initializer list
        HashTable(std::initializer_list<std::pair<Key, T>> lst);


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

        //- Return true if hashed entry is found in table
        inline bool found(const Key& key) const;

        //- Find and return an iterator set at the hashed entry
        //  If not found iterator = end()
        inline iterator find(const Key& key);

        //- Find and return an const_iterator set at the hashed entry
        //  If not found iterator = end()
        inline const_iterator find(const Key& key) const;

        //- Find and return an const_iterator set at the hashed entry
        //  If not found iterator = end()
        inline const_iterator cfind(const Key& key) const;

        //- Return hashed entry if it exists, or return the given default
        inline const T& lookup(const Key& key, const T& deflt) const;


      // Table of contents

        //- The table of contents (the keys) in unsorted order.
        List<Key> toc() const;

        //- The table of contents (the keys) in sorted order
        List<Key> sortedToc() const;

        //- The table of contents (the keys) sorted according to the
        //- specified comparator
        template<class Compare>
        List<Key> sortedToc(const Compare& comp) const;

        //- The table of contents (the keys) selected according to the
        //- unary predicate applied to the \b keys.
        //  \param invert changes the logic to select when the predicate
        //      is false
        //  \return sorted list of selected keys
        template<class UnaryPredicate>
        List<Key> tocKeys
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- The table of contents (the keys) selected according to the
        //- unary predicate applied to the \b values.
        //  \param invert changes the logic to select when the predicate
        //      is false
        //  \return sorted list of selected keys
        template<class UnaryPredicate>
        List<Key> tocValues
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- The table of contents (the keys) selected according to the
        //- binary predicate applied to the \b keys and \b values.
        //  \param invert changes the logic to select when the predicate
        //      is false
        //  \return sorted list of selected keys
        template<class BinaryPredicate>
        List<Key> tocEntries
        (
            const BinaryPredicate& pred,
            const bool invert = false
        ) const;


      // Counting

        //- Count the number of keys that satisfy the unary predicate
        //  \param invert changes the logic to select when the predicate
        //      is false
        template<class UnaryPredicate>
        label countKeys
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- Count the number of values that satisfy the unary predicate
        //  \param invert changes the logic to select when the predicate
        //      is false
        template<class UnaryPredicate>
        label countValues
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- Count the number of entries that satisfy the binary predicate.
        //  \param invert changes the logic to select when the predicate
        //      is false
        template<class BinaryPredicate>
        label countEntries
        (
            const BinaryPredicate& pred,
            const bool invert = false
        ) const;


      // Edit

        //- Insert a new entry, not overwriting existing entries.
        //  \return True if the entry inserted, which means that it did
        //  not previously exist in the table.
        inline bool insert(const Key& key, const T& obj);

        //- Assign a new entry, overwriting existing entries.
        //
        //  \return True, since it always overwrites any entries.
        inline bool set(const Key& key, const T& obj);

        //- Erase an entry specified by given iterator
        //  This invalidates the iterator until the next ++ operation.
        //
        //  Includes a safeguard against the end-iterator such that the
        //  following is safe:
        //  \code
        //      auto iter = table.find(unknownKey);
        //      table.erase(iter);
        //  \endcode
        //  which is what \code table.erase(unknownKey) \endcode does anyhow.
        //
        //  \return True if the corresponding entry existed and was removed
        bool erase(const iterator& iter);

        //- Erase an entry specified by the given key
        //  \return True if the entry existed and was removed
        bool erase(const Key& key);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        label erase(const UList<Key>& keys);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        template<unsigned Size>
        label erase(const FixedList<Key, Size>& keys);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        label erase(std::initializer_list<Key> keys);

        //- Remove table entries given by keys of the other hash-table.
        //
        //  The other hash-table must have the same type of key, but the
        //  type of values held and the hashing function are arbitrary.
        //
        //  \return The number of items removed
        template<class AnyType, class AnyHash>
        label erase(const HashTable<AnyType, Key, AnyHash>& other);

        //- Retain table entries given by keys of the other hash-table.
        //
        //  The other hash-table must have the same type of key, but the
        //  type of values held and the hashing function are arbitrary.
        //
        //  \return The number of items changed (removed)
        template<class AnyType, class AnyHash>
        label retain(const HashTable<AnyType, Key, AnyHash>& other);

        //- Generalized means to filter table entries based on their keys.
        //  Keep (or optionally prune) entries with keys that satisfy
        //  the unary predicate, which has the following signature:
        //  \code
        //  bool operator()(const Key& k);
        //  \endcode
        //
        //  For example,
        //  \code
        //  wordRes goodFields = ...;
        //  allFieldNames.filterKeys
        //  (
        //      [&goodFields](const word& k){ return goodFields.match(k); }
        //  );
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterKeys
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Generalized means to filter table entries based on their values.
        //  Keep (or optionally prune) entries with values that satisfy
        //  the unary predicate, which has the following signature:
        //  \code
        //  bool operator()(const T& v);
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterValues
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Generalized means to filter table entries based on their key/value.
        //  Keep (or optionally prune) entries with keys/values that satisfy
        //  the binary predicate, which has the following signature:
        //  \code
        //  bool operator()(const Key& k, const T& v);
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class BinaryPredicate>
        label filterEntries
        (
            const BinaryPredicate& pred,
            const bool pruning = false
        );


        //- Resize the hash table for efficiency
        void resize(const label sz);

        //- Clear all entries from table
        void clear();

        //- Clear the table entries and the table itself.
        //  Equivalent to clear() followed by resize(0)
        void clearStorage();

        //- Swap contents of the argument table into this table
        void swap(HashTable<T, Key, Hash>& ht);

        //- Transfer the contents of the argument table into this table
        //  and annul the argument table.
        void transfer(HashTable<T, Key, Hash>& ht);

        //- Transfer contents to the Xfer container
        inline Xfer<HashTable<T, Key, Hash>> xfer();


    // Member Operators

        //- Find and return a hashed entry. FatalError if it does not exist.
        inline T& operator[](const Key& key);

        //- Find and return a hashed entry. FatalError if it does not exist.
        inline const T& operator[](const Key& key) const;

        //- Return existing entry or create a new entry.
        //  A newly created entry is created as a nameless T() and is thus
        //  value-initialized. For primitives, this will be zero.
        inline T& operator()(const Key& key);

        //- Return existing entry or insert a new entry.
        inline T& operator()(const Key& key, const T& deflt);

        //- Return hashed entry if it exists, or return the given default
        inline const T& operator()(const Key& key, const T& deflt) const;

        //- Assignment
        void operator=(const HashTable<T, Key, Hash>& rhs);

        //- Assignment from an initializer list
        void operator=(std::initializer_list<std::pair<Key, T>> lst);

        //- Move assign
        void operator=(HashTable<T, Key, Hash>&& rhs);

        //- Equality. Tables are equal if all keys and values are equal,
        //- independent of order or underlying storage size.
        bool operator==(const HashTable<T, Key, Hash>& rhs) const;

        //- The opposite of the equality operation.
        bool operator!=(const HashTable<T, Key, Hash>& rhs) const;


protected:

    // Iterators and helpers

        //- The iterator base for HashTable (internal use only).
        //  Note: data and functions are protected, to allow reuse by iterator
        //  and prevent most external usage.
        //  iterator and const_iterator have the same size, allowing
        //  us to reinterpret_cast between them (if desired)

        template<bool Const>
        class Iterator
        {
        public:
          // Typedefs
            using iterator_category = std::forward_iterator_tag;
            using difference_type = this_type::difference_type;

            //- The HashTable container type
            using table_type = typename std::conditional
            <
                Const,
                const this_type,
                this_type
            >::type;

            //- The node-type being addressed
            using node_type = typename std::conditional
            <
                Const,
                const this_type::node_type,
                this_type::node_type
            >::type;

            //- The key type
            using key_type = this_type::key_type;

            //- The object type being addressed
            using mapped_type = typename std::conditional
            <
                Const,
                const this_type::mapped_type,
                this_type::mapped_type
            >::type;


          // Member Functions

            //- True if iterator points to an entry
            //  This can be used directly instead of comparing to end()
            inline bool found() const;

            //- The key associated with the iterator
            inline const Key& key() const;

          // Member Operators

            //- Compare hash-entry element pointers.
            //  Independent of const/non-const access
            inline bool operator==(const Iterator<true>& iter) const;
            inline bool operator!=(const Iterator<true>& iter) const;

            inline bool operator==(const Iterator<false>& iter) const;
            inline bool operator!=(const Iterator<false>& iter) const;


        protected:
            friend class HashTable;  // For begin/find constructors

          // Protected Data

            //- The selected entry.
            //  MUST be the first member for easy comparison between iterators
            //  and for reinterpret_cast from nullObject
            node_type* entry_;

            //- The hash-table container being iterated on.
            //  Using a pointer allows default bitwise copy/assignment
            table_type* container_;

            //- Index within the hash-table data.
            //  A signed value, since iterator_erase() needs a negative value
            //  to mark the position.
            label index_;


          // Protected Constructors

            //- Construct null (end iterator)
            inline Iterator();

            //- Construct from begin of hash-table
            inline Iterator(bool, table_type* tbl);

            //- Construct by finding key in hash table
            Iterator(table_type* tbl, const Key& key);


          // Protected Member Functions

            //- Increment to the next position
            inline void increment();

            //- The object associated with the iterator
            inline mapped_type& object() const
            {
                return entry_->mapped();
            }

            //- Permit an explicit cast to the other (const/non-const) searcher
            inline explicit operator const Iterator<!Const>&() const
            {
                return *reinterpret_cast<const Iterator<!Const>*>(this);
            }
        };


        //- Low-level entry erasure using iterator internals.
        //  This invalidates the iterator until the next ++ operation.
        //  \return True if the corresponding entry existed and was removed
        bool iterator_erase(node_type*& entry, label& index);

public:

        //- Forward iterator with non-const access
        class iterator
        :
            public Iterator<false>
        {
        public:
          // Typedefs
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = this_type::difference_type;

            using key_type    = this_type::key_type;
            using mapped_type = this_type::mapped_type;
            using value_type  = this_type::value_type;
            using pointer     = this_type::pointer;
            using reference   = this_type::reference;

          // Constructors

            //- Construct null (end iterator)
            inline iterator() {}

            //- Copy construct from similar access type
            inline explicit iterator(const Iterator<false>& iter)
            :
                Iterator<false>(iter)
            {}


          // Member functions/operators

            //- Non-const access to referenced object
            using Iterator<false>::object;

            //- Non-const access to referenced object
            inline reference operator*() const  { return this->object(); }
            inline reference operator()() const { return this->object(); }

            inline iterator& operator++();
            inline iterator operator++(int);
        };


    // STL const_iterator

        //- Forward iterator with const access
        class const_iterator
        :
            public Iterator<true>
        {
        public:
          // Typedefs
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = this_type::difference_type;

            using key_type    = this_type::key_type;
            using mapped_type = const this_type::mapped_type;
            using value_type  = const this_type::value_type;
            using pointer     = this_type::const_pointer;
            using reference   = this_type::const_reference;

          // Constructors

            //- Construct null (end iterator)
            inline const_iterator() {}

            //- Copy construct from similar access type
            inline explicit const_iterator(const Iterator<true>& iter)
            :
                Iterator<true>(iter)
            {}

            //- Copy construct from dissimilar access type
            inline explicit const_iterator(const Iterator<false>& iter)
            :
                Iterator<true>
                (
                    static_cast<const Iterator<true>&>(iter)
                )
            {}

            //- Implicit conversion from dissimilar access type
            inline const_iterator(const iterator& iter)
            :
                const_iterator(reinterpret_cast<const const_iterator&>(iter))
            {}


          // Member functions/operators

            //- Const access to referenced object
            using Iterator<true>::object;

            //- Const access to referenced object
            inline reference operator*() const  { return this->object(); }
            inline reference operator()() const { return this->object(); }

            inline const_iterator& operator++();
            inline const_iterator operator++(int);

          // Assignment

            const_iterator& operator=(const const_iterator&) = default;

            // Allow assign from iterator to const_iterator
            const_iterator& operator=(const iterator& iter)
            {
                return this->operator=
                (
                    reinterpret_cast<const const_iterator&>(iter)
                );
            }
        };


    //- Iterating over keys only

        //- An iterator wrapper for returning a reference to the key
        template<class Iter>
        class key_iterator_base
        :
            public Iter
        {
        public:
            using value_type = this_type::key_type;
            using pointer    = const Key*;
            using reference  = const Key&;

            //- Implicit conversion
            inline key_iterator_base(const Iter& iter)
            :
                Iter(iter)
            {}

            //- Return the key
            inline reference operator*() const  { return this->key(); }
            inline reference operator()() const { return this->key(); }

            inline key_iterator_base& operator++()
            {
                this->increment();
                return *this;
            }

            inline key_iterator_base operator++(int)
            {
                key_iterator_base iter(*this);
                this->increment();
                return iter;
            }
        };


        //- Forward iterator returning the key
        using key_iterator = key_iterator_base<iterator>;

        //- Forward const iterator returning the key
        using const_key_iterator = key_iterator_base<const_iterator>;

        //- A const iterator begin/end pair for iterating over keys
        const_iterator_pair<const_key_iterator, this_type> keys() const
        {
            return
                const_iterator_pair<const_key_iterator,this_type>(*this);
        }


    // Iterator access

        //- iterator set to the beginning of the HashTable
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


    // Writing

        //- Print information
        Ostream& printInfo(Ostream& os) const;

        //- Write the unordered keys as a list, with line-breaks if list length
        //- exceeds shortListLen.
        //  Using '0' suppresses line-breaks entirely.
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


// Global Functions

// Exchange contents of hash tables - see HashTable::swap().
template<class T, class Key, class Hash>
inline void Swap
(
    HashTable<T, Key, Hash>& a,
    HashTable<T, Key, Hash>& b
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HashTableI.H"
#include "HashTableIterI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef NoHashTableC
#ifdef NoRepository
    #include "HashTable.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
