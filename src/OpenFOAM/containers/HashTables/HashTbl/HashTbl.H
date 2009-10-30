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

Class
    Foam::HashTbl

Description
    An STL-conforming hash table.

Note
    Hashing index collisions are handled via chaining using a singly-linked
    list with the colliding entry being added to the head of the linked
    list. Thus copying the hash table (or indeed even resizing it) will
    often result in a different hash order. Use a sorted table-of-contents
    when the hash order is important.

SourceFiles
    HashTblI.H
    HashTbl.C
    HashTblIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashTbl_H
#define HashTbl_H

#include "label.H"
#include "uLabel.H"
#include "word.H"
#include "Xfer.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class List;
template<class T> class UList;
template<class T, class Key, class Hash> class HashTbl;
template<class T, class Key, class Hash> class HashPtrTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream&, HashTbl<T, Key, Hash>&);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream&, const HashTbl<T, Key, Hash>&);


/*---------------------------------------------------------------------------*\
                        Class HashTblName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(HashTbl);


/*---------------------------------------------------------------------------*\
                          Class HashTbl Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class HashTbl
:
    public HashTblName
{
    // Private data type for table entries

        struct hashedEntry
        {
            //- The lookup key
            Key key_;

            //- Pointer to next hashedEntry in sub-list
            hashedEntry* next_;

            //- The data object
            T obj_;

            //- Constructors

                //- Construct given key, next pointer and object
                inline hashedEntry
                (
                    const Key&,
                    hashedEntry* next,
                    const T& newEntry
                );

                //- Dissallow construction as copy
                hashedEntry(const hashedEntry&);
        };


    // Private data: size of table, the table and current number of elements

        //- The current number of elements in table
        label nElmts_;

        //- Number of primary entries allocated in table (not necessarily used)
        label tableSize_;

        //- The table of primary entries
        hashedEntry** table_;


    // Private Member Functions

        //- Return a canonical (power-of-two) size
        static label canonicalSize(const label);

        //- Return the hash index of the Key within the current table size.
        //  No checks for zero-sized tables.
        inline label hashKeyIndex(const Key&) const;

        //- Assign a new hashedEntry to a possibly already existing key
        bool set(const Key&, const T& newElmt, bool protect);

public:

        //- Declare friendship with the HashPtrTable class
        template<class T2, class Key2, class Hash2>
        friend class HashPtrTable;


    // Forward declaration of STL iterators

        class iterator;
        friend class iterator;

        class const_iterator;
        friend class const_iterator;


    // Constructors

        //- Construct given initial table size
        HashTbl(const label size = 128);

        //- Construct from Istream
        HashTbl(Istream&, const label size = 128);

        //- Construct as copy
        HashTbl(const HashTbl<T, Key, Hash>&);

        //- Construct by transferring the parameter contents
        HashTbl(const Xfer<HashTbl<T, Key, Hash> >&);


    // Destructor

        ~HashTbl();


    // Member Functions

        // Access

            //- The size of the underlying table
            inline label capacity() const;

            //- Return number of elements in table.
            inline label size() const;

            //- Return true if the hash table is empty
            inline bool empty() const;

            //- Return true if hashedEntry is found in table
            bool found(const Key&) const;

            //- Find and return an iterator set at the hashedEntry
            //  If not found iterator = end()
            iterator find(const Key&);

            //- Find and return an const_iterator set at the hashedEntry
            //  If not found iterator = end()
            const_iterator find(const Key&) const;

            //- Return the table of contents
            List<Key> toc() const;

            //- Return the table of contents as a sorted list
            List<Key> sortedToc() const;

            //- Print information
            Ostream& printInfo(Ostream&) const;

        // Edit

            //- Insert a new hashedEntry
            inline bool insert(const Key&, const T& newElmt);

            //- Assign a new hashedEntry, overwriting existing entries
            inline bool set(const Key&, const T& newElmt);

            //- Erase an hashedEntry specified by given iterator
            bool erase(const iterator&);

            //- Erase an hashedEntry specified by given key if in table
            bool erase(const Key&);

            //- Remove entries given by the listed keys from this HashTbl
            //  Return the number of elements removed
            label erase(const UList<Key>&);

            //- Remove entries given by the given keys from this HashTbl
            //  Return the number of elements removed.
            //  The parameter HashTbl needs the same type of keys, but
            //  but the type of values held is arbitrary.
            template<class AnyType>
            label erase(const HashTbl<AnyType, Key, Hash>&);

            //- Resize the hash table for efficiency
            void resize(const label newSize);

            //- Clear all entries from table
            void clear();

            //- Clear the table entries and the table itself.
            //  Equivalent to clear() followed by resize(0)
            void clearStorage();

            //- Shrink the allocated table to approx. twice number of elements
            void shrink();

            //- Transfer the contents of the argument table into this table
            //  and annull the argument table.
            void transfer(HashTbl<T, Key, Hash>&);

            //- Transfer contents to the Xfer container
            inline Xfer<HashTbl<T, Key, Hash> > xfer();


    // Member Operators

        //- Find and return an hashedEntry
        inline T& operator[](const Key&);

        //- Find and return an hashedEntry
        inline const T& operator[](const Key&) const;

        //- Find and return an hashedEntry, create it null if not present.
        inline T& operator()(const Key&);

        //- Assignment
        void operator=(const HashTbl<T, Key, Hash>&);

        //- Equality. Two hash tables are equal if all contents of first are
        //  also in second and vice versa. So does not depend on table size or
        //  order!
        bool operator==(const HashTbl<T, Key, Hash>&) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const HashTbl<T, Key, Hash>&) const;



    // STL type definitions

        //- Type of values the HashTbl contains.
        typedef T value_type;

        //- Type that can be used for storing into HashTbl::value_type
        //  objects.  This type is usually List::value_type&.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  HashTbl::value_type objects.  This type is usually const
        //  HashTbl::value_type&.
        typedef const T& const_reference;

        //- The type that can represent the size of a HashTbl.
        typedef label size_type;


    // STL iterator

        //- An STL-conforming iterator
        class iterator
        {
            friend class HashTbl;
            friend class const_iterator;

            // Private data

                //- Reference to the HashTbl this is an iterator for
                HashTbl<T, Key, Hash>& hashTable_;

                //- Current element
                hashedEntry* elmtPtr_;

                //- Current hash index
                label hashIndex_;

        public:

            // Constructors

                //- Construct from hash table, element and hash index
                inline iterator
                (
                    HashTbl<T, Key, Hash>& curHashTbl,
                    hashedEntry* elmt,
                    label hashIndex
                );

            // Member operators

                inline void operator=(const iterator&);

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                inline bool operator==(const const_iterator&) const;
                inline bool operator!=(const const_iterator&) const;

                inline T& operator*();
                inline T& operator()();

                inline const T& operator*() const;
                inline const T& operator()() const;

                inline iterator& operator++();
                inline iterator operator++(int);

                inline const Key& key() const;
        };


        //- iterator set to the begining of the HashTbl
        inline iterator begin();

        //- iterator set to beyond the end of the HashTbl
        inline const iterator& end();


    // STL const_iterator

        //- An STL-conforming const_iterator
        class const_iterator
        {
            friend class iterator;

            // Private data

                //- Reference to the HashTbl this is an iterator for
                const HashTbl<T, Key, Hash>& hashTable_;

                //- Current element
                const hashedEntry* elmtPtr_;

                //- Current hash index
                label hashIndex_;


        public:

            // Constructors

                //- Construct from hash table, element and hash index
                inline const_iterator
                (
                    const HashTbl<T, Key, Hash>& curHashTbl,
                    const hashedEntry* elmt,
                    label hashIndex
                );

                //- Construct from the non-const iterator
                inline const_iterator(const iterator&);


            // Member operators

                inline void operator=(const const_iterator&);

                inline bool operator==(const const_iterator&) const;
                inline bool operator!=(const const_iterator&) const;

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                inline const T& operator*() const;
                inline const T& operator()() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);

                inline const Key& key() const;
        };


        //- const_iterator set to the beginning of the HashTbl
        inline const_iterator cbegin() const;

        //- const_iterator set to beyond the end of the HashTbl
        inline const const_iterator& cend() const;

        //- const_iterator set to the beginning of the HashTbl
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the HashTbl
        inline const const_iterator& end() const;


    // IOstream Operator

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream&,
            HashTbl<T, Key, Hash>&
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream&,
            const HashTbl<T, Key, Hash>&
        );


private:

        //- iterator returned by end()
        iterator endIter_;

        //- const_iterator returned by end()
        const_iterator endConstIter_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "HashTblI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef NoHashTblC
#ifdef NoRepository
#   include "HashTbl.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
