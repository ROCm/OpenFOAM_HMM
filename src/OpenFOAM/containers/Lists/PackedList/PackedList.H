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
    Foam::PackedList

Description
    A Dynamically allocatable list of packed unsigned ints.

    The list resizing is similar to DynamicList, thus the methods clear()
    and setSize() behave like their DynamicList counterparts and the methods
    reserve() and setCapacity() can be used to influence the allocation.

    The number of bits per item is specified by the template parameter nBits.

Note
    In a const context, the '[]' operator simply returns the stored value,
    with out-of-range elements returned as zero.
    In a non-const context, the '[]' operator returns an iteratorBase, which
    may not have a valid reference for out-of-range elements.
    The iteratorBase class handles the assignment of new values.

    Using the iteratorBase as a proxy allows assignment of values
    between list elements. Thus the following bit of code works as expected:
    @code
        blist[1] = blist[5];     // value assignment, not iterator position
        blist[2] = blist[5] = 4; // propagates value
        blist[1] = blist[5] = blist[6];  // propagates value
    @endcode

    Using get() or the '[]' operator are similarly fast. Looping and reading
    with an iterator is approx. 15% slower, but can be more flexible.

    Using the set() operator (and the '[]' operator) are marginally slower
    (approx. 5%) than using an iterator, but the set() method has an
    advantage that it also returns a bool if the value changed.  This can be
    useful for branching on changed values.

    @code
        blist[5] = 4;
        changed = blist.set(5, 8);
        if (changed) ...
    @endcode

SeeAlso
    Foam::DynamicList

ToDo
    Checks for bad template parameters (ie, nBits=0, nBits too large)?

SourceFiles
    PackedListI.H
    PackedList.C

\*---------------------------------------------------------------------------*/

#ifndef PackedList_H
#define PackedList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<int nBits> class PackedList;

// template<int nBits>
// Ostream& operator<<(Ostream&, const PackedList<nBits>&);


/*---------------------------------------------------------------------------*\
                      Class PackedListName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(PackedList);

/*---------------------------------------------------------------------------*\
                         Class PackedList Declaration
\*---------------------------------------------------------------------------*/

template <int nBits=1>
class PackedList
:
    private List<unsigned int>
{
    typedef unsigned int      StorageType;
    typedef List<StorageType> StorageList;

    // Private data

        //- Number of nBits entries
        label size_;

    // Private Member Functions

        //- Calculate the list length when packed
        inline static label packedLength(const label);

        //- Check index I is within valid range [ 0 .. max_value() ]
        inline void checkIndex(const label) const;

public:

    // Public data

        //- The max. number of bits that can be templated.
        //  Might someday be useful for a template assert.
        inline static unsigned int max_bits();

        //- The max. value for an entry, which simultaneously the bit-mask
        //  eg, ((1 << 2) - 1) yields 0b0011
        inline static unsigned int max_value();

        //- The number of entries per packed storage element
        inline static unsigned int packing();

        //- Masking for all bits below the offset
        inline static unsigned int maskLower(unsigned offset);

    // Forward declaration of iteratorBase

        class iteratorBase;

    // Constructors

        //- Null constructor
        inline PackedList();

        //- Construct with given size, initializes list to 0.
        inline PackedList(const label size);

        //- Construct with given size and value for all elements.
        PackedList(const label size, const unsigned val);

        //- Copy constructor.
        inline PackedList(const PackedList<nBits>&);

        //- Construct by transferring the parameter contents
        inline PackedList(const Xfer<PackedList<nBits> >&);

        //- Construct from a list of labels
        PackedList(const UList<label>&);

        //- Clone
        inline autoPtr< PackedList<nBits> > clone() const;

    // Member Functions

        // Access

        //- The number of elements that can be stored before resizing
        inline label capacity() const;

        //- Number of entries.
        inline label size() const;

        //- Return true if the list is empty (i.e., if size() == 0).
        inline bool empty() const;

        //- Get value at index I.
        //  Does not auto-vivify entries.
        inline unsigned int get(const label) const;

        //- Set value at index I. Return true if value changed.
        //  Does not auto-vivify entries.
        inline bool set(const label, const unsigned int val);

        //- Return the underlying packed storage
        inline List<unsigned int>& storage();

        //- Return the underlying packed storage
        inline const List<unsigned int>& storage() const;

        //- Count number of bits set, O(log(n))
        //  Uses the Hamming weight (population count) method
        //  http://en.wikipedia.org/wiki/Hamming_weight
        unsigned int count() const;

        //- Trim any trailing zero elements
        bool trim();

        //- Return the values as a labelList
        labelList values() const;

        //- Print values and information
        Ostream& print(Ostream&) const;

        // Edit

        //- Alter the size of the underlying storage.
        //  The addressed size will be truncated if needed to fit, but will
        //  remain otherwise untouched.
        inline void setCapacity(const label);

        //- Reset addressable list size, does not shrink the allocated size.
        //  Optionally specify a value for new elements.
        inline void resize(const label, const unsigned int& val = 0);

        //- Alias for resize()
        inline void setSize(const label, const unsigned int& val = 0);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size.
        //  Optionally provide an initialization value for new elements.
        inline void reserve(const label);

        //- Clear the list, i.e. set addressable size to zero.
        //  Does not adjust the underlying storage
        inline void clear();

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Shrink the allocated space to what is used.
        inline void shrink();

        //- Transfer the contents of the argument list into this list
        //  and annull the argument list.
        inline void transfer(PackedList<nBits>&);

        //- Transfer contents to the Xfer container
        inline Xfer<PackedList<nBits> > xfer();


    // Member operators

        //- Append a value at the end of the list
        inline void append(const unsigned int val);

        //- Remove and return the last element
        inline unsigned int remove();

        //- Get value at index I
        //  Does not auto-vivify entries.
        inline unsigned int operator[](const label) const;

        //- Set value at index I.
        //  Returns iterator to perform the actual operation.
        //  Does not auto-vivify entries.
        inline iteratorBase operator[](const label);

        //- Assignment of all entries to the given value. Takes linear time.
        inline void operator=(const unsigned int val);

        //- Assignment operator. Takes linear time.
        void operator=(const PackedList<nBits>&);

        //- Assignment operator. Takes linear time.
        void operator=(const UList<label>&);

    // Ostream operator

//      // Write PackedList to Ostream.
//      friend Ostream& operator<< <nBits> (Ostream&, const PackedList<nBits>&);

    // Iterators and helpers

        //- The iterator base for PackedList
        //  Note: data and functions are protected, to allow reuse by iterator
        //  and prevent most external usage.
        class iteratorBase
        {
            friend class PackedList;

        protected:

            // Protected Data

                //- Pointer to original list
                //  This also lets us use the default bitwise copy/assignment
                PackedList* list_;

                //- Element index within storage
                unsigned index_;

                //- Offset within storage element
                unsigned offset_;

            // Protected Member Functions

                //- Get value as unsigned
                inline unsigned int get() const;

                //- Set value, returning true if changed
                inline bool set(unsigned int);

                //- Increment to new position
                inline void incr();

                //- Decrement to new position
                inline void decr();

                //- Move to new position, but not beyond end()
                inline void seek(const iteratorBase&);


            // Constructors

                //- Construct null
                inline iteratorBase();

                //- Copy construct
                inline iteratorBase(const iteratorBase&);

                //- Construct from base list and position index
                inline iteratorBase(const PackedList*, const label);

        public:

            // Member Functions

                //- Return true if the element is within addressable range
                inline bool valid() const;

                //- Move iterator to end() if it would otherwise be out-of-range
                //  Returns true if the element was already ok
                inline bool validate();

            // Member Operators

                //- Compare positions
                inline bool operator==(const iteratorBase&) const;
                inline bool operator!=(const iteratorBase&) const;

                //- Assign value, not position.
                //  This allows packed[0] = packed[3] for assigning values
                inline unsigned int operator=(const iteratorBase&);

                //- Assign value
                inline unsigned int operator=(const unsigned int val);

                //- Conversion operator
                inline operator unsigned int () const;

            //- Print value and information
            Ostream& print(Ostream&) const;
        };


        //- The iterator class used for PackedList
        class iterator
        :
            public iteratorBase
        {

            //- Should never violate const-ness!
            void operator=(const const_iterator&);

        public:

            // Constructors

                //- Construct null
                inline iterator();

                //- Construct from iterator base, eg iter(packedlist[i])
                inline iterator(const iteratorBase&);

                //- Construct from base list and position index
                inline iterator(const PackedList*, const label);

            // Member Operators

                //- Assign from iteratorBase, eg iter = packedlist[i]
                inline iterator& operator=(const iteratorBase&);

                //- Return value
                inline unsigned int operator*() const;

                //- Return value
                inline unsigned int operator()() const;

                //- Return iteratorBase for assigning values
                inline iteratorBase& operator*();

                //- Return iteratorBase for assigning values
                inline iteratorBase& operator()();

                inline iterator& operator++();
                inline iterator operator++(int);

                inline iterator& operator--();
                inline iterator operator--(int);

        };

        //- iterator set to the beginning of the PackedList
        inline iterator begin();

        //- iterator set to beyond the end of the HashTable
        inline iterator end();


        //- The const_iterator for PackedList
        class const_iterator
        :
            public iteratorBase
        {
        public:

            // Constructors

                //- Construct null
                inline const_iterator();

                //- Construct from iterator base
                inline const_iterator(const iteratorBase&);

                //- Construct from base list and position index
                inline const_iterator(const PackedList*, const label);

                //- Construct from non-const iterator
                inline const_iterator(const iterator&);

            // Member operators

                //- Assign from iteratorBase or derived
                //  eg, iter = packedlist[i] or iter = [non-const]list.begin()
                inline const_iterator& operator=(const iteratorBase&);

                //- Return referenced value directly
                inline unsigned int operator*() const;

                //- Return referenced value directly
                inline unsigned int operator()() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);

                inline const_iterator& operator--();
                inline const_iterator operator--(int);

        };


        //- const_iterator set to the beginning of the PackedList
        inline const_iterator cbegin() const;

        //- const_iterator set to beyond the end of the PackedList
        inline const_iterator cend() const;

        //- const_iterator set to the beginning of the PackedList
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the PackedList
        inline const_iterator end() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PackedListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PackedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
