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
    The iterator '()' dereferencing operator returns a proxy class that handles
    assignment. However, the iterator itself can also handle assignment
    directly and also returns a bool as per the set() methods.
    Thus the following bit of code works as expected:

    @code
        blist[5] = 4;
        changed = (blist[5] = 10);
        changed = blist.set(5, 8);
    @endcode

    Using set() might be more readable (as expected) though.

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

template <int nBits>
class PackedList
:
    private List<unsigned int>
{
    typedef unsigned int PackedStorage;

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

    // Forward declaration of iterator and const_iterator
        class iterator;
        class const_iterator;

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
        inline const List<unsigned int>& storage() const;

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
        //  Does not auto-vivify elements.
        inline unsigned int operator[](const label) const;

        //- Set value at index I.
        //  Returns iterator to perform the actual operation.
        //  Auto-vivifies any new values to zero.
        inline iterator operator[](const label);

        //- Assignment of all entries to the given value. Takes linear time.
        inline void operator=(const unsigned int val);

        //- Assignment operator. Takes linear time.
        void operator=(const PackedList<nBits>&);

        //- Assignment operator. Takes linear time.
        void operator=(const UList<label>&);

    // Ostream operator

    //  // Write PackedList to Ostream.
    //  friend Ostream& operator<< <nBits> (Ostream&, const PackedList<nBits>&);

    //- The const_iterator for PackedList
    //  Note: make data and functions protected, to allow reuse by iterator.
    //  Try not to be disturbed by non-const methods such as set(), they are
    //  just inherited by iterator.
    class const_iterator
    {
        friend class PackedList;

    protected:

        // Protected Data

        //- Reference to original list
        PackedList& list_;

        //- Element index within storage
        unsigned index_;

        //- Offset within storage element
        unsigned offset_;

        // Protected Member Functions

        //- Get value as unsigned
        inline unsigned int get() const;

        //- Set value, returning true if changed
        //  This is obviously used by iterator and not by const_iterator
        inline bool set(unsigned int);

    public:

        // Constructors

        //- Construct from base list and position index
        inline const_iterator(const PackedList&, const label);

        //- Construct from non-const iterator
        explicit inline const_iterator(const iterator&);

        // Member operators

        inline bool operator==(const const_iterator&) const;
        inline bool operator!=(const const_iterator&) const;

        //- Return referenced value directly
        inline unsigned int operator*() const;

        //- Return referenced value directly
        inline unsigned int operator()() const;

        inline const_iterator& operator++();
        inline const_iterator operator++(int);

        //- Conversion operator
        inline operator unsigned int () const;

        //- Print value and information
        Ostream& print(Ostream&) const;
    };


    //- const_iterator set to the beginning of the PackedList
    inline const_iterator cbegin() const;

    //- const_iterator set to beyond the end of the PackedList
    inline const_iterator cend() const;

    //- const_iterator set to the beginning of the PackedList
    inline const_iterator begin() const;

    //- const_iterator set to beyond the end of the PackedList
    inline const_iterator end() const;


    //- The iterator class used for PackedList
    class iterator
    :
        public const_iterator
    {
        friend class PackedList;

        //- A private proxy reference to the iterator, used for setting values
        //  Also provides a uniform target for the iterator 'operator*'
        //  and 'operator()' methods
        class iteratorProxy
        {
            iterator* iter_;
        public:

            //- Construct from pointer to iterator
            inline iteratorProxy(iterator*);

            //- Assign value, returns true if the value changed
            inline bool operator=(const unsigned int val);

            //- Conversion operator
            inline operator unsigned int () const;
        };

    public:

        // Constructors

        //- Construct from base list and position index
        inline iterator(const PackedList&, const label);

        // Member Operators
        inline void operator=(const iterator&);

        //- Assign value directly, returns true if the value changed
        inline bool operator=(unsigned int);

        //- Return value
        inline unsigned int operator*() const;

        //- Return value
        inline unsigned int operator()() const;

        //- Return proxy for assigning the referenced value
        inline iteratorProxy operator*();

        //- Return proxy for assigning the referenced value
        inline iteratorProxy operator()();

        inline iterator& operator++();
        inline iterator operator++(int);
    };

    //- iterator set to the beginning of the PackedList
    inline iterator begin();

    //- iterator set to beyond the end of the HashTable
    inline iterator end();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "PackedListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PackedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
