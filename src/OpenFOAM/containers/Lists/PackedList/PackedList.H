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

    Gets given the number of bits per item.

Note
    The list resizing is similar to DynamicList, thus the methods clear()
    and setSize() behave like their DynamicList counterparts and the methods
    reserve() and setCapacity() can be used to influence the allocation.

SeeAlso
    Foam::DynamicList

ToDo
    Add checks for bad template parameters (ie, nBits=0, nBits too large).
    The missing const_iterator might eventually still be needed.

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

    // Forward declaration of iterator
        class iterator;
        friend class iterator;


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
        //  Does not auto-vivifies entries.
        inline unsigned int get(const label) const;

        //- Set value at index I. Return true if value changed.
        //  Does not auto-vivifies entries.
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
        inline void reserve(const label, const unsigned int& val = 0);

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

        //- Get value at index I
        //  Auto-vivifies any new values to zero.
        inline unsigned int operator[](const label) const;

        //- Set value at index I.
        //  Returns proxy to perform the actual operation.
        //  Auto-vivifies any new values to zero.
        inline iterator operator[](const label);

        //- Assignment of all entries to the given value.
        //  Does set on all elements.
        inline void operator=(const unsigned int val);

        //- Assignment operator. Takes linear time.
        void operator=(const PackedList<nBits>&);

        //- Assignment operator. Takes linear time.
        void operator=(const UList<label>&);

    // Ostream operator

    //  // Write PackedList to Ostream.
    //  friend Ostream& operator<< <nBits> (Ostream&, const PackedList<nBits>&);


    //- The iterator-like class used for PackedList
    class iterator
    {
        friend class PackedList;

        // Private Data

        //- Reference to original list
        PackedList& list_;

        //- Element index within storage
        unsigned elem_;

        //- Offset within storage element
        unsigned offset_;

        //- Return the raw storage element
        inline PackedStorage& chunk() const;

        //- Return the position in the PackedList
        inline label position() const;

    public:

        // Constructors

        //- Construct from base list and position index
        inline iterator(const PackedList&, const label i);

        // Members

        //- Assign value, returns true if the value changed
        inline bool operator=(const unsigned int val);

        //- Conversion operator
        inline operator unsigned int () const;

        //- Conversion operator
        inline operator bool() const;

        // Member operators

        inline void operator=(const iterator&);
        inline bool operator==(const iterator&) const;
        inline bool operator!=(const iterator&) const;

        //- Return referenced value
        inline unsigned int operator()() const;

        inline iterator& operator++();
        inline iterator operator++(int);

        //- Print value and information
        Ostream& print(Ostream&) const;

    };


    //- iterator set to the begining of the PackedList
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
