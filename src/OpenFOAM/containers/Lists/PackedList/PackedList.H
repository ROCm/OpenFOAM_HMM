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
    Could make PackedBitRef an iterator and use for traversing as well.

SourceFiles
    PackedListI.H
    PackedList.C

\*---------------------------------------------------------------------------*/

#ifndef PackedList_H
#define PackedList_H

#include "labelList.H"
#include "List.H"

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
                        Class PackedBitRef Declaration
\*---------------------------------------------------------------------------*/

//- The PackedBitRef is used for PackedList
template <int nBits>
class PackedBitRef
{
private:

    // private data
    // const mutable List<unsigned int>& list_;

    unsigned int& elem_;

    const label startBit_;

public:

    // Static Public Members

        //- The max. number of bits that can be templated.
        //  Might someday be useful for a template assert.
        inline static unsigned int max_bits();

        //- The max. value for an entry, can also be used as the mask
        //  eg, ((1 << 2) - 1) yields 0b0011
        inline static unsigned int max_value();

        //- The number of entries per storage entry
        inline static unsigned int packing();


    // Constructors
    inline PackedBitRef
    (
        // const List<unsigned int>& lst
        const unsigned int& elem, label startBit)
    :
        // list_(lst),
        elem_(const_cast<unsigned int&>(elem)),
        startBit_(startBit)
    {}


    // Members

    // Assign value
    inline void operator=(const unsigned int val)
    {
        unsigned int shiftedMask = max_value() << startBit_;
        unsigned int shiftedVal  = (val & max_value()) << startBit_;

        elem_ = (elem_ & ~shiftedMask) | shiftedVal;
    }

    // return value
    inline unsigned int operator()() const
    {
        return ((elem_ >> startBit_) & max_value());
    }

    //- Conversion operator
    inline operator unsigned int () const
    {
        return ((elem_ >> startBit_) & max_value());
    }

    //- Conversion operator
    inline operator bool() const
    {
        return !!((elem_ >> startBit_) & max_value());
    }

};


/*---------------------------------------------------------------------------*\
                         Class PackedList Declaration
\*---------------------------------------------------------------------------*/

template <int nBits>
class PackedList
:
    private List<unsigned int>
{
    // Private data

        //- Number of nBits entries
        label size_;

    // Private Member Functions

        //- Calculate underlying list size
        inline static label storageSize(const label);

        //- Calculate element index and offset (start) bit within storage
        inline static label location(const label, label& offset);

        //- Check if value is representable in nBits
        inline static void checkValue(const unsigned int);

        //- Check index i is within valid range (0 ... size-1).
        inline void checkIndex(const label) const;

public:

    // Public data

        //- The max. number of bits that can be templated.
        //  Might someday be useful for a template assert.
        inline static unsigned int max_bits()
        {
            return PackedBitRef<nBits>::max_bits();
        }

        //- The max. value for an entry, can also be used as the mask
        //  eg, ((1 << 2) - 1) yields 0b0011
        inline static unsigned int max_value()
        {
            return PackedBitRef<nBits>::max_value();
        }

        //- The number of entries per storage entry
        inline static unsigned int packing()
        {
            return PackedBitRef<nBits>::packing();
        }


    // Constructors

        //- Null constructor
        inline PackedList();

        //- Construct with given size. Note: initializes list to 0.
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

        //- Number of packed elements.
        inline label size() const;

        //- Return true if the list is empty (i.e., if size() == 0).
        inline bool empty() const;

        //- Get value at index I.
        //  Does not auto-vivifies elements.
        inline unsigned int get(const label i) const;

        //- Set value at index I. Return true if value changed.
        //  Does not auto-vivifies elements.
        inline bool set(const label i, const unsigned int val);

        //- Underlying storage
        inline const List<unsigned int>& storage() const;

        //- Return as labelList
        labelList values() const;

        // Edit

        //- Alter the size of the underlying storage.
        //  The addressed size will be truncated if needed to fit, but will
        //  remain otherwise untouched.
        //  Use this or reserve() in combination with append().
        inline void setCapacity(const label);

        //- Reset size of List, optionally specify a value for new elements.
        inline void resize(const label, const unsigned int& val = 0);

        //- Reset size of List, optionally specify a value for new elements.
        inline void setSize(const label, const unsigned int& val = 0);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size.
        //  Optionally provide an initialization value for new elements.
        inline void reserve(const label, const unsigned int& val = 0);

        //- Clear the list, i.e. set size to zero.
        //- Does not adjust the underlying storage
        inline void clear();

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Shrink the allocated space to what is used.
        inline void shrink();

        //- Transfer the contents of the argument List into this List
        //  and annull the argument list.
        inline void transfer(PackedList<nBits>&);

        //- Transfer contents to the Xfer container
        inline Xfer<PackedList<nBits> > xfer();

    // Member operators

        //- Append a value at the end of the list. Return true if value changed.
        inline bool append(const unsigned int val);

        //- Get value at index I
        //  Does not auto-vivifies elements.
        inline unsigned int operator[](const label i) const;

        //- Set value at index I.
        //  Returns proxy to perform the actual operation.
        //  Auto-vivifies any new values to zero.
        inline ::Foam::PackedBitRef<nBits> operator[](const label i);

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
