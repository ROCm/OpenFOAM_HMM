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
    List of packed unsigned ints.

    Gets given the number of bits per item.

ToDo
    Add checks for bad template parameters (ie, nBits=0, nBits too large).
    Could make PackedBitRef an iterator and use for traversing as well.
    It could be useful to make PackedList behave a bit like DynamicList.

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

/*---------------------------------------------------------------------------*\
                      Class PackedListName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(PackedList);


/*---------------------------------------------------------------------------*\
                        Class PackedBitRef Declaration
\*---------------------------------------------------------------------------*/

//- The PackedBitRef is used for PackedList
class PackedBitRef
{
private:

    // private data
    unsigned int& elem_;

    const label startBit_;

    const unsigned int mask_;

public:

    inline PackedBitRef(unsigned int& elem, label startBit, unsigned int mask)
    :
        elem_(elem),
        startBit_(startBit),
        mask_(mask)
    {}

    inline void operator=(const unsigned int val)
    {
        unsigned int shiftedMask = mask_ << startBit_;
        unsigned int shiftedVal  = (val & mask_) << startBit_;

        elem_ = (elem_ & ~shiftedMask) | shiftedVal;
    }

    inline operator unsigned int () const
    {
        return ((elem_ >> startBit_) & mask_);
    }

    inline operator bool() const
    {
        return !!((elem_ >> startBit_) & mask_);
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
        inline static unsigned int max_bits();

        //- The max. value for an entry, can also be used as the mask
        //  eg, ((1 << 2) - 1) yields 0b0011
        inline static unsigned int max_value();

        //- The number of entries per storage entry
        inline static unsigned int packing();


    // Constructors

        //- Null constructor
        inline PackedList();

        //- Construct with given size. Note: initializes intList to 0.
        inline PackedList(const label size);

        //- Construct with given size and value for all elements.
        PackedList(const label size, const unsigned val);

        //- Copy constructor.
        PackedList(const PackedList<nBits>& PList);

        //- Construct by transferring the parameter contents
        PackedList(const Xfer<PackedList<nBits> >&);

        //- Construct from a list of labels
        PackedList(const UList<label>&);

        //- Clone
        inline autoPtr<PackedList<nBits> > clone() const;

    // Member Functions

        // Access

            //- The number of elements that can be stored before resizing
            inline label capacity() const;

            //- Number of packed elements
            inline label size() const;

            //- Return true if the list is empty (i.e., if size() == 0).
            inline bool empty() const;

            //- Get value at index I
            inline unsigned int get(const label i) const;

            //- Set value at index I. Return true if value changed.
            inline bool set(const label i, const unsigned int val);

            //- Underlying storage
            inline List<unsigned int>& storage();

            //- Underlying storage
            inline const List<unsigned int>& storage() const;


        // Edit

            //- Reset size of List, setting zero for any new elements.
            void setSize(const label);

            //- Reset size of List and value for new elements.
            void setSize(const label, const unsigned int& val);

            //- Reset size of List, setting zero for any new elements.
            inline void resize(const label);

            //- Reset size of List and value for new elements.
            inline void resize(const label, const unsigned int& val);

        //- Construct with given size and value for all elements.

            //- Clear the list, i.e. set size to zero.
            void clear();

            //- Transfer the contents of the argument List into this List
            //  and annull the argument list.
            void transfer(PackedList<nBits>&);

            //- Transfer contents to the Xfer container
            inline Xfer<PackedList<nBits> > xfer();


    // Member operators

            //- Get value at index i
            inline unsigned int operator[](const label i) const;

            //- Set value at index i.
            //  Returns proxy to perform the actual operation
            inline ::Foam::PackedBitRef operator[](const label i);

            //- Assignment operator. Takes linear time.
            void operator=(const PackedList<nBits>&);

            //- Assignment of all entries to the given value.
            //  Does set on all elements.
            inline void operator=(const unsigned int val);

            //- Return as labelList
            labelList operator()() const;


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
