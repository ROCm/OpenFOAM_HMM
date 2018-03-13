/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::PackedList

Description
    A dynamically list of packed unsigned integers, with the number of bits
    per item specified by the \<nBits\> template parameter.

    Resizing is similar to DynamicList so that clear() and resize() affect
    the addessed size, but not the allocated size. The reserve() and
    setCapacity() methods can be used to influence the allocation.

Note
    In a const context, the '[]' operator simply returns the stored value,
    with out-of-range elements returned as zero.

    In a non-const context, the '[]' operator returns a reference to an
    existing value. Some caution with out-of-range elements to ensure
    that the const version of the [] operator is being called.
    The get() method is functionally identical the the '[]' operator, but
    is always const access.

    The set() and unset() methods return a bool if the value changed.
    This can be useful for branching on changed values.

    \code
        list.set(5, 4);
        changed = list.set(5, 8);
        if (changed) ...
    \endcode

    In a const context, reading an out-of-range element returns zero without
    affecting the list size.
    For example,
    \code
        list.resize(4);
        Info<< list.get(10) << "\n";  // print zero, but doesn't adjust list
    \endcode

    Also note that all unused internal storage elements are guaranteed to
    always be bit-wise zero. This property must not be violated by any
    inheriting classes.

See also
    Foam::DynamicList

SourceFiles
    PackedListI.H
    PackedList.C

\*---------------------------------------------------------------------------*/

#ifndef PackedList_H
#define PackedList_H

#include "labelList.H"
#include "UIndirectList.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<unsigned nBits> class PackedList;
class Istream;
class Ostream;

template<unsigned nBits>
Istream& operator>>(Istream& is, PackedList<nBits>& lst);
template<unsigned nBits>
Ostream& operator<<(Ostream& os, const PackedList<nBits>& lst);


/*---------------------------------------------------------------------------*\
                       Class PackedListCore Declaration
\*---------------------------------------------------------------------------*/

//- Template-invariant bits for PackedList
struct PackedListCore
{
    //- Construct null
    PackedListCore()
    {}

    //- Define template name and debug
    ClassName("PackedList");
};


/*---------------------------------------------------------------------------*\
                         Class PackedList Declaration
\*---------------------------------------------------------------------------*/

template<unsigned nBits=1>
class PackedList
:
    public PackedListCore,
    private List<unsigned int>
{
protected:

    typedef unsigned int      StorageType;
    typedef List<StorageType> StorageList;

    // Protected Member Functions

        //- Calculate the list length when packed
        inline static constexpr label packedLength(const label nElem);

        //- Read a list entry (allows for specialization)
        inline static unsigned int readValue(Istream& is);

        //- Read an index/value pair and set accordingly.
        //  For bool specialization, read a single index value
        inline void setPair(Istream& is);


private:

    //- nBits must be positive (non-zero) and fit within the storage.
    //  For efficiency, however, require packing at least 2 items otherwise
    //  it is more efficient to use a normal list.
    //  Thus max nBits is 1/2 of the base storage size.
    //  For simplicity, assume 8-bit bytes in the assert.
    static_assert
    (
        nBits && nBits <= (sizeof(StorageType) << 2),
        "nBits must be positive (non-zero) and fit within the storage"
    );

    // Private data

        //- Number of nBits entries
        label size_;


protected:

    // Protected Member Functions

        //- Write as a dictionary entry
        void writeEntry(Ostream& os) const;


public:

    // Forward declaration of access classes

        class reference;
        typedef unsigned int const_reference;


    // Public data

        //- The max. number of bits that can be templated.
        //  Might someday be useful for a template assert.
        inline static constexpr unsigned int max_bits();

        //- The max. value for an entry, which simultaneously the bit-mask
        //  eg, ((1 << 2) - 1) yields 0b0011
        inline static constexpr unsigned int max_value();

        //- The number of entries per packed storage element
        inline static constexpr unsigned int packing();

        //- Masking for all bits below the offset
        inline static constexpr unsigned int maskLower(unsigned offset);


    // Constructors

        //- Null constructor
        inline PackedList();

        //- Construct with given size, initializes list to 0
        explicit inline PackedList(const label size);

        //- Construct with given size and value for all elements
        inline PackedList(const label size, const unsigned val);

        //- Construct from Istream
        inline PackedList(Istream& is);

        //- Copy construct
        inline PackedList(const PackedList<nBits>& lst);

        //- Move construct
        inline PackedList(PackedList<nBits>&& lst);

        //- Construct from a list of labels
        explicit inline PackedList(const labelUList& lst);

        //- Construct from an indirect list of labels
        explicit inline PackedList(const labelUIndList& lst);

        //- Clone
        inline autoPtr<PackedList<nBits>> clone() const;


    // Member Functions

    // Access

        //- The number of elements that can be stored before reallocating
        inline label capacity() const;

        //- Number of entries.
        inline label size() const;

        //- Return true if the list is empty (ie, size() is zero).
        inline bool empty() const;

        //- Get value at index I.
        //  Never auto-vivify entries.
        inline unsigned int get(const label i) const;

        //- Set value at index I, default value set is the max_value.
        //  Does auto-vivify for non-existent, non-zero entries.
        //  \return true if value changed.
        inline bool set(const label i, const unsigned int val = ~0u);

        //- Unset the entry at index I.
        //  Never auto-vivify entries.
        //  \return true if value changed.
        inline bool unset(const label i);

        //- Return the underlying packed storage
        //  Manipulate with utmost caution
        inline List<unsigned int>& storage();

        //- Return the underlying packed storage
        inline const List<unsigned int>& storage() const;

        //- The list length when packed
        inline label packedLength() const;

        //- The number of bytes used in the underlying storage
        inline std::streamsize byteSize() const;

        //- Count number of bits set, O(log(n))
        //  Uses the Hamming weight (population count) method
        //  http://en.wikipedia.org/wiki/Hamming_weight
        unsigned int count() const;

        //- Return the values as a list of labels
        labelList values() const;

        //- Print bit patterns, optionally output unused elements
        //
        // addressable bits:
        //   on: '1', off: '-'
        //
        // non-addressable bits:
        //   on: '!', off: '.'
        //
        Ostream& printBits(Ostream& os, const bool fullOutput=false) const;

        //- Print information and bit patterns (with printBits)
        Ostream& printInfo(Ostream& os, const bool fullOutput=false) const;

    // Check

        //- Check index is within valid range [0,size)
        inline void checkIndex(const label i) const;


    // Edit

        //- Trim any trailing zero elements
        bool trim();

        //- Invert the bits in the addressable region
        void flip();

        //- Clear all bits
        inline void reset();

        //- Alter the size of the underlying storage.
        //  The addressed size will be truncated if needed to fit, but will
        //  remain otherwise untouched.
        inline void setCapacity(const label nElem);

        //- Reset addressable list size, does not shrink the allocated size.
        //  Optionally specify a value for new elements.
        inline void resize(const label nElem, const unsigned int val = 0u);

        //- Alias for resize()
        inline void setSize(const label nElem, const unsigned int val = 0u);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size.
        //  The list size is adjusted as per DynamicList with
        //  SizeInc=0, SizeMult=2, SizeDiv=1
        inline void reserve(const label nElem);

        //- Clear the list, i.e. set addressable size to zero.
        //  Does not adjust the underlying storage
        inline void clear();

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Shrink the allocated space to what is actually used.
        inline void shrink();

        //- Transfer the contents of the argument list into this list
        //  and annul the argument list.
        inline void transfer(PackedList<nBits>& lst);


    // IO

        //- Clear list and read from stream
        Istream& read(Istream& is);

        //- Write the List, with line-breaks in ASCII if the list length
        //- exceeds shortListLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortListLen=0) const;

        //- Write as a dictionary entry with keyword
        void writeEntry(const word& keyword, Ostream& os) const;


    // Member Operators

        //- Append a value at the end of the list
        inline PackedList<nBits>& append(const unsigned int val);

        //- Remove and return the last element
        inline unsigned int remove();

        //- Identical to get() - get value at index.
        //  Never auto-vivify entries.
        inline unsigned int operator[](const label i) const;

        //- Non-const access to value at index.
        //  Fatal for out-of-range indices
        inline reference operator[](const label i);

        //- Assignment of all entries to the given value. Takes linear time.
        inline void operator=(const unsigned int val);

        //- Copy assignment.
        void operator=(const PackedList<nBits>& lst);

        //- Move assignment.
        void operator=(PackedList<nBits>&& lst);

        //- Assignment operator.
        void operator=(const labelUList& lst);

        //- Assignment operator.
        void operator=(const labelUIndList& lst);


    // Iterators and helpers

        //- A reference for read/write access to an entry
        class reference
        {
            friend class PackedList;    // Access for PackedList
            void operator&() = delete;  // Refuse to provide its address

            unsigned int& ref_;   //!< Reference to underlying block value
            unsigned int shift_;  //!< The bit shift for the sub-portion

            //- Construct by taking reference of block from within the list
            //- at the specified index.
            inline reference(PackedList* lst, const label index);

            //- Get value as unsigned, no range-checking
            inline unsigned int get() const;

            //- Set value, returning true if changed, no range-checking
            inline bool set(unsigned int val);

        public:

            //- Value assignment
            inline void operator=(const reference& other);

            //- Value assignment
            inline void operator=(const unsigned int val);

            //- Conversion operator
            inline operator unsigned int () const;
        };


    // IOstream Operators

        friend Istream& operator>> <nBits>
        (
            Istream& is,
            PackedList<nBits>& lst
        );

        friend Ostream& operator<< <nBits>
        (
            Ostream& os,
            const PackedList<nBits>& lst
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PackedListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PackedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
