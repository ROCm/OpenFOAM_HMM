/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::FixedList

Description
    A 1D vector of objects of type \<T\> with a fixed length \<N\>.

SourceFiles
    FixedList.C
    FixedListI.H
    FixedListIO.C

\*---------------------------------------------------------------------------*/

#ifndef FixedList_H
#define FixedList_H

#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "zero.H"
#include "contiguous.H"
#include "autoPtr.H"
#include "Swap.H"
#include "HashFwd.H"
#include "SLListFwd.H"
#include "ListPolicy.H"

#include <initializer_list>
#include <iterator>
#include <type_traits>
#include <limits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class T, unsigned N> class FixedList;
template<class T> class UList;

template<class T, unsigned N>
Istream& operator>>(Istream& is, FixedList<T, N>& list);


/*---------------------------------------------------------------------------*\
                           Class FixedList Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned N>
class FixedList
{
    static_assert
    (
        N && N <= std::numeric_limits<int>::max(),
        "Size must be positive (non-zero) and fit as a signed int value"
    );

    // Private Data

        //- Vector of values of type T of length N.
        T v_[N];


protected:

    // Protected Member Functions

        //- Write the FixedList with its compound type
        void writeEntry(Ostream& os) const;


public:

    // STL Type Definitions

        //- The value type the FixedList contains
        typedef T value_type;

        //- The pointer type for non-const access to value_type items
        typedef T* pointer;

        //- The pointer type for const access to value_type items
        typedef const T* const_pointer;

        //- The type used for storing into value_type objects
        typedef T& reference;

        //- The type used for reading from constant value_type objects.
        typedef const T& const_reference;

        //- Random access iterator for traversing FixedList
        typedef T* iterator;

        //- Random access iterator for traversing FixedList
        typedef const T* const_iterator;

        //- The type to represent the size of a FixedList
        typedef label size_type;

        //- The difference between iterator objects
        typedef label difference_type;

        //- Reverse iterator (non-const access)
        typedef std::reverse_iterator<iterator> reverse_iterator;

        //- Reverse iterator (const access)
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;


    // Static Member Functions

        //- Return a null FixedList
        inline static const FixedList<T, N>& null();


    // Constructors

        //- Null constructor
        inline FixedList() = default;

        //- Construct and initialize all entries to given value
        inline explicit FixedList(const T& val);

        //- Construct and initialize all entries to zero
        inline explicit FixedList(const zero);

        //- Copy construct from C-array
        inline explicit FixedList(const T list[N]);

        //- Copy constructor
        inline FixedList(const FixedList<T, N>& list);

        //- Move construct by using move assignment for the individual
        //- list elements
        inline FixedList(FixedList<T, N>&& list);

        //- Construct given begin/end iterators
        //  Uses std::distance when verifying the size.
        template<class InputIterator>
        inline FixedList(InputIterator begIter, InputIterator endIter);

        //- Construct from an initializer list
        inline FixedList(std::initializer_list<T> list);

        //- Construct from UList
        inline explicit FixedList(const UList<T>& list);

        //- Construct from SLList
        inline explicit FixedList(const SLList<T>& list);

        //- Construct from Istream
        explicit FixedList(Istream& is);

        //- Clone
        inline autoPtr<FixedList<T, N>> clone() const;


    // Member Functions

    // Access

        //- Return a const pointer to the first data element.
        //  Similar to the STL front() method and the string::data() method
        //  This can be used (with caution) when interfacing with C code
        inline const T* cdata() const noexcept;

        //- Return a pointer to the first data element.
        //  Similar to the STL front() method and the string::data() method
        //  This can be used (with caution) when interfacing with C code
        inline T* data() noexcept;

        //- The first element of the list, position [0]
        inline T& first() noexcept;

        //- The first element of the list, position [0]
        inline const T& first() const noexcept;

        //- The last element of the list, position [N-1]
        inline T& last() noexcept;

        //- The last element of the list, position [N-1]
        inline const T& last() const noexcept;


        //- Return the forward circular index, i.e. next index
        //- which returns to the first at the end of the list
        inline label fcIndex(const label i) const;

        //- Return forward circular value (ie, next value in the list)
        inline const T& fcValue(const label i) const;

        //- Return forward circular value (ie, next value in the list)
        inline T& fcValue(const label i);

        //- Return the reverse circular index, i.e. previous index
        //- which returns to the last at the beginning of the list
        inline label rcIndex(const label i) const;

        //- Return reverse circular value (ie, previous value in the list)
        inline const T& rcValue(const label i) const;

        //- Return reverse circular value (ie, previous value in the list)
        inline T& rcValue(const label i);


    // Check

        //- Check start is within valid range [0,size)
        inline void checkStart(const label start) const;

        //- Check size is identical to template parameter N
        inline void checkSize(const label size) const;

        //- Check index is within valid range [0,N)
        inline void checkIndex(const label i) const;

        //- True if all entries have identical values, and list is non-empty
        inline bool uniform() const;


    // Search

        //- Find index of the first occurence of the value.
        //  Linear search.
        //  \return -1 if not found.
        label find(const T& val, const label start=0) const;

        //- True if the value if found in the list. Linear search.
        inline bool found(const T& val, const label start=0) const;


    // Edit

        //- Dummy resize function, to make FixedList consistent with List
        inline void resize(const label n);

        //- Dummy setSize function, to make FixedList consistent with List
        inline void setSize(const label n);

        //- Move element to the first position.
        void moveFirst(const label i);

        //- Move element to the last position.
        void moveLast(const label i);

        //- Swap element with the first element.
        void swapFirst(const label i);

        //- Swap element with the last element.
        void swapLast(const label i);

        //- Transfer by swapping using a move assignment for the content
        //- of the individual list elements
        inline void transfer(FixedList<T, N>& list);


    // Member Operators

        //- Return element of FixedList
        inline T& operator[](const label i);

        //- Return element of constant FixedList
        inline const T& operator[](const label i) const;

        //- Assignment to array operator. Takes linear time
        inline void operator=(const T list[N]);

        //- Assignment to UList operator. Takes linear time
        inline void operator=(const UList<T>& list);

        //- Assignment to SLList operator. Takes linear time
        inline void operator=(const SLList<T>& list);

        //- Assignment to an initializer list. Takes linear time
        inline void operator=(std::initializer_list<T> list);

        //- Assignment of all entries to the given value
        inline void operator=(const T& val);

        //- Copy assignment
        inline void operator=(const FixedList<T, N>& list);

        //- Move assignment
        inline void operator=(FixedList<T, N>&& list);


    // Random access iterator (non-const)

        //- Return an iterator to begin traversing the FixedList
        inline iterator begin();

        //- Return an iterator to end traversing the FixedList
        inline iterator end();


    // Random access iterator (const)

        //- Return const_iterator to begin traversing the constant FixedList
        inline const_iterator cbegin() const;

        //- Return const_iterator to end traversing the constant FixedList
        inline const_iterator cend() const;

        //- Return const_iterator to begin traversing the constant FixedList
        inline const_iterator begin() const;

        //- Return const_iterator to end traversing the constant FixedList
        inline const_iterator end() const;


    // Reverse iterator (non-const)

        //- Return reverse_iterator to begin reverse traversing the FixedList
        inline reverse_iterator rbegin();

        //- Return reverse_iterator to end reverse traversing the FixedList
        inline reverse_iterator rend();


    // Reverse iterator (const)

        //- Return const_reverse_iterator to begin reverse traversing FixedList
        inline const_reverse_iterator crbegin() const;

        //- Return const_reverse_iterator to end reverse traversing FixedList
        inline const_reverse_iterator crend() const;

        //- Return const_reverse_iterator to begin reverse traversing FixedList
        inline const_reverse_iterator rbegin() const;

        //- Return const_reverse_iterator to end reverse traversing FixedList
        inline const_reverse_iterator rend() const;


    // STL Member Functions

        //- Always false since zero-sized FixedList is compile-time disabled.
        static constexpr bool empty() noexcept { return !N; }

        //- Return the number of elements in the FixedList
        static constexpr label size() noexcept { return N; }

        //- The dimensioned size (template parameter N) of the FixedList
        static constexpr unsigned max_size() noexcept { return N; }

        //- Swap lists by swapping the content of the individual list elements
        inline void swap(FixedList<T, N>& list);


    // STL Member Operators

        //- Equality operation on FixedLists of the same type.
        //  Returns true when the FixedLists are element-wise equal
        //  (using FixedList::value_type::operator==). Takes linear time
        bool operator==(const FixedList<T, N>& list) const;

        //- The opposite of the equality operation. Takes linear time
        bool operator!=(const FixedList<T, N>& list) const;

        //- Compare two FixedLists lexicographically. Takes linear time
        bool operator<(const FixedList<T, N>& list) const;

        //- Compare two FixedLists lexicographically. Takes linear time
        bool operator>(const FixedList<T, N>& list) const;

        //- Return true if !(a > b). Takes linear time
        bool operator<=(const FixedList<T, N>& list) const;

        //- Return true if !(a < b). Takes linear time
        bool operator>=(const FixedList<T, N>& list) const;


    // Writing

        //- Write the list as a dictionary entry with keyword
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write List, with line-breaks in ASCII when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const;


    // IOstream Operators

        //- Read from Istream, discarding contents of existing List
        friend Istream& operator>> <T, N>
        (
            Istream& is,
            FixedList<T, N>& list
        );


    // Hashing

        //- Hashing function class for FixedList
        //  Normally use the global Hash specialization, but can also use
        //  this one for inheritance in sub-classes
        template<class HashT=Foam::Hash<T>>
        struct Hash
        {
            inline unsigned operator()
            (
                const FixedList<T, N>& obj,
                unsigned seed=0
            ) const
            {
                if (is_contiguous<T>::value)
                {
                    return Hasher(obj.cdata(), N*sizeof(T), seed);
                }

                for (const T& val : obj)
                {
                    seed = HashT()(val, seed);
                }
                return seed;
            }
        };
};

// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- FixedList is contiguous if the type is contiguous
template<class T, unsigned N>
struct is_contiguous<FixedList<T, N>> : is_contiguous<T> {};

//- Check for FixedList of labels
template<class T, unsigned N>
struct is_contiguous_label<FixedList<T, N>> : is_contiguous_label<T> {};

//- Check for FixedList of scalars
template<class T, unsigned N>
struct is_contiguous_scalar<FixedList<T, N>> : is_contiguous_scalar<T> {};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Swap FixedList contents - see FixedList::swap().
//  Internally this actually swaps the individual list elements
template<class T, unsigned N>
inline void Swap(FixedList<T, N>& lhs, FixedList<T, N>& rhs);


//- Hashing for FixedList data, which uses Hasher for contiguous data and
//- element-wise incrementally hashing otherwise.
template<class T, unsigned N>
struct Hash<FixedList<T, N>>
{
    inline unsigned operator()
    (
        const FixedList<T, N>& obj,
        unsigned seed=0
    ) const
    {
        if (is_contiguous<T>::value)
        {
            return Hasher(obj.cdata(), N*sizeof(T), seed);
        }

        for (const T& val : obj)
        {
            seed = Hash<T>()(val, seed);
        }
        return seed;
    }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Write List to Ostream, as per FixedList::writeList() with default length.
//  The default short-length is given by Detail::ListPolicy::short_length
template<class T, unsigned N>
Ostream& operator<<(Ostream& os, const FixedList<T, N>& list)
{
    return list.writeList(os, Detail::ListPolicy::short_length<T>::value);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FixedListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FixedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
