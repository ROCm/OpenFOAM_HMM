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
    Foam::UList

Description
    A 1D vector of objects of type \<T\>, where the size of the vector is
    known and can be used for subscript bounds checking, etc.

    Storage is not allocated during construction or use but is supplied to
    the constructor as an argument.  This type of list is particularly useful
    for lists that refer to parts of existing lists such as SubList.

SourceFiles
    UList.C
    UListI.H
    UListIO.C

\*---------------------------------------------------------------------------*/

#ifndef UList_H
#define UList_H

#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "nullObject.H"
#include "zero.H"
#include "stdFoam.H"
#include "Swap.H"

#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class labelRange;
template<class T> class List;
template<class T> class SubList;

// Forward declaration of friend functions and operators
template<class T> class UList;
template<class T> Ostream& operator<<(Ostream&, const UList<T>&);
template<class T> Istream& operator>>(Istream&, UList<T>&);

// Commonly required list types
typedef UList<char> charUList;
typedef UList<label> labelUList;

/*---------------------------------------------------------------------------*\
                           Class UList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UList
{
    // Private data

        //- Number of elements in UList
        label size_;

        //- Vector of values of type T
        T* __restrict__ v_;


    // Private Member Functions

        //- Disallow default shallow-copy assignment
        //
        //  Assignment of UList<T> may need to be either shallow (copy pointer)
        //  or deep (copy elements) depending on context or the particular type
        //  of list derived from UList and it is confusing and prone to error
        //  for the default assignment to be either.  The solution is to
        //  disallow default assignment and provide separate 'shallowCopy' and
        //  'deepCopy' member functions
        void operator=(const UList<T>&) = delete;


protected:

    // Protected Member Functions

        //- Override size to be inconsistent with allocated storage.
        //  Use with care
        inline void size(const label n);

        //- Write the UList with its compound type
        void writeEntry(Ostream& os) const;

        //- Return a validated (start,size) subset range, which means that it
        //- always addresses a valid section of the list.
        labelRange validateRange(const labelRange& range) const;

        //- Return a validated (start,size) subset range, which means that it
        //- always addresses a valid section of the list.
        labelRange validateRange
        (
            std::initializer_list<label> start_size_pair
        ) const;

public:

    // Related types

        //- Declare friendship with the List class
        friend class List<T>;

        //- Declare friendship with the SubList class
        friend class SubList<T>;


    // Static Member Functions

        //- Return a null UList
        inline static const UList<T>& null();


    // Public classes

        //- A list compare binary predicate for normal sort
        struct less
        {
            const UList<T>& values;

            less(const UList<T>& list)
            :
                values(list)
            {}

            bool operator()(const label a, const label b) const
            {
                return values[a] < values[b];
            }
        };

        //- A list compare binary predicate for reverse sort
        struct greater
        {
            const UList<T>& values;

            greater(const UList<T>& list)
            :
                values(list)
            {}

            bool operator()(const label a, const label b) const
            {
                return values[a] > values[b];
            }
        };


    // Constructors

        //- Null constructor
        inline UList();

        //- Construct from components
        inline UList(T* __restrict__ v, label size);


    // Member Functions

        // Access

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

            //- Return the binary size in number of characters of the UList
            //- if the element is a primitive type
            //  i.e. contiguous<T>() == true.
            //  Note that is of type streamsize since used in stream ops
            std::streamsize byteSize() const;


            //- Return a const pointer to the first data element.
            //  Similar to the STL front() method and the string::data() method
            //  This can be used (with caution) when interfacing with C code
            inline const T* cdata() const;

            //- Return a pointer to the first data element.
            //  Similar to the STL front() method and the string::data() method
            //  This can be used (with caution) when interfacing with C code
            inline T* data();

            //- Return the first element of the list
            inline T& first();

            //- Return first element of the list
            inline const T& first() const;

            //- Return the last element of the list
            inline T& last();

            //- Return the last element of the list
            inline const T& last() const;


        // Check

            //- Check start is within valid range (0 ... size-1)
            inline void checkStart(const label start) const;

            //- Check size is within valid range (0 ... size)
            inline void checkSize(const label size) const;

            //- Check index i is within valid range (0 ... size-1)
            inline void checkIndex(const label i) const;


      // Search

        //- Find index of the first occurence of the value.
        //  When start is specified, any occurences before start are ignored.
        //  Linear search.
        //  \return position in list or -1 if not found.
        label find(const T& val, const label start=0) const;

        //- Find index of the last occurence of the value.
        //  When pos is specified, any occurences after pos are ignored.
        //  Linear search.
        //  \return position in list or -1 if not found.
        label rfind(const T& val, const label pos=-1) const;

        //- True if the value if found in the list.
        //  When start is specified, any occurences before start are ignored.
        //  Linear search.
        //  \return true if found.
        inline bool found(const T& val, const label start=0) const;


      // Edit

        //- Move element to the first position.
        void moveFirst(const label i);

        //- Move element to the last position.
        void moveLast(const label i);

        //- Swap element with the first element. Fatal on an empty list.
        void swapFirst(const label i);

        //- Swap element with the last element. Fatal on an empty list.
        void swapLast(const label i);


      // Copy

        //- Copy the pointer held by the given UList
        inline void shallowCopy(const UList<T>& a);

        //- Copy elements of the given UList
        void deepCopy(const UList<T>& a);


    // Member operators

        //- Return element of UList
        inline T& operator[](const label i);

        //- Return element of constant UList
        //  Note that the bool specialization adds lazy evaluation so reading
        //  an out-of-range element returns false without any ill-effects
        inline const T& operator[](const label i) const;

        //- Return (start,size) subset from UList with non-const access.
        //  The range is subsetted with the list size itself to ensure that the
        //  result always addresses a valid section of the list.
        UList<T> operator[](const labelRange& range);

        //- Return (start,size) subset from UList with const access.
        //  The range is subsetted with the list size itself to ensure that the
        //  result always addresses a valid section of the list.
        const UList<T> operator[](const labelRange& range) const;

        //- Return (start,size) subset from UList with non-const access.
        //  The range is subsetted with the list size itself to ensure that the
        //  result always addresses a valid section of the list.
        UList<T> operator[](std::initializer_list<label> start_size);

        //- Return (start,size) subset from UList with const access.
        //  The range is subsetted with the list size itself to ensure that the
        //  result always addresses a valid section of the list.
        const UList<T> operator[]
        (
            std::initializer_list<label> start_size
        ) const;

        //- Allow cast to a const List<T>&
        inline operator const Foam::List<T>&() const;

        //- Assignment of all entries to the given value
        void operator=(const T& val);

        //- Assignment of all entries to zero
        void operator=(const zero);


    // STL type definitions

        //- Type of values the UList contains
        typedef T value_type;

        //- The type used for storing into value_type objects
        typedef T& reference;

        //- The type used for reading from constant value_type objects.
        typedef const T& const_reference;

        //- The type that can represent the difference between any two
        //- UList iterator objects
        typedef label difference_type;

        //- The type that can represent the size of a UList
        typedef label size_type;


    // STL iterator

        //- Random access iterator for traversing UList
        typedef T* iterator;

        //- Return an iterator to begin traversing the UList
        inline iterator begin();

        //- Return an iterator to end traversing the UList
        inline iterator end();


    // STL const_iterator

        //- Random access iterator for traversing UList
        typedef const T* const_iterator;

        //- Return const_iterator to begin traversing the constant UList
        inline const_iterator cbegin() const;

        //- Return const_iterator to end traversing the constant UList
        inline const_iterator cend() const;

        //- Return const_iterator to begin traversing the constant UList
        inline const_iterator begin() const;

        //- Return const_iterator to end traversing the constant UList
        inline const_iterator end() const;

    // Reverse iterators

        //- Generic const/non-const reverse iterator
        template<bool Const>
        class reverse_iterator_base
        {
        public:
            //- The const/non-const type for entries
            typedef typename std::conditional
                <Const, const T, T>::type value_type;

            //- A pointer to a const/non-const entry
            typedef value_type* pointer;

            //- A reference to a const/non-const entry
            typedef value_type& reference;


        private:

            //- The element pointer
            pointer ptr_;


        public:

            //- Construct null or from list element pointer
            inline reverse_iterator_base(pointer ptr = nullptr)
            :
                ptr_(ptr)
            {}


            //- Copy construct
            inline reverse_iterator_base(const reverse_iterator_base& iter)
            :
                ptr_(iter.ptr_)
            {}


            //- Reverse increment
            inline void operator++()
            {
                --ptr_;
            }

            //- Reverse increment
            inline reverse_iterator_base operator++(int)
            {
                reverse_iterator_base old(*this);
                --ptr_;
                return old;
            }

            //- Reverse increase
            inline void operator+=(int n)
            {
                ptr_ -= n;
            }

            //- Dereference iterator
            reference operator*() const
            {
                return *ptr_;
            }

            //- Dereference iterator
            pointer operator->() const
            {
                return ptr_;
            }

            //- Equality
            bool operator==(const reverse_iterator_base& iter) const
            {
                return ptr_ == iter.ptr_;
            }

            //- inequality
            bool operator!=(const reverse_iterator_base& iter) const
            {
                return ptr_ != iter.ptr_;
            }
        };


        //- STL reverse_iterator
        typedef reverse_iterator_base<false> reverse_iterator;

        //- Return reverse_iterator to begin reverse traversing the UList
        inline reverse_iterator rbegin();

        //- Return reverse_iterator to end reverse traversing the UList
        inline reverse_iterator rend();


        //- STL const reverse iterator
        typedef reverse_iterator_base<true> const_reverse_iterator;

        //- Return const_reverse_iterator to begin reverse traversing the UList
        inline const_reverse_iterator crbegin() const;

        //- Return const_reverse_iterator to end reverse traversing the UList
        inline const_reverse_iterator crend() const;

        //- Return const_reverse_iterator to begin reverse traversing the UList
        inline const_reverse_iterator rbegin() const;

        //- Return const_reverse_iterator to end reverse traversing the UList
        inline const_reverse_iterator rend() const;


    // STL member functions

        //- Return the number of elements in the UList
        inline label size() const;

        //- Return size of the largest possible UList
        inline label max_size() const;

        //- Return true if the UList is empty (ie, size() is zero)
        inline bool empty() const;

        //- Swap content with another UList of the same type in constant time
        inline void swap(UList<T>& lst);


    // STL member operators

        //- Equality operation on ULists of the same type.
        //  Returns true when the ULists are element-wise equal
        //  (using UList::value_type::operator==).  Takes linear time
        bool operator==(const UList<T>& a) const;

        //- The opposite of the equality operation. Takes linear time
        bool operator!=(const UList<T>& a) const;

        //- Compare two ULists lexicographically. Takes linear time
        bool operator<(const UList<T>& a) const;

        //- Compare two ULists lexicographically. Takes linear time
        bool operator>(const UList<T>& a) const;

        //- Return true if !(a > b). Takes linear time
        bool operator<=(const UList<T>& a) const;

        //- Return true if !(a < b). Takes linear time
        bool operator>=(const UList<T>& a) const;


    // Writing

        //- Write the List as a dictionary entry with keyword
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write the List, with line-breaks in ASCII if the list length
        //- exceeds shortListLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortListLen=0) const;


    // IOstream operators

        //- Write List to Ostream, as per writeList() with shortListLen=10
        friend Ostream& operator<< <T>
        (
            Ostream& os,
            const UList<T>& lst
        );

        //- Read List contents from Istream.
        //  Requires size to have been set before
        friend Istream& operator>> <T>
        (
            Istream& os,
            UList<T>& L
        );
};


// Global Functions

template<class T>
void sort(UList<T>& a);

template<class T, class Compare>
void sort(UList<T>& a, const Compare& comp);

template<class T>
void stableSort(UList<T>& a);

template<class T, class Compare>
void stableSort(UList<T>& a, const Compare& comp);

template<class T>
void shuffle(UList<T>& a);

// Reverse the first n elements of the list
template<class T>
inline void reverse(UList<T>& lst, const label n);

// Reverse all the elements of the list
template<class T>
inline void reverse(UList<T>& lst);

// Exchange contents of lists - see UList::swap().
template<class T>
inline void Swap(UList<T>& a, UList<T>& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
