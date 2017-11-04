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
    Foam::FixedList

Description
    A 1D vector of objects of type \<T\> with a fixed size \<Size\>.

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
#include "Hash.H"
#include "autoPtr.H"
#include "Swap.H"

#include <type_traits>
#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned Size> class FixedList;

template<class T, unsigned Size>
Istream& operator>>(Istream&, FixedList<T, Size>&);

template<class T, unsigned Size>
Ostream& operator<<(Ostream&, const FixedList<T, Size>&);

template<class T> class UList;

class SLListBase;
template<class LListBase, class T> class LList;
template<class T>
using SLList = LList<SLListBase, T>;

/*---------------------------------------------------------------------------*\
                           Class FixedList Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned Size>
class FixedList
{
    static_assert
    (
        Size && Size <= INT_MAX,
        "Size must be positive (non-zero) and also fit as a signed value"
    );

    // Private data

        //- Vector of values of type T of size Size.
        T v_[Size];


protected:

    // Protected Member Functions

        //- Write the FixedList with its compound type
        void writeEntry(Ostream& os) const;


public:

    //- Hashing function class.
    //  Use Hasher directly for contiguous data. Otherwise hash incrementally.
    template<class HashT=Hash<T>>
    class Hash
    {
    public:
        Hash()
        {}

        inline unsigned operator()
        (
            const FixedList<T, Size>&,
            unsigned seed = 0
        ) const;
    };


    // Static Member Functions

        //- Return a null FixedList
        inline static const FixedList<T, Size>& null();


    // Constructors

        //- Null constructor
        inline FixedList();

        //- Construct from value
        explicit inline FixedList(const T& val);

        //- Copy construct from C-array
        explicit inline FixedList(const T lst[Size]);

        //- Copy constructor
        inline FixedList(const FixedList<T, Size>& lst);

        //- Move constructor
        inline FixedList(FixedList<T, Size>&& lst);

        //- Construct given begin/end iterators
        //  Uses std::distance when verifying the size.
        template<class InputIterator>
        inline FixedList(InputIterator begIter, InputIterator endIter);

        //- Construct from an initializer list
        inline FixedList(std::initializer_list<T> lst);

        //- Construct from UList
        explicit inline FixedList(const UList<T>& lst);

        //- Construct from SLList
        explicit inline FixedList(const SLList<T>& lst);

        //- Construct from Istream
        FixedList(Istream& is);

        //- Clone
        inline autoPtr<FixedList<T, Size>> clone() const;


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
        //  Linear search.
        //  \return -1 if not found.
        label find(const T& val, const label start=0) const;

        //- True if the value if found in the list. Linear search.
        inline bool found(const T& val, const label start=0) const;


        // Edit

        //- Dummy resize function
        //  needed to make FixedList consistent with List
        inline void resize(const label s);

        //- Dummy setSize function
        //  needed to make FixedList consistent with List
        inline void setSize(const label s);

        //- Move element to the first position.
        void moveFirst(const label i);

        //- Move element to the last position.
        void moveLast(const label i);

        //- Swap element with the first element.
        void swapFirst(const label i);

        //- Swap element with the last element.
        void swapLast(const label i);

        //- Copy (not transfer) the argument contents
        //  needed to make FixedList consistent with List
        void transfer(const FixedList<T, Size>& lst);


    // Member operators

        //- Return element of FixedList
        inline T& operator[](const label i);

        //- Return element of constant FixedList
        inline const T& operator[](const label i) const;

        //- Assignment to array operator. Takes linear time
        inline void operator=(const T lst[Size]);

        //- Assignment to UList operator. Takes linear time
        inline void operator=(const UList<T>& lst);

        //- Assignment to SLList operator. Takes linear time
        inline void operator=(const SLList<T>& lst);

        //- Assignment to an initializer list. Takes linear time
        inline void operator=(std::initializer_list<T> lst);

        //- Assignment of all entries to the given value
        inline void operator=(const T& val);

        //- Copy assignment
        inline void operator=(const FixedList<T, Size>& lst);

        //- Move assignment
        inline void operator=(FixedList<T, Size>&& lst);


    // STL type definitions

        //- Type of values the FixedList contains
        typedef T value_type;

        //- The type used for storing into value_type objects
        typedef T& reference;

        //- The type used for reading from constant value_type objects.
        typedef const T& const_reference;

        //- The type that can represent the difference between any two
        //- FixedList iterator objects
        typedef label difference_type;

        //- The type that can represent the size of a FixedList
        typedef label size_type;


    // STL iterator

        //- Random access iterator for traversing FixedList
        typedef T* iterator;

        //- Return an iterator to begin traversing the FixedList
        inline iterator begin();

        //- Return an iterator to end traversing the FixedList
        inline iterator end();


    // STL const_iterator

        //- Random access iterator for traversing FixedList
        typedef const T* const_iterator;

        //- Return const_iterator to begin traversing the constant FixedList
        inline const_iterator cbegin() const;

        //- Return const_iterator to end traversing the constant FixedList
        inline const_iterator cend() const;

        //- Return const_iterator to begin traversing the constant FixedList
        inline const_iterator begin() const;

        //- Return const_iterator to end traversing the constant FixedList
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

        //- Return reverse_iterator to begin reverse traversing the FixedList
        inline reverse_iterator rbegin();

        //- Return reverse_iterator to end reverse traversing the FixedList
        inline reverse_iterator rend();


        //- STL const reverse iterator
        typedef reverse_iterator_base<true> const_reverse_iterator;

        //- Return const_reverse_iterator to begin reverse traversing FixedList
        inline const_reverse_iterator crbegin() const;

        //- Return const_reverse_iterator to end reverse traversing FixedList
        inline const_reverse_iterator crend() const;

        //- Return const_reverse_iterator to begin reverse traversing FixedList
        inline const_reverse_iterator rbegin() const;

        //- Return const_reverse_iterator to end reverse traversing FixedList
        inline const_reverse_iterator rend() const;


    // STL member functions

        //- Return the number of elements in the FixedList
        inline label size() const;

        //- Return size of the largest possible FixedList
        inline label max_size() const;

        //- Always false since zero-sized FixedList is compile-time disabled.
        inline bool empty() const;

        //- Swap content with another FixedList of the same type.
        void swap(FixedList<T, Size>& lst);


    // STL member operators

        //- Equality operation on FixedLists of the same type.
        //  Returns true when the FixedLists are element-wise equal
        //  (using FixedList::value_type::operator==).  Takes linear time
        bool operator==(const FixedList<T, Size>& a) const;

        //- The opposite of the equality operation. Takes linear time
        bool operator!=(const FixedList<T, Size>& a) const;

        //- Compare two FixedLists lexicographically. Takes linear time
        bool operator<(const FixedList<T, Size>& a) const;

        //- Compare two FixedLists lexicographically. Takes linear time
        bool operator>(const FixedList<T, Size>& a) const;

        //- Return true if !(a > b). Takes linear time
        bool operator<=(const FixedList<T, Size>& a) const;

        //- Return true if !(a < b). Takes linear time
        bool operator>=(const FixedList<T, Size>& a) const;


    // Writing

        //- Write the List as a dictionary entry with keyword
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write the List, with line-breaks in ASCII if the list length
        //- exceeds shortListLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortListLen=0) const;


    // IOstream operators

        //- Read List from Istream, discarding contents of existing List
        friend Istream& operator>> <T, Size>
        (
            Istream& is,
            FixedList<T, Size>& L
        );

        //- Write List to Ostream, as per writeList() with shortListLen=10
        friend Ostream& operator<< <T, Size>
        (
            Ostream& os,
            const FixedList<T, Size>& L
        );
};


// Global Functions

// Exchange contents of lists - see FixedList::swap().
template<class T, unsigned Size>
inline void Swap(FixedList<T,Size>& a, FixedList<T,Size>& b);


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
