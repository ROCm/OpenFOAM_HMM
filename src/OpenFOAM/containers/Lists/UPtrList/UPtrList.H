/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::UPtrList

Description
    A list of pointers to objects of type \<T\>, without allocation/deallocation
    management of the pointers - this is to be done elsewhere.
    The operator[] returns a reference to the object, not the pointer.

See Also
    Foam::PtrList

SourceFiles
    UPtrListI.H
    UPtrList.C
    UPtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef UPtrList_H
#define UPtrList_H

#include "List.H"
#include "PtrList.H"
#include <iterator>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

template<class T> class PtrList;
template<class T> class UPtrList;

template<class T> Ostream& operator<<(Ostream& os, const UPtrList<T>& lst);


/*---------------------------------------------------------------------------*\
                          Class UPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UPtrList
{
protected:

    // Protected data

        List<T*> ptrs_;


public:

    // STL type definitions

        //- Type of values the list contains
        typedef T value_type;

        //- A non-const reference to the value_type
        typedef T& reference;

        //- A const reference to the value_type
        typedef const T& const_reference;

        //- Random-access iterator with non-const access
        class iterator;

        //- Random-access iterator with const access
        class const_iterator;


    // Constructors

        //- Construct null
        UPtrList() = default;

        //- Construct with specified size, each element initialized to nullptr
        explicit inline UPtrList(const label len);

        //- Construct from PtrList, copying addresses of each list element.
        //  The argument is non-const access since this UPtrList can be used
        //  to change its values.
        explicit inline UPtrList(PtrList<T>& lst);

        //- Construct from UList, taking addresses of each list element
        explicit UPtrList(UList<T>& lst);

        //- Construct as copy or re-use as specified
        inline UPtrList(UPtrList<T>& lst, bool reuse);

        //- Copy construct
        inline UPtrList(const UPtrList<T>& lst) = default;

        //- Move construct
        inline UPtrList(UPtrList<T>&& lst);


    // Member functions

    // Access

        //- Return the number of elements in the UPtrList
        inline label size() const;

        //- Return true if the UPtrList is empty (ie, size() is zero)
        inline bool empty() const;

        //- Return reference to the first element of the list
        inline T& first();

        //- Return reference to first element of the list
        inline const T& first() const;

        //- Return reference to the last element of the list
        inline T& last();

        //- Return reference to the last element of the list
        inline const T& last() const;


    // Edit

        //- Set list size to zero.
        inline void clear();

        //- Reset size of UPtrList.
        //  New entries are initialized to nullptr.
        inline void setSize(const label newLen);

        //- Reset size of UPtrList.
        //  New entries are initialized to nullptr.
        inline void resize(const label newLen);

        //- Append an element to the end of the list
        inline void append(T* ptr);

        //- Swap content
        inline void swap(UPtrList<T>& lst);

        //- Transfer contents into this list and annul the argument
        inline void transfer(UPtrList<T>& lst);

        //- Return true if element is set (ie, not a nullptr)
        inline bool set(const label i) const;

        //- Set element to given pointer and return old element (can be null).
        //  No checks on new element
        inline T* set(const label i, T* ptr);

        //- Reorder elements. Reordering must be unique (ie, shuffle).
        void reorder(const labelUList& oldToNew);


    // Member operators

        //- Return element const reference
        inline const T& operator[](const label i) const;

        //- Return element reference
        inline T& operator[](const label i);

        //- Return element const pointer
        inline const T* operator()(const label i) const;

        //- Copy assignment
        UPtrList<T>& operator=(const UPtrList<T>& lst) = default;

        //- Move assignment
        inline void operator=(UPtrList<T>&& lst);


    // Iterators

        //- Random-access iterator with non-const access
        class iterator
        {
            T** ptr_;

        public:

            using iterator_category = std::random_access_iterator_tag;
            using value_type = T;
            using difference_type = label;
            using pointer = T*;
            using reference = T&;
            friend class const_iterator;

            //- Construct for a given entry
            inline iterator(T** ptr);

            // Member operators

            inline bool operator==(const iterator& iter) const;
            inline bool operator!=(const iterator& iter) const;

            inline reference operator*() const;
            inline reference operator()() const;

            // Forward iteration
            inline iterator& operator++();
            inline iterator operator++(int);

            inline iterator& operator--();
            inline iterator operator--(int);

            // Random-access
            inline iterator& operator+=(difference_type n);
            inline iterator& operator-=(difference_type n);
            inline iterator operator+(difference_type n) const;
            inline iterator operator-(difference_type n) const;

            inline difference_type operator-(const iterator& iter) const;

            inline reference operator[](difference_type n) const;

            inline bool operator<(const iterator& iter) const;
            inline bool operator>(const iterator& iter) const;

            inline bool operator<=(const iterator& iter) const;
            inline bool operator>=(const iterator& iter) const;
        };


        //- Random-access iterator with const access
        class const_iterator
        {
            const T* const* ptr_;

        public:

            using iterator_category = std::random_access_iterator_tag;
            using value_type = const T;
            using difference_type = label;
            using pointer = const T*;
            using reference = const T&;

            //- Construct for a given entry
            inline const_iterator(const T* const* ptr);

            //- Copy construct from non-const iterator
            inline const_iterator(const iterator& iter);

            // Member operators

            inline bool operator==(const const_iterator& iter) const;
            inline bool operator!=(const const_iterator& iter) const;

            inline reference operator*() const;
            inline reference operator()() const;

            // Forward iteration
            inline const_iterator& operator++();
            inline const_iterator operator++(int);

            inline const_iterator& operator--();
            inline const_iterator operator--(int);

            // Random-access
            inline const_iterator& operator+=(difference_type n);
            inline const_iterator& operator-=(difference_type n);
            inline const_iterator operator+(difference_type n) const;
            inline const_iterator operator-(difference_type n) const;

            inline difference_type operator-(const const_iterator& iter) const;

            inline reference operator[](difference_type n) const;

            inline bool operator<(const const_iterator& iter) const;
            inline bool operator>(const const_iterator& iter) const;

            inline bool operator<=(const const_iterator& iter) const;
            inline bool operator>=(const const_iterator& iter) const;
        };


        //- Return an iterator to begin traversing the UPtrList
        inline iterator begin();

        //- Return an iterator to end traversing the UPtrList
        inline iterator end();

        //- Return an const_iterator to begin traversing the UPtrList
        inline const_iterator cbegin() const;

        //- Return an const_iterator to end traversing the UPtrList
        inline const_iterator cend() const;

        //- Return an const_iterator to begin traversing the UPtrList
        inline const_iterator begin() const;

        //- Return an const_iterator to end traversing the UPtrList
        inline const_iterator end() const;


    // IOstream operator

        //- Write UPtrList to Ostream
        friend Ostream& operator<< <T>(Ostream& os, const UPtrList<T>& lst);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UPtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPtrList.C"
    #include "UPtrListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
