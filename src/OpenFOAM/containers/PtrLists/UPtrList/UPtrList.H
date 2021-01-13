/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

Note
    The class definition is such that it contains a list of pointers, but
    itself does not inherit from a list of pointers since this would
    wreak havoc later with inheritance resolution.

See Also
    Foam::PtrList
    Foam::PtrDynList

SourceFiles
    UPtrListI.H
    UPtrList.C

\*---------------------------------------------------------------------------*/

#ifndef UPtrList_H
#define UPtrList_H

#include "PtrListDetail.H"
#include <iterator>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class PtrList;
template<class T> class UPtrList;
template<class T> Ostream& operator<<(Ostream& os, const UPtrList<T>& list);

/*---------------------------------------------------------------------------*\
                          Class UPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UPtrList
{
protected:

    // Protected Member Data

        //- The list of pointers
        Detail::PtrListDetail<T> ptrs_;


    // Protected Member Functions

        //- Adjust addressable size
        inline void setAddressableSize(const label n) noexcept;


    // Constructors

        //- Low-level move construct
        inline explicit UPtrList(Detail::PtrListDetail<T>&& ptrs);


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

        //- Default construct
        inline constexpr UPtrList() noexcept;

        //- Construct with specified size, each element initialized to nullptr
        inline explicit UPtrList(const label len);

        //- Copy construct (shallow copies addresses)
        inline UPtrList(const UPtrList<T>& list);

        //- Move construct
        inline UPtrList(UPtrList<T>&& list);

        //- Construct as shallow copy or re-use as specified
        inline UPtrList(UPtrList<T>& list, bool reuse);

        //- Shallow copy from PtrList.
        //  The argument is non-const to reflect that the UPtrList can change
        //  the values (not the addresses) within the original list.
        explicit UPtrList(PtrList<T>& list);

        //- Construct from UList of pointers
        inline explicit UPtrList(const UList<T*>& list);

        //- Construct from UList, taking the address of each list element
        //  The argument is non-const to reflect that the UPtrList can change
        //  the values of the original list.
        inline explicit UPtrList(UList<T>& list);


    // Member Functions

    // Access

        //- The number of elements in the list
        inline label size() const noexcept;

        //- True if the list is empty (ie, size() is zero)
        inline bool empty() const noexcept;

        //- Return reference to the first element of the list
        inline T& first();

        //- Return reference to first element of the list
        inline const T& first() const;

        //- Return reference to the last element of the list
        inline T& last();

        //- Return reference to the last element of the list
        inline const T& last() const;

        //- Return pointer to element (can be nullptr),
        //- without bounds checking.
        inline T* get(const label i);

        //- Return const pointer to element (can be nullptr),
        //- without bounds checking.
        inline const T* get(const label i) const;

        //- Return const pointer to element (can be nullptr),
        //- without bounds checking - same as get().
        //  The return value can also be tested as a bool.
        const T* set(const label i) const { return this->get(i); }


    // Edit

        //- Set list size to zero.
        inline void clear();

        //- Change the size of the list.
        //  New entries are initialized to nullptr.
        inline void resize(const label newLen);

        //- Alias for resize()
        void setSize(const label n) { this->resize(n); }

        //- Squeeze out intermediate nullptr entries in the list of pointers
        //  \return the number of non-null entries
        label squeezeNull();

        //- Append an element to the end of the list
        inline void append(T* ptr);

        //- Swap content
        inline void swap(UPtrList<T>& list);

        //- Transfer contents into this list and annul the argument
        inline void transfer(UPtrList<T>& list);

        //- Set element to specified pointer and return the old list element,
        //- which can be a nullptr.
        //  No-op if the new pointer value is identical to the current content.
        inline T* set(const label i, T* ptr);

        //- Reorder elements.
        //- Reordering must be unique (ie, shuffle).
        //  Optionally verify that all pointers have been set.
        void reorder(const labelUList& oldToNew, const bool testNull = true);

        //- Reorder elements according to new order mapping (newToOld).
        //- Reordering must be unique (ie, shuffle).
        //  Optionally verify that all pointers have been set.
        void sortOrder(const labelUList& order, const bool testNull = true);


    // Member Operators

        //- Return const reference to the element
        inline const T& operator[](const label i) const;

        //- Return reference to the element
        inline T& operator[](const label i);

        //- Return const pointer to the element - same as get().
        inline const T* operator()(const label i) const;

        //- Copy assignment (shallow copies addresses)
        inline void operator=(const UPtrList<T>& list);

        //- Move assignment
        inline void operator=(UPtrList<T>&& list);


    // Iterators

        //- Random-access iterator with non-const access
        class iterator
        {
            // Pointer to parent
            T** ptr_;

        public:

            using iterator_category = std::random_access_iterator_tag;
            using value_type = T;
            using difference_type = label;
            using pointer = T*;
            using reference = T&;
            friend class const_iterator;

            //- Construct for a given entry
            inline iterator(T** ptr) noexcept;

            // Member Functions

            //- Return pointer, can be nullptr.
            inline pointer get() const;

            // Member Operators

            inline pointer operator->() const;
            inline reference operator*() const;

            reference operator()() const { return this->operator*(); }
            inline reference operator[](difference_type n) const;

            // Forward iteration
            inline iterator& operator++() noexcept;
            inline iterator operator++(int) noexcept;

            inline iterator& operator--() noexcept;
            inline iterator operator--(int) noexcept;

            // Random-access
            inline iterator& operator+=(difference_type n) noexcept;
            inline iterator& operator-=(difference_type n) noexcept;
            inline iterator operator+(difference_type n) const noexcept;
            inline iterator operator-(difference_type n) const noexcept;

            inline difference_type operator-(const iterator& iter)
            const noexcept;

            inline bool operator==(const iterator& iter) const noexcept;
            inline bool operator!=(const iterator& iter) const noexcept;

            inline bool operator<(const iterator& iter) const noexcept;
            inline bool operator>(const iterator& iter) const noexcept;

            inline bool operator<=(const iterator& iter) const noexcept;
            inline bool operator>=(const iterator& iter) const noexcept;
        };


        //- Random-access iterator with const access
        class const_iterator
        {
            // Pointer to parent
            const T* const* ptr_;

        public:

            using iterator_category = std::random_access_iterator_tag;
            using value_type = const T;
            using difference_type = label;
            using pointer = const T*;
            using reference = const T&;

            //- Construct for a given entry
            inline const_iterator(const T* const* ptr) noexcept;

            //- Copy construct from non-const iterator
            inline const_iterator(const iterator& iter) noexcept;

            // Member Functions

            //- Return pointer, can be nullptr.
            inline pointer get() const;

            // Member Operators

            inline pointer operator->() const;
            inline reference operator*() const;

            reference operator()() const { return this->operator*(); }
            inline reference operator[](difference_type n) const;

            // Forward iteration
            inline const_iterator& operator++() noexcept;
            inline const_iterator operator++(int) noexcept;

            inline const_iterator& operator--() noexcept;
            inline const_iterator operator--(int) noexcept;

            // Random-access
            inline const_iterator& operator+=(difference_type n) noexcept;
            inline const_iterator& operator-=(difference_type n) noexcept;
            inline const_iterator operator+(difference_type n) const noexcept;
            inline const_iterator operator-(difference_type n) const noexcept;

            inline difference_type operator-(const const_iterator& iter)
            const noexcept;

            inline bool operator==(const const_iterator& iter) const noexcept;
            inline bool operator!=(const const_iterator& iter) const noexcept;

            inline bool operator<(const const_iterator& iter) const noexcept;
            inline bool operator>(const const_iterator& iter) const noexcept;

            inline bool operator<=(const const_iterator& iter) const noexcept;
            inline bool operator>=(const const_iterator& iter) const noexcept;
        };


        //- Return an iterator to begin traversing the UPtrList
        inline iterator begin() noexcept;

        //- Return an iterator to end traversing the UPtrList
        inline iterator end() noexcept;

        //- Return an const_iterator to begin traversing the UPtrList
        inline const_iterator cbegin() const noexcept;

        //- Return an const_iterator to end traversing the UPtrList
        inline const_iterator cend() const noexcept;

        //- Return an const_iterator to begin traversing the UPtrList
        inline const_iterator begin() const noexcept;

        //- Return an const_iterator to end traversing the UPtrList
        inline const_iterator end() const noexcept;


    // IOstream operator

        //- Write UPtrList to Ostream
        friend Ostream& operator<< <T>
        (
            Ostream& os,
            const UPtrList<T>& list
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UPtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
