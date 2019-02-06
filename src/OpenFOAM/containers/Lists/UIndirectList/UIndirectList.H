/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::UIndirectList

Description
    A List with indirect addressing.

    Like IndirectList but does not store addressing.

    Note the const_cast of the list values. This is so we can use it both
    on const and non-const lists. Alternative would be to have a const_
    variant etc.

SourceFiles
    UIndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef UIndirectList_H
#define UIndirectList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class T> class UIndirectList;

// Commonly required list types
typedef UIndirectList<bool> boolUIndList;
typedef UIndirectList<label> labelUIndList;

/*---------------------------------------------------------------------------*\
                        Class UIndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UIndirectList
{
    // Private data

        UList<T>& values_;
        const labelUList& addressing_;

public:

    // STL type definitions

        //- Type of values the list contains.
        typedef T value_type;

        //- The pointer type for non-const access to value_type items
        typedef T* pointer;

        //- The pointer type for const access to value_type items
        typedef const T* const_pointer;

        //- The type used for storing into value_type objects
        typedef T& reference;

        //- The type used for reading from constant value_type objects.
        typedef const T& const_reference;

        //- The type to represent the size of a UList
        typedef label size_type;

        //- The difference between iterator objects
        typedef label difference_type;

        //- Forward iterator with non-const access
        class iterator;

        //- Forward iterator with const access
        class const_iterator;


    // Constructors

        //- Copy construct from the values list and the addressing array
        inline UIndirectList(const UList<T>& values, const labelUList& addr);


    // Member Functions

    // Access

        //- Return the number of elements in the list
        inline label size() const;

        //- Return true if the list is empty (ie, size() is zero).
        inline bool empty() const;

        //- True if all entries have identical values, and list is non-empty
        inline bool uniform() const;

        //- Return the first element of the list.
        inline T& first();

        //- Return first element of the list.
        inline const T& first() const;

        //- Return the last element of the list.
        inline T& last();

        //- Return the last element of the list.
        inline const T& last() const;

        //- Return the complete list of values
        inline const UList<T>& completeList() const;

        //- Return the list addressing
        inline const List<label>& addressing() const;


    // Search

        //- Find index of the first occurrence of the value.
        //  When start is specified, any occurences before start are ignored.
        //  Linear search.
        //  \return -1 if not found.
        label find(const T& val, const label start=0) const;

        //- Find index of the last occurrence of the value.
        //  When pos is specified, any occurrences after pos are ignored.
        //  Linear search.
        //  \return -1 if not found.
        label rfind(const T& val, const label pos=-1) const;

        //- True if the value if found in the list. Linear search.
        inline bool found(const T& val, const label start=0) const;


    // Member Operators

        //- Return the addressed elements as a List
        inline List<T> operator()() const;

        //- Return non-const access to an element
        inline T& operator[](const label i);

        //- Return const access to an element
        inline const T& operator[](const label i) const;

        //- Copy assignment from a list of the addressed elements
        inline void operator=(const UList<T>& rhs);

        //- Copy assignment from a indirect list of the addressed elements
        inline void operator=(const UIndirectList<T>& rhs);

        //- Assignment of all entries to the given value
        inline void operator=(const T& val);

        //- An iterator for an indirect list
        class iterator
        {
            typename UList<T>::pointer data_;
            labelUList::const_iterator base_;

        public:

            using difference_type = label;
            using value_type = T;
            using pointer = T*;
            using reference = T&;
            using iterator_category = std::forward_iterator_tag;

            iterator
            (
                UList<T>& list,
                labelUList::const_iterator baseIter
            )
            :
                data_(list.begin()),
                base_(baseIter)
            {}

            reference operator*() const
            {
                return data_[*base_];
            }

            iterator& operator++()
            {
                ++base_;
                return *this;
            }

            bool operator==(iterator& rhs) const
            {
                return base_ == rhs.base_;
            }

            bool operator!=(iterator& rhs) const
            {
                return (base_ != rhs.base_);
            }
        };


        //- A const iterator for an indirect list
        class const_iterator
        {
            typename UList<T>::const_pointer data_;
            labelUList::const_iterator base_;

        public:

            using difference_type = label;
            using value_type = const T;
            using pointer = const T*;
            using reference = const T&;
            using iterator_category = std::forward_iterator_tag;

            const_iterator
            (
                const UList<T>& list,
                labelUList::const_iterator baseIter
            )
            :
                data_(list.begin()),
                base_(baseIter)
            {}

            reference operator*() const
            {
                return data_[*base_];
            }

            const_iterator& operator++()
            {
                ++base_;
                return *this;
            }

            bool operator==(const_iterator& rhs) const
            {
                return base_ == rhs.base_;
            }

            bool operator!=(const_iterator& rhs) const
            {
                return base_ != rhs.base_;
            }
        };


    // iterator (non-const)

        //- Return an iterator at begin of list
        inline iterator begin()
        {
            return iterator(values_, addressing_.cbegin());
        }

        //- Return an iterator at end of list
        inline iterator end()
        {
            return iterator(values_, addressing_.cend());
        }


    // iterator (const)

        //- Return a const_iterator at begin of list
        inline const_iterator cbegin() const
        {
            return const_iterator(values_, addressing_.cbegin());
        }

        //- Return a const_iterator at end of list
        inline const_iterator cend() const
        {
            return const_iterator(values_, addressing_.cend());
        }

        //- Return a const_iterator at end of list
        inline const_iterator begin() const
        {
            return cbegin();
        }

        //- Return a const_iterator at end of list
        inline const_iterator end() const
        {
            return cend();
        }


    // Writing

        //- Write List, with line-breaks in ASCII when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const;
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Write List to Ostream, as per UList::writeList() with default length.
//  The default short-length is given by Detail::ListPolicy::short_length
template<class T>
Ostream& operator<<(Ostream& os, const UIndirectList<T>& list)
{
    return list.writeList(os, Detail::ListPolicy::short_length<T>::value);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UIndirectListI.H"

#ifdef NoRepository
    #include "UIndirectList.C"
    #include "UIndirectListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
