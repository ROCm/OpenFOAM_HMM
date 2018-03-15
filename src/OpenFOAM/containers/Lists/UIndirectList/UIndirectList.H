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
    Foam::UIndirectList

Description
    A List with indirect addressing.

    Like IndirectList but does not store addressing.

    Note the const_cast of the completeList. This is so we can use it both
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

// Forward declaration of friend functions and operators
template<class T> class UIndirectList;
template<class T> Ostream& operator<<(Ostream&, const UIndirectList<T>&);

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

        UList<T>& completeList_;
        const labelUList& addressing_;


public:

    // Constructors

        //- Construct given the complete list and the addressing array
        inline UIndirectList
        (
            const UList<T>& completeList,
            const labelUList& addr
        );


    // Member Functions

      // Access

        //- Return the number of elements in the list
        inline label size() const;

        //- Return true if the list is empty (ie, size() is zero).
        inline bool empty() const;

        //- Return the first element of the list.
        inline T& first();

        //- Return first element of the list.
        inline const T& first() const;

        //- Return the last element of the list.
        inline T& last();

        //- Return the last element of the list.
        inline const T& last() const;

        //- Return the complete list
        inline const UList<T>& completeList() const;

        //- Return the list addressing
        inline const List<label>& addressing() const;


      // Search

        //- Find index of the first occurence of the value.
        //  When start is specified, any occurences before start are ignored.
        //  Linear search.
        //  \return -1 if not found.
        label find(const T& val, const label start=0) const;

        //- Find index of the last occurence of the value.
        //  When pos is specified, any occurences after pos are ignored.
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

        //- Assignment to UList of addressed elements
        inline void operator=(const UList<T>& ae);

        //- Assignment to UIndirectList of addressed elements
        inline void operator=(const UIndirectList<T>& ae);

        //- Assignment of all entries to the given value
        inline void operator=(const T& t);


    // STL type definitions

        //- Type of values the list contains.
        typedef T value_type;

        //- The type used for storing into UList::value_type objects.
        typedef T& reference;

        //- The type used for reading from constant UList::value_type objects
        typedef const T& const_reference;

        //- The type that can represent the difference between any two
        //- UList iterator objects.
        typedef label difference_type;

        //- The type that can represent the size of a UList.
        typedef label size_type;


    // Writing

        //- Write the List, with line-breaks in ASCII if the list length
        //- exceeds shortListLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortListLen=0) const;


    // Ostream operator

        //- Write List to Ostream, as per writeList() with shortListLen=10
        friend Ostream& operator<< <T>
        (
            Ostream& os,
            const UIndirectList<T>& lst
        );
};


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
