/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::DynamicList

Description
    A 1D vector of objects of type \<T\> that resizes itself as necessary to
    accept the new objects.

    Internal storage is a compact array and the list can be shrunk to compact
    storage. The increase of list size is controlled by three template
    parameters, which allows the list storage to either increase by the given
    increment or by the given multiplier and divider (allowing non-integer
    multiples).

SourceFiles
    DynamicListI.H
    DynamicList.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicList_H
#define DynamicList_H

#include "List.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
class DynamicList;

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Ostream& operator<<
(
    Ostream& os,
    const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
);
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Istream& operator>>
(
    Istream& is,
    DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
);


/*---------------------------------------------------------------------------*\
                           Class DynamicList Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned SizeInc=0, unsigned SizeMult=2, unsigned SizeDiv=1>
class DynamicList
:
    public List<T>
{
    static_assert
    (
        (SizeInc || SizeMult) && SizeDiv,
        "Invalid sizing parameters"
    );

    // Private data

        //- The capacity (allocated size) of the underlying list.
        label capacity_;


public:

    // Related types

        //- Declare friendship with the List class
        friend class List<T>;


    // Constructors

        //- Construct null
        inline DynamicList();

        //- Construct given size.
        explicit inline DynamicList(const label nElem);

        //- Construct with given size and value for all elements.
        inline DynamicList(const label nElem, const T& a);

        //- Construct copy.
        inline DynamicList
        (
            const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
        );

        //- Construct from UList. Size set to UList size.
        //  Also constructs from DynamicList with different sizing parameters.
        explicit inline DynamicList(const UList<T>& lst);

        //- Construct given begin/end iterators.
        //  Uses std::distance to determine the size.
        template<class InputIterator>
        inline DynamicList(InputIterator begIter, InputIterator endIter);

        //- Construct from an initializer list. Size set to list size.
        explicit inline DynamicList(std::initializer_list<T> lst);

        //- Construct from UIndirectList. Size set to UIndirectList size.
        explicit inline DynamicList(const UIndirectList<T>& lst);

        //- Construct by transferring the parameter contents
        explicit inline DynamicList(const Xfer<List<T>>& lst);

        //- Construct from Istream. Size set to size of list read.
        explicit DynamicList(Istream& is);


    // Member Functions

        // Access

            //- Size of the underlying storage.
            inline label capacity() const;

        // Edit

            //- Alter the size of the underlying storage.
            //  The addressed size will be truncated if needed to fit, but will
            //  remain otherwise untouched.
            //  Use this or reserve() in combination with append().
            inline void setCapacity(const label nElem);

            //- Alter the addressed list size.
            //  New space will be allocated if required.
            //  Use this to resize the list prior to using the operator[] for
            //  setting values (as per List usage).
            inline void setSize(const label nElem);

            //- Alter the addressed list size and fill new space with a
            //  constant.
            inline void setSize(const label nElem, const T& t);

            //- Alter the addressed list size.
            //  New space will be allocated if required.
            //  Use this to resize the list prior to using the operator[] for
            //  setting values (as per List usage).
            inline void resize(const label nElem);

            //- Alter the addressed list size and fill new space with a
            //  constant.
            inline void resize(const label nElem, const T& t);

            //- Reserve allocation space for at least this size.
            //  Never shrinks the allocated size, use setCapacity() for that.
            inline void reserve(const label nElem);

            //- Clear the addressed list, i.e. set the size to zero.
            //  Allocated size does not change
            inline void clear();

            //- Clear the list and delete storage.
            inline void clearStorage();

            //- Shrink the allocated space to the number of elements used.
            //  Returns a reference to the DynamicList.
            inline DynamicList<T, SizeInc, SizeMult, SizeDiv>& shrink();

            //- Transfer contents of the argument List into this.
            inline void transfer(List<T>& lst);

            //- Transfer contents of the argument DynamicList into this.
            inline void transfer
            (
                DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
            );

            //- Transfer contents to the Xfer container as a plain List
            inline Xfer<List<T>> xfer();


        // Member Operators

            //- Append an element at the end of the list
            inline DynamicList<T, SizeInc, SizeMult, SizeDiv>& append
            (
                const T& t
            );

            //- Append a List at the end of this list
            inline DynamicList<T, SizeInc, SizeMult, SizeDiv>& append
            (
                const UList<T>& lst
            );

            //- Append an initializer list at the end of this list.
            inline DynamicList<T, SizeInc, SizeMult, SizeDiv>& append
            (
                std::initializer_list<T> lst
            );

            //- Append a UIndirectList at the end of this list
            inline DynamicList<T, SizeInc, SizeMult, SizeDiv>& append
            (
                const UIndirectList<T>& lst
            );

            //- Remove and return the top element
            inline T remove();

            //- Return non-const access to an element, resizing list if
            //  necessary
            inline T& operator()(const label elemI);

            //- Assignment of all addressed entries to the given value
            inline void operator=(const T& t);

            //- Assignment to DynamicList
            inline void operator=
            (
                const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
            );

            //- Assignment to UList
            inline void operator=(const UList<T>& lst);

            //- Assignment from initializer list
            inline void operator=(std::initializer_list<T> lst);

            //- Assignment to UIndirectList
            inline void operator=(const UIndirectList<T>& lst);


        // STL member functions

            //- Erase an element, move the remaining elements to fill the gap
            //  and resize the List
            typename UList<T>::iterator erase
            (
                typename UList<T>::iterator curIter
            );


        // IOstream operators

            // Write DynamicList to Ostream.
            friend Ostream& operator<< <T, SizeInc, SizeMult, SizeDiv>
            (
                Ostream& os,
                const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
            );

            //- Read from Istream, discarding contents of existing DynamicList.
            friend Istream& operator>> <T, SizeInc, SizeMult, SizeDiv>
            (
                Istream& is,
                DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DynamicList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
