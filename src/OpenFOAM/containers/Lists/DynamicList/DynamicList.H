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
template<class T, int SizeMin> class DynamicList;

template<class T, int SizeMin>
Ostream& operator<<
(
    Ostream& os,
    const DynamicList<T, SizeMin>& lst
);

template<class T, int SizeMin>
Istream& operator>>
(
    Istream& is,
    DynamicList<T, SizeMin>& lst
);


/*---------------------------------------------------------------------------*\
                           Class DynamicList Declaration
\*---------------------------------------------------------------------------*/

template<class T, int SizeMin = 16>
class DynamicList
:
    public List<T>
{
    static_assert(SizeMin > 0, "Invalid min size parameter");

    // Private data

        //- The capacity (allocated size) of the underlying list.
        label capacity_;

private:

    // Private Member Functions

        //- Remove elements in range
        label removeElements(const labelRange& slice);

        //- Subset elements in range
        label subsetElements(const labelRange& slice);


protected:

    // Protected Member Functions

        //- Copy assignment from another list
        template<class ListType>
        inline void assignDynList(const ListType& lst);

public:

    // Related types

        //- Declare friendship with the List class
        friend class List<T>;


    // Constructors

        //- Construct null
        inline DynamicList();

        //- Construct an empty list with given reserve size.
        explicit inline DynamicList(const label nElem);

        //- Construct with given size and value for all elements.
        inline DynamicList(const label nElem, const T& val);

        //- Construct with given size initializing all elements to zero
        inline DynamicList(const label nElem, const zero);

        //- Copy construct.
        inline DynamicList(const DynamicList<T, SizeMin>& lst);

        //- Copy construct from DynamicList with different sizing parameters
        template<int AnySizeMin>
        inline DynamicList(const DynamicList<T, AnySizeMin>& lst);

        //- Construct from UList. Size set to UList size.
        //  Also constructs from DynamicList with different sizing parameters.
        explicit inline DynamicList(const UList<T>& lst);

        //- Construct from a FixedList
        template<unsigned FixedSize>
        inline DynamicList(const FixedList<T, FixedSize>& lst);

        //- Construct given begin/end iterators.
        //  Uses std::distance to determine the size.
        template<class InputIterator>
        inline DynamicList(InputIterator begIter, InputIterator endIter);

        //- Construct from an initializer list. Size set to list size.
        explicit inline DynamicList(std::initializer_list<T> lst);

        //- Construct from UIndirectList. Size set to UIndirectList size.
        explicit inline DynamicList(const UIndirectList<T>& lst);

        //- Transfer (move) construct
        explicit inline DynamicList(const Xfer<List<T>>& lst);

        //- Move construct.
        inline DynamicList(DynamicList<T, SizeMin>&& lst);

        //- Move construct with different sizing parameters
        template<int AnySizeMin>
        inline DynamicList(DynamicList<T, AnySizeMin>&& lst);

        //- Move construct from List
        inline DynamicList(List<T>&& lst);

        //- Move construct from SortableList
        DynamicList(SortableList<T>&& lst);

        //- Construct from Istream. Size set to size of list read.
        explicit DynamicList(Istream& is);


    // Member Functions

      // Access

        //- Normal lower capacity limit - the SizeMin template parameter
        inline label min_size() const;

        //- Size of the underlying storage.
        inline label capacity() const;


        // Edit

            //- Alter the size of the underlying storage.
            //  The addressed size will be truncated if needed to fit, but will
            //  remain otherwise untouched.
            //  Use this or reserve() in combination with append().
            inline void setCapacity(const label nElem);

            //- Alter addressable list size.
            //  New space will be allocated if required.
            //  Use this to resize the list prior to using the operator[] for
            //  setting values (as per List usage).
            inline void setSize(const label nElem);

            //- Alter addressable list size and fill new space with constant.
            inline void setSize(const label nElem, const T& val);

            //- Alter addressable list size.
            //  New space will be allocated if required.
            //  Use this to resize the list prior to using the operator[] for
            //  setting values (as per List usage).
            inline void resize(const label nElem);

            //- Alter addressable list size and fill new space with constant.
            inline void resize(const label nElem, const T& val);

            //- Reserve allocation space for at least this size.
            //  Never shrinks the allocated size, use setCapacity() for that.
            inline void reserve(const label nElem);

            //- Clear the addressed list, i.e. set the size to zero.
            //  Allocated size does not change
            inline void clear();

            //- Clear the list and delete storage.
            inline void clearStorage();

            //- Expand the addressable size to fit the allocated capacity.
            //  Returns the previous addressable size.
            inline label expandStorage();

            //- Shrink the allocated space to the number of elements used.
            //  Returns a reference to the DynamicList.
            inline DynamicList<T, SizeMin>& shrink();

            //- Swap content with any sized DynamicList
            template<int AnySizeMin>
            inline void swap(DynamicList<T, AnySizeMin>& lst);

            //- Transfer contents of the argument List into this.
            inline void transfer(List<T>& lst);

            //- Transfer contents of any sized DynamicList into this.
            template<int AnySizeMin>
            inline void transfer(DynamicList<T, AnySizeMin>& lst);

            //- Transfer contents  of the argument SortableList into this.
            inline void transfer(SortableList<T>& lst);

            //- Transfer contents to the Xfer container as a plain List
            inline Xfer<List<T>> xfer();


            //- Append an element to the end of this list.
            inline DynamicList<T, SizeMin>& append(const T& val);

            //- Move append an element
            inline DynamicList<T, SizeMin>& append(T&& val);

            //- Append another list to the end of this list.
            inline DynamicList<T, SizeMin>& append(const UList<T>& lst);

            //- Append a FixedList to the end of this list.
            template<unsigned FixedSize>
            inline DynamicList<T, SizeMin>&
            append(const FixedList<T, FixedSize>& lst);

            //- Append an initializer list at the end of this list.
            inline DynamicList<T, SizeMin>&
            append(std::initializer_list<T> lst);

            //- Append a UIndirectList at the end of this list
            inline DynamicList<T, SizeMin>&
            append(const UIndirectList<T>& lst);

            //- Move append list
            inline DynamicList<T, SizeMin>& append(List<T>&& lst);

            //- Move append list
            inline DynamicList<T, SizeMin>&
            append(DynamicList<T, SizeMin>&& lst);

            //- Move append list
            template<int AnySizeMin>
            inline DynamicList<T, SizeMin>&
            append(DynamicList<T, AnySizeMin>&& lst);

            //- Move append list
            inline DynamicList<T, SizeMin>&
            append(SortableList<T>&& lst);

            //- Remove and return the last element. Fatal on an empty list.
            inline T remove();

            //- Remove and return the specified element. Fatal on an empty list.
            //  With fast=true (operates in constant time), the place of the
            //  removed element is swapped with the last one in the list, which
            //  changes the ordering.
            //  With fast=false (operates in linear time), the elements
            //  are swapped down in the list to preserve ordering.
            inline T remove(const label idx, const bool fast=false);

            //- Remove a (start,size) subset from the list.
            //  The range is subsetted with the list size itself to ensure
            //  result always addresses a valid section of the list.
            //  Remaining elements are moved down.
            inline label remove(const labelRange& range);

            //- Remove a (start,size) subset from the list.
            inline label remove(std::initializer_list<label> start_size);

            //- Retain a (start,size) subset from the list.
            //  The range is subsetted with the list size itself to ensure
            //  result always addresses a valid section of the list.
            //  Remaining elements are moved down.
            inline label subset(const labelRange& range);

            //- Retain a (start,size) subset from List.
            inline label subset(std::initializer_list<label> start_size);


        // Member Operators

            //- Return non-const access to an element, resizing list if
            //  necessary
            inline T& operator()(const label i);

            //- Assignment of all addressed entries to the given value
            inline void operator=(const T& val);

            //- Assignment of all entries to zero
            inline void operator=(const zero);

            //- Assignment to UList
            inline void operator=(const UList<T>& lst);

            //- Assignment to FixedList
            template<unsigned FixedSize>
            inline void operator=(const FixedList<T, FixedSize>& lst);

            //- Assignment to DynamicList
            inline void operator=(const DynamicList<T, SizeMin>& lst);

            //- Assignment from DynamicList with different sizing parameters
            template<int AnySizeMin>
            inline void operator=(const DynamicList<T, AnySizeMin>& lst);

            //- Assignment from initializer list
            inline void operator=(std::initializer_list<T> lst);

            //- Assignment to UIndirectList
            inline void operator=(const UIndirectList<T>& lst);

            //- Move assignment
            inline void operator=(List<T>&& lst);

            //- Move assignment
            inline void operator=(DynamicList<T, SizeMin>&& lst);

            //- Move assignment
            template<int AnySizeMin>
            inline void operator=(DynamicList<T, AnySizeMin>&& lst);

            //- Move assignment
            inline void operator=(SortableList<T>&& lst);


        // IOstream operators

            // Write DynamicList to Ostream.
            friend Ostream& operator<< <T, SizeMin>
            (
                Ostream& os,
                const DynamicList<T, SizeMin>& lst
            );

            //- Read from Istream, discarding contents of existing DynamicList.
            friend Istream& operator>> <T, SizeMin>
            (
                Istream& is,
                DynamicList<T, SizeMin>& lst
            );
};


// Global Functions

// Exchange contents of lists - see DynamicList::swap().
template<class T, int SizeMin1, int SizeMin2>
inline void Swap(DynamicList<T, SizeMin1>& a, DynamicList<T, SizeMin2>& b);


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
