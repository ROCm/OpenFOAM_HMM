/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    storage. The increase of list size uses a doubling strategy, with the
    SizeMin template parameter dictating a lower bound.

SourceFiles
    DynamicList.C
    DynamicListI.H
    DynamicListIO.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicList_H
#define DynamicList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T, int SizeMin> class DynamicList;

template<class T, int SizeMin>
Istream& operator>>(Istream& is, DynamicList<T, SizeMin>& list);

template<class T, int SizeMin>
Ostream& operator<<(Ostream& os, const DynamicList<T, SizeMin>& list);


/*---------------------------------------------------------------------------*\
                           Class DynamicList Declaration
\*---------------------------------------------------------------------------*/

template<class T, int SizeMin = 16>
class DynamicList
:
    public List<T>
{
    static_assert(SizeMin > 0, "Invalid min size parameter");

    // Private Data

        //- The capacity (allocated size) of the underlying list.
        label capacity_;


    // Private Member Functions

        //- Remove elements in range
        label removeElements(const labelRange& slice);

        //- Subset elements in range
        label subsetElements(const labelRange& slice);

        //- Copy assignment from another list
        template<class ListType>
        inline void assignDynList(const ListType& list);

        //- Alter the size of the underlying storage
        //  The 'nocopy' option will not attempt to recover old content
        inline void doCapacity(const bool nocopy, const label len);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size, use setCapacity() for that.
        //  The 'nocopy' option will not attempt to recover old content
        inline void doReserve(const bool nocopy, const label len);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size, use setCapacity() for that.
        //  The 'nocopy' option will not attempt to recover old content
        inline void doResize(const bool nocopy, const label len);


public:

    // Constructors

        //- Default construct, an empty list without allocation.
        inline constexpr DynamicList() noexcept;

        //- Construct an empty list with given reserve size.
        inline explicit DynamicList(const label len);

        //- Construct with given size and value for all elements.
        inline DynamicList(const label len, const T& val);

        //- Construct with given size initializing all elements to zero
        inline DynamicList(const label len, const Foam::zero);

        //- Copy construct.
        inline DynamicList(const DynamicList<T, SizeMin>& lst);

        //- Copy construct from DynamicList with different sizing parameters
        template<int AnySizeMin>
        inline DynamicList(const DynamicList<T, AnySizeMin>& lst);

        //- Construct from UList. Size set to UList size.
        //  Also constructs from DynamicList with different sizing parameters.
        inline explicit DynamicList(const UList<T>& lst);

        //- Construct from a FixedList
        template<unsigned N>
        inline explicit DynamicList(const FixedList<T, N>& lst);

        //- Construct from an initializer list. Size set to list size.
        inline explicit DynamicList(std::initializer_list<T> lst);

        //- Construct from IndirectList. Size set to addressing size.
        template<class Addr>
        inline explicit DynamicList(const IndirectListBase<T, Addr>& lst);

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
        static constexpr label min_size() noexcept { return SizeMin; }

        //- Size of the underlying storage.
        inline label capacity() const noexcept;

        //- Number of contiguous bytes of the underlying storage.
        //  \note Only meaningful for contiguous data
        inline std::streamsize capacity_bytes() const noexcept;


    // Sizing

        //- Alter the size of the underlying storage.
        //  The addressed size will be truncated if needed to fit, but will
        //  remain otherwise untouched.
        //  Use this or reserve() in combination with append().
        inline void setCapacity(const label len);

        //- Alter the size of the underlying storage,
        //- \em without retaining old content.
        //  The addressed size will be truncated if needed to fit, but will
        //  remain otherwise untouched.
        inline void setCapacity_nocopy(const label len);

        //- Change the value for the list capacity directly (ADVANCED, UNSAFE)
        //- Does not perform any memory management or resizing.
        inline void setCapacity_unsafe(const label len) noexcept;

        //- Reserve allocation space for at least this size, allocating new
        //- space if required and \em retaining old content.
        //  Never shrinks the allocated size, use setCapacity() for that.
        inline void reserve(const label len);

        //- Reserve allocation space for at least this size, allocating new
        //- space if required \em without retaining old content.
        //  Never shrinks the allocated size, use setCapacity() for that.
        inline void reserve_nocopy(const label len);

        //- Alter addressable list size, allocating new space if required
        //- while \em recovering old content.
        //  If no reallocation is required, the contents remain untouched.
        //  Otherwise new entries will be uninitialized.
        //  Use this to resize the list prior to using the operator[] for
        //  setting values (as per List usage).
        inline void resize(const label len);

        //- Alter addressable size and fill new entries with constant value
        inline void resize(const label len, const T& val);

        //- Alter addressable list size, allocating new space if required
        //- \em without necessarily recovering old content.
        //  If no reallocation is required, the contents remain untouched.
        //  Otherwise all entries will be uninitialized.
        inline void resize_nocopy(const label len);

        //- Same as resize()
        void setSize(const label n) { this->resize(n); }

        //- Same as resize()
        void setSize(const label n, const T& val) { this->resize(n, val); }

        //- Clear the addressed list, i.e. set the size to zero.
        //  Allocated size does not change
        inline void clear() noexcept;

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Expand the addressable size to fit the allocated capacity.
        //  Returns the previous addressable size.
        inline label expandStorage() noexcept;

        //- Shrink the allocated space to the number of elements used.
        inline void shrinkStorage();

        //- Shrink the allocated space to the number of elements used.
        //  Returns a reference to the DynamicList.
        inline DynamicList<T, SizeMin>& shrink();


    // Edit

            //- Swap content, independent of sizing parameter
            template<int AnySizeMin>
            inline void swap(DynamicList<T, AnySizeMin>& other);

            //- Transfer contents of the argument List into this.
            inline void transfer(List<T>& list);

            //- Transfer contents of any sized DynamicList into this.
            template<int AnySizeMin>
            inline void transfer(DynamicList<T, AnySizeMin>& list);

            //- Transfer contents  of the argument SortableList into this.
            inline void transfer(SortableList<T>& list);

            //- Append an element to the end of this list.
            inline DynamicList<T, SizeMin>& append(const T& val);

            //- Move append an element
            inline DynamicList<T, SizeMin>& append(T&& val);

            //- Append another list to the end of this list.
            inline DynamicList<T, SizeMin>& append(const UList<T>& lst);

            //- Append a FixedList to the end of this list.
            template<unsigned N>
            inline DynamicList<T, SizeMin>&
            append(const FixedList<T, N>& lst);

            //- Append an initializer list at the end of this list.
            inline DynamicList<T, SizeMin>&
            append(std::initializer_list<T> lst);

            //- Append a IndirectList at the end of this list
            template<class Addr>
            inline DynamicList<T, SizeMin>&
            append(const IndirectListBase<T, Addr>& lst);

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

            //- Append an element if not already in the list.
            //  \return the change in list length
            inline label appendUniq(const T& val);

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

            //- Return non-const access to an element, resizing list if needed
            inline T& operator()(const label i);

            //- Assignment of all addressed entries to the given value
            inline void operator=(const T& val);

            //- Assignment of all entries to zero
            inline void operator=(const Foam::zero);

            //- Assignment to UList
            inline void operator=(const UList<T>& lst);

            //- Assignment to FixedList
            template<unsigned N>
            inline void operator=(const FixedList<T, N>& lst);

            //- Assignment to DynamicList
            inline void operator=(const DynamicList<T, SizeMin>& lst);

            //- Assignment from DynamicList with different sizing parameters
            template<int AnySizeMin>
            inline void operator=(const DynamicList<T, AnySizeMin>& lst);

            //- Assignment from initializer list
            inline void operator=(std::initializer_list<T> lst);

            //- Assignment from IndirectList
            template<class Addr>
            inline void operator=(const IndirectListBase<T, Addr>& lst);

            //- Move assignment
            inline void operator=(List<T>&& lst);

            //- Move assignment
            inline void operator=(DynamicList<T, SizeMin>&& lst);

            //- Move assignment
            template<int AnySizeMin>
            inline void operator=(DynamicList<T, AnySizeMin>&& lst);

            //- Move assignment
            inline void operator=(SortableList<T>&& lst);


    // Reading/writing

        //- Read from Istream, discarding existing contents
        Istream& readList(Istream& is);


    // IOstream Operators

        //- Use the readList() method to read contents from Istream.
        friend Istream& operator>> <T, SizeMin>
        (
            Istream& is,
            DynamicList<T, SizeMin>& list
        );

        //- Write to Ostream
        friend Ostream& operator<< <T, SizeMin>
        (
            Ostream& os,
            const DynamicList<T, SizeMin>& list
        );
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read List contents from Istream
template<class T, int SizeMin>
Istream& operator>>(Istream& is, DynamicList<T, SizeMin>& list)
{
    return list.readList(is);
}


//- Write List to Ostream, as per UList::writeList() with default length.
template<class T, int SizeMin>
Ostream& operator<<(Ostream& os, const DynamicList<T, SizeMin>& list)
{
    return (os << static_cast<const UList<T>&>(list));
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Exchange contents of lists - see DynamicList::swap().
template<class T, int SizeMinA, int SizeMinB>
inline void Swap(DynamicList<T, SizeMinA>& a, DynamicList<T, SizeMinB>& b)
{
    a.swap(b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hashing for List data
template<class T, int SizeMin>
struct Hash<DynamicList<T, SizeMin>> : List<T>::hasher {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DynamicList.C"
    #include "DynamicListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
