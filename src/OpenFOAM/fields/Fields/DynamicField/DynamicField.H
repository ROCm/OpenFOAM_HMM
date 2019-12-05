/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::DynamicField

Description
    Dynamically sized Field.

SourceFiles
    DynamicFieldI.H
    DynamicField.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicField_H
#define DynamicField_H

#include "Field.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

template<class T, int SizeMin> class DynamicField;

template<class T, int SizeMin>
Ostream& operator<<
(
    Ostream& os,
    const DynamicField<T, SizeMin>& fld
);

template<class T, int SizeMin>
Istream& operator>>
(
    Istream& is,
    DynamicField<T, SizeMin>& fld
);


/*---------------------------------------------------------------------------*\
                        Class DynamicField Declaration
\*---------------------------------------------------------------------------*/

template<class T, int SizeMin=64>
class DynamicField
:
    public Field<T>
{
    static_assert(SizeMin > 0, "Invalid min size parameter");

    // Private data

        //- The capacity (allocated size) of the underlying field.
        label capacity_;


    // Private Member Functions

        //- Copy assignment from another list
        template<class ListType>
        inline void assignDynField(const ListType& list);

public:

    // Static Member Functions

        //- Return a null field
        inline static const DynamicField<T, SizeMin>& null()
        {
            return NullObjectRef<DynamicField<T, SizeMin>>();
        }


    // Constructors

        //- Construct null
        inline constexpr DynamicField() noexcept;

        //- Construct empty field with given reserve size.
        inline explicit DynamicField(const label len);

        //- Construct given size and initial value
        inline DynamicField(const label len, const T& val);

        //- Construct given size and initial value of zero
        inline DynamicField(const label len, const zero);

        //- Copy construct
        inline DynamicField(const DynamicField<T, SizeMin>& list);

        //- Copy construct with different sizing parameters
        template<int AnySizeMin>
        inline DynamicField(const DynamicField<T, AnySizeMin>& list);

        //- Copy construct from UList. Size set to UList size.
        //  Also constructs from DynamicField with different sizing parameters.
        inline explicit DynamicField(const UList<T>& list);

        //- Copy construct from IndirectList
        template<class Addr>
        inline explicit DynamicField(const IndirectListBase<T, Addr>& list);

        //- Move construct from List contents
        inline explicit DynamicField(List<T>&& content);

        //- Move construct from dynamic Field contents
        inline DynamicField(DynamicField<T, SizeMin>&& content);

        //- Move construct with different sizing parameters
        template<int AnySizeMin>
        inline DynamicField(DynamicField<T, AnySizeMin>&& content);

        //- Construct by 1 to 1 mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelUList& mapAddressing
        );

        //- Construct by interpolative mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Construct by mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const FieldMapper& map
        );

        //- Construct from Istream. Size set to size of list read.
        explicit DynamicField(Istream& is);

        //- Clone
        tmp<DynamicField<T, SizeMin>> clone() const;


    // Member Functions

    // Access

        //- Size of the underlying storage.
        inline label capacity() const noexcept;


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

        //- Alter the addressed list size and fill new space with a constant.
        inline void setSize(const label nElem, const T& val);

        //- Alter the addressed list size.
        //  New space will be allocated if required.
        //  Use this to resize the list prior to using the operator[] for
        //  setting values (as per List usage).
        inline void resize(const label nElem);

        //- Alter the addressed list size and fill new space with a
        //  constant.
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
        //  Returns a reference to the DynamicField.
        inline DynamicField<T, SizeMin>& shrink();

        //- Swap content with any sized DynamicField
        template<int AnySizeMin>
        inline void swap(DynamicField<T, AnySizeMin>& list);

        //- Transfer the parameter contents into this
        inline void transfer(List<T>& list);

        //- Transfer the parameter contents into this
        template<int AnySizeMin>
        inline void transfer(DynamicList<T, AnySizeMin>& list);

        //- Transfer the parameter contents into this
        template<int AnySizeMin>
        inline void transfer(DynamicField<T, AnySizeMin>& list);


        //- Append an element at the end of the list
        inline DynamicField<T, SizeMin>&
        append(const T& val);

        //- Append a List at the end of this list
        inline DynamicField<T, SizeMin>&
        append(const UList<T>& list);

        //- Remove and return the top element
        inline T remove();


    // Member Operators

        //- Return non-const access to an element, resizing list if necessary
        inline T& operator()(const label i);

        //- Assign addressed entries to the given value
        inline void operator=(const T& val);

        //- Copy assignment
        inline void operator=(const UList<T>& list);

        //- Copy assignment
        inline void operator=(const DynamicField<T, SizeMin>& list);

        //- Move assignment
        inline void operator=(List<T>&& list);

        //- Move assignment
        inline void operator=(DynamicField<T, SizeMin>&& list);

        //- Move assignment
        template<int AnySizeMin>
        inline void operator=(DynamicField<T, AnySizeMin>&& list);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DynamicField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
