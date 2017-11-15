/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

// Forward declaration of friend functions and operators

template<class T, int SizeMin> class DynamicField;

template<class T, int SizeMin>
Ostream& operator<<
(
    Ostream&,
    const DynamicField<T, SizeMin>&
);

template<class T, int SizeMin>
Istream& operator>>
(
    Istream&,
    DynamicField<T, SizeMin>&
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
        inline void assignDynField(const ListType& lst);

public:

    // Static Member Functions

        //- Return a null field
        inline static const DynamicField<T, SizeMin>& null()
        {
            return *reinterpret_cast<DynamicField<T, SizeMin>*>(0);
        }


    // Constructors

        //- Construct null
        inline DynamicField();

        //- Construct given size.
        explicit inline DynamicField(const label nElem);

        //- Construct from UList. Size set to UList size.
        //  Also constructs from DynamicField with different sizing parameters.
        explicit inline DynamicField(const UList<T>& lst);

        //- Construct by transferring the parameter contents
        explicit inline DynamicField(const Xfer<List<T>>& lst);

        //- Construct by transferring the parameter contents
        explicit inline DynamicField(const Xfer<Field<T>>& lst);

        //- Construct by 1 to 1 mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelList& mapAddressing
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

        //- Construct copy
        inline DynamicField(const DynamicField<T, SizeMin>& lst);

        //- Construct by transferring the Field contents
        inline DynamicField
        (
            const Xfer<DynamicField<T, SizeMin>>& lst
        );

        //- Construct from Istream. Size set to size of list read.
        explicit DynamicField(Istream& is);

        //- Clone
        tmp<DynamicField<T, SizeMin>> clone() const;


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

            //- Transfer contents to the Xfer container as a plain List
            inline Xfer<List<T>> xfer();


        // Member Operators

            //- Append an element at the end of the list
            inline DynamicField<T, SizeMin>&
            append(const T& val);

            //- Append a List at the end of this list
            inline DynamicField<T, SizeMin>&
            append(const UList<T>& lst);

            //- Remove and return the top element
            inline T remove();

            //- Return non-const access to an element, resizing list if
            //  necessary
            inline T& operator()(const label i);

            //- Assignment of all addressed entries to the given value
            inline void operator=(const T& val);

            //- Assignment to DynamicField
            inline void operator=
            (
                const DynamicField<T, SizeMin>& lst
            );

            //- Assignment to UList
            inline void operator=(const UList<T>& lst);
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
