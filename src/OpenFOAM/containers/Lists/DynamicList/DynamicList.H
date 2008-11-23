/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::DynamicList

Description
    A 1D vector of objects of type \<T\> that resizes itself as necessary to
    accept the new objects.

    Internal storage is a compact array and the list can be shrunk to compact
    storage. The increase of list size is controlled by three template
    parameters, which allows the list storage to either increase by the given
    increment or the given multiplier and divider (allowing non-integer
    multiples).

SourceFiles
    DynamicListI.H
    DynamicList.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicList_H
#define DynamicList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
class DynamicList;

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Ostream& operator<<
(
    Ostream&,
    const DynamicList<T, SizeInc, SizeMult, SizeDiv>&
);
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Istream& operator>>
(
    Istream&,
    DynamicList<T, SizeInc, SizeMult, SizeDiv>&
);


/*---------------------------------------------------------------------------*\
                           Class DynamicList Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned SizeInc=0, unsigned SizeMult=2, unsigned SizeDiv=1>
class DynamicList
:
    public List<T>
{
    // Private data

        //- Allocated size for underlying List.
        label allocSize_;

    // Private Member Functions

        // Disabled, since the usefulness and semantics are not quite clear
        void setSize(const label, const T&);

public:

    // Related types

        //- Declare friendship with the List class
        friend class List<T>;

    // Constructors

        //- Construct null
        inline DynamicList();

        //- Construct given size.
        explicit inline DynamicList(const label);

        //- Construct from UList. Size set to UList size.
        explicit inline DynamicList(const UList<T>&);

        //- Construct from Istream. Size set to size of read list.
        explicit DynamicList(Istream&);


    // Member Functions

        // Access

        //- Size of the underlying storage.
        inline label allocSize() const;

        // Edit

        //- Alter the size of the underlying storage.
        //  The addressed size will be truncated if needed to fit, but will
        //  otherwise remain untouched.
        inline void allocSize(const label);

        // CURRENT BEHAVIOUR
        //- When the new size is greater than the addressed list size, the
        //  allocated list sizes is adjusted and the
        //  addressed size does not change.
        //  Otherwise the addressed list size is just reduced and the
        //  allocated size does not change.
        //
        // PROPOSED BEHAVIOUR
        //- Alter the addressed list size.
        //  New space will be allocated if required.
        inline void setSize(const label);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size, use allocSize() for that.
        inline void reserve(const label);

        //- Clear the addressed list, i.e. set the size to zero.
        //  Allocated size does not change
        inline void clear();

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Shrink the allocated space to the number of elements used.
        //  Returns a reference to the DynamicList.
        inline DynamicList<T, SizeInc, SizeMult, SizeDiv>& shrink();

        //- Transfer contents of the argument List into this DynamicList
        inline void transfer(List<T>&);

        //- Transfer contents of the argument DynamicList into this DynamicList
        inline void transfer(DynamicList<T, SizeInc, SizeMult, SizeDiv>&);


    // Member Operators

        //- Append an element at the end of the list
        inline void append(const T& e);

        //- Remove and return the top element
        inline T remove();

        //- Return non-const access to an element, resizing list if necessary
        inline T& operator()(const label);

        //- Assignment of all addressed entries to the given value
        inline void operator=(const T&);

        //- Assignment from List<T>. Also handles assignment from DynamicList.
        inline void operator=(const List<T>&);


    // IOstream operators

        // Write DynamicList to Ostream.
        friend Ostream& operator<<
        #ifndef __CINT__
        <T, SizeInc, SizeMult, SizeDiv>
        #endif
        (
            Ostream&,
            const DynamicList<T, SizeInc, SizeMult, SizeDiv>&
        );

        //- Read from Istream, discarding contents of existing DynamicList.
        friend Istream& operator>>
        #ifndef __CINT__
        <T, SizeInc, SizeMult, SizeDiv>
        #endif
        (
            Istream&,
            DynamicList<T, SizeInc, SizeMult, SizeDiv>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynamicList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
