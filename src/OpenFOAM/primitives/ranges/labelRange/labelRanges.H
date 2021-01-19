/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::labelRanges

Description
    A list of labelRange.

SourceFiles
    labelRanges.C

\*---------------------------------------------------------------------------*/

#ifndef labelRanges_H
#define labelRanges_H

#include "labelRange.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
class labelRanges;
Istream& operator>>(Istream& is, labelRanges& ranges);
Ostream& operator<<(Ostream& is, const labelRanges& ranges);

/*---------------------------------------------------------------------------*\
                         Class labelRanges Declaration
\*---------------------------------------------------------------------------*/

class labelRanges
:
    private Foam::DynamicList<Foam::labelRange>
{
    // Private Typedefs

        typedef Foam::DynamicList<Foam::labelRange> StorageContainer;


    // Private Member Functions

        //- Insert range before specified insertion index, by copying up
        void insertBefore(const label insert, const labelRange& range);

        //- Purge empty ranges, by copying down
        void purgeEmpty();

public:

    // STL type definitions

        //- Input iterator with const access
        class const_iterator;


    // Constructors

        //- Default construct
        labelRanges() = default;

        //- Construct given size
        inline explicit labelRanges(const label nElem);

        //- Construct from Istream.
        labelRanges(Istream& is);


    // Member Functions

        //- Clear the addressed list
        using DynamicList<labelRange>::clear;

        //- Return true if the list is empty
        using DynamicList<labelRange>::empty;

        //- Return true if the value is found any of the sub-ranges
        inline bool found(const label value) const;

        //- Add the range to the list
        bool add(const labelRange& range);

        //- Remove the range from the list
        bool remove(const labelRange& range);


    // STL iterator

        //- Forward input iterator with const access
        class const_iterator
        {
            //- The list for which this is an iterator
            const labelRanges* list_;

            //- The index into the list
            label index_;

            //- Index of current element at list-index
            label subIndex_;

        public:

        // Constructors

            //- Construct from range list at given index
            inline explicit constexpr const_iterator
            (
                const labelRanges* list,
                const label i = 0
            ) noexcept;


        // Member Operators

            //- Return the current label
            inline label operator*() const;

            inline const_iterator& operator++();
            inline const_iterator operator++(int);

            inline constexpr bool operator==(const const_iterator& iter)
            const noexcept;

            inline constexpr bool operator!=(const const_iterator& iter)
            const noexcept;
        };


        //- A const_iterator set to the beginning of the list
        inline const_iterator cbegin() const noexcept;

        //- A const_iterator set to beyond the end of the list
        inline const const_iterator cend() const noexcept;

        //- A const_iterator set to the beginning of the list
        inline const_iterator begin() const noexcept;

        //- A const_iterator set to beyond the end of the list
        inline const const_iterator end() const noexcept;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, labelRanges& ranges);
        friend Ostream& operator<<(Ostream& os, const labelRanges& ranges);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelRangesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
