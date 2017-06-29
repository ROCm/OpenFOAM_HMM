/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class labelRanges;
Istream& operator>>(Istream& is, labelRanges& ranges);
Ostream& operator<<(Ostream& is, const labelRanges& ranges);

/*---------------------------------------------------------------------------*\
                         Class labelRanges Declaration
\*---------------------------------------------------------------------------*/

class labelRanges
:
    private DynamicList<labelRange>
{
    // Private typedefs for convenience

        typedef DynamicList<labelRange> StorageContainer;


    // Private Member Functions

        //- Insert range before specified insertion index, by copying up
        void insertBefore(const label insert, const labelRange& range);

        //- Purge empty ranges, by copying down
        void purgeEmpty();

        //- Print the range for debugging purposes
        Ostream& printRange(Ostream& os, const labelRange& range) const;


public:

    // Constructors

        //- Construct null
        inline labelRanges();

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

        //- An STL const_iterator
        class const_iterator
        {
            friend class labelRanges;

          // Private data

            //- Reference to the list for which this is an iterator
            const labelRanges* list_;

            //- Current list-index
            label index_;

            //- Index of current element at list-index
            label subIndex_;

          // Constructors

            //- Construct from ranges at given index.
            //  A negative index signals the 'end' position
            inline const_iterator(const labelRanges* lst, const label i);

        public:

          // Member operators

            inline bool operator==(const const_iterator& iter) const;
            inline bool operator!=(const const_iterator& iter) const;

            //- Return the current label
            inline label operator*();

            inline const_iterator& operator++();
            inline const_iterator operator++(int);
        };


        //- A const_iterator set to the beginning of the list
        inline const_iterator cbegin() const;

        //- A const_iterator set to beyond the end of the list
        inline const const_iterator cend() const;

        //- A const_iterator set to the beginning of the list
        inline const_iterator begin() const;

        //- A const_iterator set to beyond the end of the list
        inline const const_iterator end() const;


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
