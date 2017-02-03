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
    Foam::labelRange

Description
    A range of labels.

SourceFiles
    labelRange.C

\*---------------------------------------------------------------------------*/
#ifndef labelRange_H
#define labelRange_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class labelRange;
Istream& operator>>(Istream& is, labelRange& range);
Ostream& operator<<(Ostream& os, const labelRange& range);

/*---------------------------------------------------------------------------*\
                         Class labelRange Declaration
\*---------------------------------------------------------------------------*/

class labelRange
{
    // Private data

        label start_;
        label size_;

public:

    static int debug;


    // Public classes

        //- Less function class for sorting labelRange
        class less
        {
        public:

            bool operator()(const labelRange& a, const labelRange& b)
            {
                return a.operator<(b);
            }
        };


    // Constructors

        //- Construct an empty range with zero as start and size.
        inline labelRange();

        //- Construct a range from start and size.
        //  Optionally adjust the start to avoid any negative indices.
        //  Always reduce a negative size to zero.
        inline labelRange
        (
            const label start,
            const label size,
            const bool adjustStart = false
        );

        //- Construct from Istream.
        labelRange(Istream& is);


    // Member Functions

        //- Alias for setSize(const label)
        inline void resize(const label n);

        //- Adjust size
        inline void setSize(const label n);

        //- Reset to zero start and zero size
        inline void clear();

        //- Is the range empty?
        inline bool empty() const;

        //- Adjust the start to avoid any negative indices
        void adjust();

        //- Is the range valid (non-empty)?
        inline bool valid() const;

        //- Return the effective size of the range
        inline label size() const;

        //- The (inclusive) lower value of the range
        inline label first() const;

        //- The (inclusive) upper value of the range
        inline label last() const;

        //- Reset start and size.
        //  Optionally adjust the start to avoid any negative indices.
        //  Always reduce a negative size to zero.
        //  Return true if the updated range valid (non-empty).
        inline bool reset
        (
            const label start,
            const label size,
            const bool adjustStart = false
        );

        //- Return true if the value is within the range
        inline bool contains(const label value) const;

        //- Return true if the ranges overlap.
        //  Optional test for ranges that also just touch each other
        bool overlaps(const labelRange& range, bool touches=false) const;

        //- Return a joined range, squashing any gaps in between
        //  A prior overlaps() check can be used to avoid squashing gaps.
        labelRange join(const labelRange& range) const;


    // Member Operators

        //- Return element in range, no bounds checking
        inline label operator[](const label i) const;

        //- Comparison function for sorting, compares the start.
        //  If the start values are equal, also compares the size.
        inline bool operator<(const labelRange& rhs) const;

        //- Join ranges, squashing any gaps in between
        //  A prior overlaps() check can be used to avoid squashing gaps.
        void operator+=(const labelRange& rhs);

        inline bool operator==(const labelRange& rhs) const;
        inline bool operator!=(const labelRange& rhs) const;


    // STL iterator

        //- An STL const_iterator
        class const_iterator
        {
            friend class labelRange;

            // Private data

                //- Reference to the range for which this is an iterator
                const labelRange& range_;

                //- Current index
                label index_;


            // Constructors

                //- Construct from range at 'begin' or 'end' position
                inline const_iterator
                (
                    const labelRange& range,
                    const bool endIter = false
                );

        public:

            // Member operators

                inline bool operator==(const const_iterator& iter) const;
                inline bool operator!=(const const_iterator& iter) const;

                inline label operator*() const;
                inline label operator()() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);
        };


        //- A const_iterator set to the beginning of the range
        inline const_iterator cbegin() const;

        //- A const_iterator set to beyond the end of the range
        inline const const_iterator cend() const;

        //- A const_iterator set to the beginning of the range
        inline const_iterator begin() const;

        //- A const_iterator set to beyond the end of the range
        inline const const_iterator end() const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, labelRange& range);
        friend Ostream& operator<<(Ostream& os, const labelRange& range);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelRangeI.H"

#endif

// ************************************************************************* //
