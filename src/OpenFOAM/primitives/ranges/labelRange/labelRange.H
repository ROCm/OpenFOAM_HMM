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
    A range of labels defined by a start and a size.

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


    // Constructors

        //- An empty range with zero for start/size.
        inline labelRange();

        //- Construct a range from start and size, enforcing non-negative size.
        //  Optionally adjust the start to avoid any negative indices.
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

        //- Is the range non-empty?
        inline bool valid() const;

        //- The (inclusive) lower value of the range
        inline label start() const;

        //- The effective size of the range
        inline label size() const;

        //- The (inclusive) lower value of the range - same as start
        inline label first() const;

        //- The (inclusive) upper value of the range
        inline label last() const;

        //- Reset start and size, enforcing non-negative size.
        //  Optionally adjust the start to avoid any negative indices.
        //  Return true if the updated range valid (non-empty).
        inline bool reset
        (
            const label start,
            const label size,
            const bool adjustStart = false
        );

        //- Return true if the value is located the range
        inline bool found(const label value) const;

        //- Return true if the ranges overlap.
        //  Optional test for ranges that also just touch each other
        bool overlaps(const labelRange& range, bool touches=false) const;

        //- Return a joined range, squashing any gaps in between
        //  A prior overlaps() check can be used to avoid squashing gaps.
        labelRange join(const labelRange& range) const;

        //- Calculate the intersection of the range with another.
        //  If there is no intersection, it returns an empty range with zero
        //  for start/size.
        labelRange subset(const labelRange& range) const;

        //- Calculate the intersection of the range with another.
        //  If there is no intersection, it returns an empty range with zero
        //  for start/size.
        labelRange subset(const label start, const label size) const;


    // Member Operators

        //- Return element in the range, no bounds checking
        inline label operator[](const label localIndex) const;


    // STL iterator

        //- An STL const_iterator
        class const_iterator
        {
            friend class labelRange;

          // Private data

            //- Reference to the range for which this is an iterator
            const labelRange* range_;

            //- Current index, relative to the start
            label index_;

          // Constructors

            //- Construct from range at given index.
            //  A negative index signals the 'end' position
            inline const_iterator(const labelRange* range, const label i);

        public:

          // Member operators

            inline bool operator==(const const_iterator& iter) const;
            inline bool operator!=(const const_iterator& iter) const;

            //- Return the current label
            inline label operator*() const;

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


// Global Operators

inline bool operator==(const labelRange& a, const labelRange& b)
{
    return (a.first() == b.first() && a.size() == b.size());
}

inline bool operator!=(const labelRange& a, const labelRange& b)
{
    return !(a == b);
}


//- Comparison function for sorting, compares the start.
//  If the start values are equal, also compares the size.
inline bool operator<(const labelRange& a, const labelRange& b)
{
    return
    (
        a.first() < b.first()
     ||
        (
            !(b.first() < a.first())
         && a.size() < b.size()
        )
    );
}

inline bool operator<=(const labelRange& a, const labelRange& b)
{
    return !(b < a);
}


inline bool operator>(const labelRange& a, const labelRange& b)
{
    return (b < a);
}

inline bool operator>=(const labelRange& a, const labelRange& b)
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelRangeI.H"

#endif

// ************************************************************************* //
