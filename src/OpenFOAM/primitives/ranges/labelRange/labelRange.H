/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    A range or interval of labels defined by a start and a size.

SourceFiles
    labelRange.C

\*---------------------------------------------------------------------------*/
#ifndef labelRange_H
#define labelRange_H

#include "label.H"
#include <iterator>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class labelRange;
class Istream;
class Ostream;
template<class T> class List;

Istream& operator>>(Istream& is, labelRange& range);
Ostream& operator<<(Ostream& os, const labelRange& range);

/*---------------------------------------------------------------------------*\
                         Class labelRange Declaration
\*---------------------------------------------------------------------------*/

class labelRange
{
    // Private Data

        //- The start point for the interval
        label start_;

        //- The length of the interval
        label size_;

public:

    // STL type definitions

        //- Type of values the range contains
        typedef label value_type;

        //- The type that can represent the size of the range
        typedef label size_type;

        //- Forward iterator with const access
        class const_iterator;


    // Static Data Members

        static int debug;

        //- An empty range with start=0, size=0.
        static const labelRange null;


    // Constructors

        //- An empty range with zero for start/size.
        inline constexpr labelRange() noexcept;

        //- Construct a range from start/size, enforcing non-negative size.
        inline labelRange(const label start, const label size) noexcept;

        //- Construct a range from start/size, enforcing non-negative size.
        //  Optionally adjust the start to avoid any negative indices.
        inline labelRange
        (
            const label start,
            const label size,
            const bool adjustStart
        ) noexcept;

        //- Construct from Istream.
        labelRange(Istream& is);


    // Static Member Functions

        //- An identity range with range[i] == i.
        inline static labelRange identity(const label len) noexcept;


    // Member Functions

        //- Change the start position
        inline void setStart(const label i) noexcept;

        //- Change the size
        inline void resize(const label n) noexcept;

        //- Change the size - alias for resize()
        inline void setSize(const label n) noexcept;

        //- Decrease the size by 1, but never below 0.
        inline void decrement() noexcept;

        //- Increase the size by 1.
        inline void increment() noexcept;

        //- Reset to zero start and zero size
        inline void clear() noexcept;

        //- Is the range empty?
        inline bool empty() const noexcept;

        //- Return the range as a list of labels
        List<label> labels() const;

        //- Adjust the start to avoid negative indices.
        //  The size is decreased accordingly, but will never become negative.
        //  Eg, adjusting (-10, 15) becomes (0,5).
        //  adjusting (-20, 15) becomes (0,0)
        void adjust() noexcept;

        //- Is the range non-empty?
        inline bool valid() const noexcept;

        //- The (inclusive) lower value of the range
        inline label start() const noexcept;

        //- The effective size of the range
        inline label size() const noexcept;

        //- The (inclusive) lower value of the range - same as start
        inline label first() const noexcept;

        //- The (inclusive) upper value of the range
        inline label last() const noexcept;

        //- The value before the start of the range
        inline label before() const noexcept;

        //- The value after the last element in the range
        inline label after() const noexcept;

        //- Reset start and size, enforcing non-negative size.
        //  \return true if the updated range is valid (non-empty).
        inline bool reset(const label start, const label size) noexcept;

        //- Reset start and size, enforcing non-negative size.
        //  Optionally adjust the start to avoid any negative indices.
        //  \return true if the updated range is valid (non-empty).
        inline bool reset
        (
            const label start,
            const label size,
            const bool adjustStart
        ) noexcept;

        //- Return true if the (global) value is located within the range
        inline bool found(const label value) const noexcept;

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

        //- Calculate the intersection with the given start/size range.
        //  If there is no intersection, it returns an empty range with zero
        //  for start/size.
        labelRange subset(const label start, const label size) const;

        //- Calculate the intersection with the given 0/size range.
        //  If there is no intersection, it returns an empty range with zero
        //  for start/size.
        labelRange subset0(const label size) const;


        //- Return const_iterator to element in the range,
        //- with bounds checking.
        //  \return iterator at the requested position, or end() if it is
        //      out of bounds
        inline const_iterator at(const label localIndex) const;


    // Member Operators

        //- Return element in the range, without bounds checking
        inline label operator[](const label localIndex) const noexcept;

        //- Return true if the (global) value is located within the range.
        //  Behaviour identical to found() - usable as a predicate
        inline bool operator()(const label value) const noexcept;

        //- Increase the size by 1.
        inline label operator++() noexcept;
        inline label operator++(int) noexcept;

        //- Decrease the size by 1, but never below 0.
        inline label operator--() noexcept;
        inline label operator--(int) noexcept;


    // STL iterator

        //- Forward iterator with const access
        class const_iterator
        :
            public std::iterator
            <
                std::input_iterator_tag,
                label,
                label,
                const label*,
                const label&
            >
        {
            //- The current (global) index value
            label value_;

        public:

        // Constructors

            //- Construct from range at given local index.
            //  If the local index is out of range (eg, negative),
            //  this creates the 'end' iterator
            inline const_iterator(const labelRange* range, const label i=0);

        // Member Operators

            //- Return the current (global) value
            inline label operator*() const;

            inline const_iterator& operator++();
            inline const_iterator operator++(int);

            inline bool operator==(const const_iterator& iter) const;
            inline bool operator!=(const const_iterator& iter) const;
        };


        //- A const_iterator set to the beginning of the range
        //  The value returned is guaranteed to be the same as start()
        inline const_iterator begin() const;

        //- A const_iterator set to the beginning of the range
        //  The value returned is guaranteed to be the same as start()
        inline const_iterator cbegin() const;

        //- A const_iterator set to beyond the end of the range
        inline const const_iterator cend() const;

        //- A const_iterator set to beyond the end of the range
        inline const const_iterator end() const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, labelRange& range);
        friend Ostream& operator<<(Ostream& os, const labelRange& range);

};


// Global Functions

//- Conversion/extraction to labelRange operation (functor).
//  Specializations shall provide a corresponding \c operator().
//  For example,
//  \code
//  template<>
//  struct labelRangeOp<polyPatch>
//  {
//      labelRange operator()(const polyPatch& pp) const
//      {
//          return labelRange(pp.start(), pp.size());
//      }
//  };
//  \endcode
template<class> struct labelRangeOp;


// Global Operators

inline bool operator==(const labelRange& a, const labelRange& b) noexcept
{
    return (a.first() == b.first() && a.size() == b.size());
}

inline bool operator!=(const labelRange& a, const labelRange& b) noexcept
{
    return (a.first() != b.first() || a.size() != b.size());
}


//- Comparison function for sorting, compares the start.
//  If the start values are equal, also compares the size.
inline bool operator<(const labelRange& a, const labelRange& b) noexcept
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

inline bool operator<=(const labelRange& a, const labelRange& b) noexcept
{
    return !(b < a);
}


inline bool operator>(const labelRange& a, const labelRange& b) noexcept
{
    return (b < a);
}

inline bool operator>=(const labelRange& a, const labelRange& b) noexcept
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelRangeI.H"

#endif

// ************************************************************************* //
