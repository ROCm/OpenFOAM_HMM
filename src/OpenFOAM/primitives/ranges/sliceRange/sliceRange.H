/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::sliceRange

Description
    A set of labels defined by a start, a length and a stride.

SourceFiles
    sliceRange.C
    sliceRangeI.H

\*---------------------------------------------------------------------------*/
#ifndef sliceRange_H
#define sliceRange_H

#include "label.H"
#include <iterator>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Ostream;
template<class T> class List;
template<class T, unsigned N> class FixedList;

/*---------------------------------------------------------------------------*\
                         Class sliceRange Declaration
\*---------------------------------------------------------------------------*/

class sliceRange
{
protected:

    // Protected Data

        //- The start point for the interval
        label start_;

        //- The length of the interval
        label size_;

        //- The stride within the interval
        label stride_;


public:

    // STL type definitions

        //- Type of values the range contains
        typedef label value_type;

        //- The type that can represent the size of the range
        typedef label size_type;

        //- Forward iterator with const access
        class const_iterator;


    // Constructors

        //- Default copy construct
        sliceRange(const sliceRange&) = default;

        //- Default move construct
        sliceRange(sliceRange&&) = default;

        //- An empty slice (0,0,0)
        inline constexpr sliceRange() noexcept;

        //- Construct slice from start/size/stride,
        //- enforcing non-negative size and stride.
        inline sliceRange(label start, label size, label stride) noexcept;

        //- Construct slice from start/size/stride coefficients,
        //- enforcing non-negative size and stride.
        explicit sliceRange(const FixedList<label,3>& coeffs);


    // Member Functions

        //- Is the range empty?
        bool empty() const noexcept
        {
            return !size_;
        }

        //- Is the range non-empty?
        bool valid() const noexcept
        {
            return size_;
        }

        //- The (inclusive) lower value of the range
        constexpr label start() const noexcept
        {
            return start_;
        }

        //- The size of the range
        constexpr label size() const noexcept
        {
            return size_;
        }

        //- The stride for the range
        constexpr label stride() const noexcept
        {
            return stride_;
        }

        //- The (inclusive) lower value of the range - same as start()
        constexpr label first() const noexcept
        {
            return start_;
        }

        //- The (inclusive) upper value of the range
        label last() const noexcept
        {
            return start_ + (size_-1) * stride_;
        }

        //- Return the slice as a list of labels
        List<label> labels() const;


    // Member Operators

        //- Default copy assignment
        sliceRange& operator=(const sliceRange&) = default;

        //- Default move assignment
        sliceRange& operator=(sliceRange&&) = default;

        //- Return element in the range, without bounds checking
        inline label operator[](const label i) const noexcept;


    // Iterators

        //- A value indexer, for iteration or generation
        class indexer
        {
            //- The stride when indexing
            const label stride_;

            //- The global value
            label value_;

        public:

        // Constructors

            //- Construct from range at given local index.
            //  An out-of-range index (eg, negative) creates an 'end' iterator
            inline indexer(const sliceRange* range, const label i=0);


        // Member Functions

            //- Forward increment, no checking
            inline void next() noexcept;

            //- Forward increment, no checking
            inline void next(const label n) noexcept;

            //- Test for equality of values, not stride
            inline bool equals(const indexer& other) const noexcept;

        public:

        // Member Operators

            //- Return the value
            inline label operator*() const noexcept;

            //- Apply a postfix increment and return the current value.
            //  This operator definition is required for a generator -
            //  see std::generate()
            inline label operator()();
        };


        //- Forward iterator with const access
        class const_iterator
        :
            protected indexer,
            public std::iterator
            <
                std::input_iterator_tag,
                label,
                label,
                const label*,
                const label&
            >
        {
        public:

        // Constructors

            //- Construct from range at given local index.
            //  An out-of-range index (eg, negative) creates an 'end' iterator
            inline const_iterator(const sliceRange* range, const label i=0);


        // Member Operators

            //- Return the (global) value
            using indexer::operator*;

            //- Prefix increment, no checking
            inline const_iterator& operator++() noexcept;

            //- Arbitrary increment, no checking
            inline const_iterator& operator+=(const label n) noexcept;

            //- Prefix decrement, no checking
            inline const_iterator& operator--() noexcept;

            //- Arbitrary decrement, no checking
            inline const_iterator& operator-=(const label n) noexcept;

            //- Test for equality of values, not stride
            inline bool operator==(const const_iterator& iter) const noexcept;

            //- Test for inequality of values, not stride
            inline bool operator!=(const const_iterator& iter) const noexcept;
        };


        //- A const_iterator set to the beginning of the range
        //  The value returned is guaranteed to be the same as start()
        inline const_iterator begin() const;

        //- A const_iterator set to the beginning of the range
        //  The value returned is guaranteed to be the same as start()
        inline const_iterator cbegin() const;

        //- A const_iterator set to 1 beyond the end of the range.
        //  The value returned is the same as after()
        inline const const_iterator cend() const;

        //- A const_iterator set to 1 beyond the end of the range.
        //  The value returned is the same as after()
        inline const const_iterator end() const;

        //- Return a forward values generator
        inline indexer generator() const;

        //- Return const_iterator to a position within the range,
        //- with bounds checking.
        //  \return iterator at the requested position, or end() for
        //      out of bounds
        inline const_iterator at(const label i) const;
};


// IOstream Operators

//- Write sliceRange to Ostream as (start size stride) tuple
Ostream& operator<<(Ostream& os, const sliceRange& range);


// Global Operators

inline bool operator==(const sliceRange& a, const sliceRange& b) noexcept
{
    return
    (
        a.first() == b.first()
     && a.size() == b.size()
     && a.stride() == b.stride()
    );
}


inline bool operator!=(const sliceRange& a, const sliceRange& b) noexcept
{
    return !(a == b);
}


//- Comparison function for sorting, compares the start.
//  If the start values are equal, also compares the size.
//  If the sizes are equal, also compares the stride.
inline bool operator<(const sliceRange& a, const sliceRange& b) noexcept
{
    return
    (
        a.first() < b.first()
     ||
        (
            !(b.first() < a.first())
         &&
            (
                a.size() < b.size()
             ||
                (
                    !(b.size() < a.size())
                 && a.stride() < b.stride()
                )
            )
        )
    );

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sliceRangeI.H"

#endif

// ************************************************************************* //
