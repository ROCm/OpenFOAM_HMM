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
    Foam::timeSeries

Description
    A list of times and values.
    The time values must be positive and monotonically increasing.

    The treatment of out-of-bounds values depends on the current setting
    of bounding.

    If @a REPEAT bounding is in effect, the final time value is treated
    as being equivalent to time=0 for the following periods.

Note
    - Accessing an empty list will result in an error.
    - Accessing a list with a single element will always return the same value.

SourceFiles
    timeSeries.C

\*---------------------------------------------------------------------------*/

#ifndef timeSeries_H
#define timeSeries_H

#include "List.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class timeSeries Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class timeSeries
:
    public List<Tuple2<scalar, T> >
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound times
        enum bounds
        {
            ERROR,          /*!< Exit with a FatalError */
            WARN,           /*!< Issue warning and clamp value (default) */
            CLAMP,          /*!< Clamp value to the start/end value */
            REPEAT          /*!< Treat as a repeating list */
        };

private:

    // Private data

        //- Enumeration for handling out-of-bound times
        bounds bounding_;

public:

    // Constructors

        //- Construct null, optionally with a given bounding
        timeSeries(const bounds = timeSeries::WARN);

        //- Construct null with a given bounding
        timeSeries(const word&);

        //- Construct from Istream, optionally with a given bounding
        timeSeries(Istream&, const bounds = timeSeries::WARN);

        //- Construct from Istream with a given bounding
        timeSeries(Istream&, const word&);

    // Destructor

        ~timeSeries();

    // Member Functions

    // Access

        //- Return the size
        label size() const
        {
            return List<Tuple2<scalar, T> >::size();
        }

        //- Return the out-of-bounds treatment as a word
        word bounding() const;

    // Check

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void check() const;

    // Edit

        //- Set the out-of-bounds treatment from enum, return previous setting
        bounds bounding(const bounds& bound)
        {
            bounds prev = bounding_;
            bounding_ = bound;
            return prev;
        }

        //- Set the out-of-bounds treatment from word
        void bounding(const word& bound);

    // Member Operators

        //- Return an element of constant Tuple2<scalar, T>
        const Tuple2<scalar, T>& operator[](const label) const;

        //- Return an interpolated value
        T operator()(const scalar) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "timeSeries.C"
#endif

#endif
// ************************************************************************* //
