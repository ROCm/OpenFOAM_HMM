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
    Foam::interpolationTable

Description
    A list of times and values.
    The time values must be positive monotonically increasing.

    The handling of out-of-bounds values depends on the current setting
    of @a outOfBounds.

    If @a REPEAT is chosen for the out-of-bounds handling, the final time
    value is treated as being equivalent to time=0 for the following periods.

Note
    - Accessing an empty list results in an error.
    - Accessing a list with a single element will always return the same value.

SourceFiles
    interpolationTable.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationTable_H
#define interpolationTable_H

#include "List.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;

/*---------------------------------------------------------------------------*\
                    Class interpolationTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationTable
:
    public List<Tuple2<scalar, Type> >
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,          /*!< Exit with a FatalError */
            WARN,           /*!< Issue warning and clamp value (default) */
            CLAMP,          /*!< Clamp value to the start/end value */
            REPEAT          /*!< Treat as a repeating list */
        };


private:

    // Private data

        //- Enumeration for handling out-of-bound values
        boundsHandling boundsHandling_;

        //- File name
        fileName fileName_;


public:

    // Constructors

        //- Construct null
        interpolationTable();

        //- Construct from objectRegistry and dictionary
        interpolationTable(const objectRegistry& obr, const dictionary& dict);

        //- Construct copy
        interpolationTable(const interpolationTable& interpTable);


    //- Destructor
    ~interpolationTable();


    // Member Functions

        // Access

            //- Return the size
            label size() const
            {
                return List<Tuple2<scalar, Type> >::size();
            }

            //- Return the out-of-bounds handling as a word
            word boundsHandlingToWord(const boundsHandling& bound) const;

            //- Return the out-of-bounds handling as an enumeration
            boundsHandling wordToBoundsHandling(const word& bound) const;


        // Check

            //- Check that list is monotonically increasing
            //  Exit with a FatalError if there is a problem
            void check() const;


        // Edit

            //- Set the out-of-bounds handling from enum,
            //  return previous setting
            boundsHandling outOfBounds(const boundsHandling& bound);


        // Member Operators

            //- Return an element of constant Tuple2<scalar, Type>
            const Tuple2<scalar, Type>& operator[](const label) const;

            //- Return an interpolated value
            Type operator()(const scalar) const;


        // I-O

            //- Write
            void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interpolationTable.C"
#endif

#endif

// ************************************************************************* //
