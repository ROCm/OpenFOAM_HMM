/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::clockValue

Description
    Access to high-resolution clock value with some basic operations.

SourceFiles
    clockValue.C

\*---------------------------------------------------------------------------*/

#ifndef clockValue_H
#define clockValue_H

#include <chrono>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class clockValue Declaration
\*---------------------------------------------------------------------------*/

class clockValue
{
    // Private Data

        //- Time structure used
        typedef std::chrono::high_resolution_clock::duration value_type;

        value_type value_;


public:

    // Constructors

        //- Construct zero initialized
        clockValue();

        //- Construct from duration with the same clock base
        explicit clockValue(const value_type& value);

        //- Construct zero initialized or with current time
        explicit clockValue(bool useNow);


    // Factory Methods

        //- Return the current time value from the system
        inline static clockValue now()
        {
            return clockValue(true);
        }


    // Member Functions

        //- Return the value
        inline const value_type& value() const
        {
            return value_;
        }

        //- Reset to zero
        void clear();

        //- Update to the current now() time from the system
        void update();

        //- The time elapsed from now() since the start time.
        clockValue elapsed() const;

        //- The value in seconds (rounded)
        long seconds() const;

        //- Format as day-hh:mm:ss string
        std::string str() const;


    // Operators

        //- Conversion operator to seconds in floating point
        operator double() const;

        //- Subtract time value
        clockValue& operator-=(const clockValue& rhs);

        //- Add time value
        clockValue& operator+=(const clockValue& rhs);
};


// Global Operators

//- Subtraction of clock values
clockValue operator-(const clockValue& a, const clockValue& b);

//- Addition of clock values
clockValue operator+(const clockValue& a, const clockValue& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
