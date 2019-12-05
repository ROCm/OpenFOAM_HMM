/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::clock

Description
    Read access to the system clock with formatting.

SourceFiles
    clock.C

\*---------------------------------------------------------------------------*/

#ifndef clock_H
#define clock_H

#include <ctime>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clock Declaration
\*---------------------------------------------------------------------------*/

class clock
{
    // Private data

        //- Time structure used
        typedef time_t value_type;

        //- Start time in seconds, at the time of construction
        value_type start_;

        //- Last time when elapsedClockTime or clockTimeIncrement was called
        mutable value_type last_;

public:

    // Constructors

        //- Construct null, storing the start time
        clock();


    // Static Member Functions

        //- Get the current clock time in seconds
        static time_t getTime();

        //- Return the current wall-clock date as a raw struct
        static const struct tm rawDate();

        //- Return the current wall-clock date/time as a string
        //  format according to ISO-8601 (yyyy-mm-ddThh:mm:ss)
        static std::string dateTime();

        //- Return the current wall-clock date as a string
        static std::string date();

        //- Return the current wall-clock time as a string
        static std::string clockTime();


    // Member Functions

        //- Returns wall-clock time from clock instantiation
        double elapsedClockTime() const;

        //- Returns wall-clock time from last call of clockTimeIncrement()
        double clockTimeIncrement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
