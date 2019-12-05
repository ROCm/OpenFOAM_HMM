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
    Foam::clockTime

Description
    Starts timing (using rtc) and returns elapsed time from start.
    Better resolution (2uSec instead of ~20mSec) than cpuTime.

SourceFiles
    clockTime.C

\*---------------------------------------------------------------------------*/

#ifndef clockTime_H
#define clockTime_H

#include "clockValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class clockTime Declaration
\*---------------------------------------------------------------------------*/

class clockTime
{
    // Private Data

        //- The values being used.
        typedef clockValue value_type;

        //- Start time, at the time of construction
        value_type start_;

        //- Last time when elapsedTime or timeIncrement was called
        mutable value_type last_;


public:

    // Constructors

        //- Construct with the current clock time for the starting value
        clockTime();

        //- Construct with the given clock value for the starting value
        clockTime(const clockValue& clockval);


    // Member Functions

        //- Reset to use the current time for the start time
        void resetTime();

        //- Return time (in seconds) from the start
        double elapsedTime() const;

        //- Return time (in seconds) since last call to timeIncrement()
        double timeIncrement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
