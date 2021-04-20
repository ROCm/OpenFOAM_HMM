/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::function1Base

Description
    Top level data entry class for use in dictionaries. Provides a mechanism
    to specify a variable as a certain type, e.g. constant or time varying, and
    provide functions to return the (interpolated) value, and integral between
    limits.

    Extends the Function1 class by adding autoMap and rMap functions

SourceFiles
    function1Base.C

\*---------------------------------------------------------------------------*/

#ifndef function1Base_H
#define function1Base_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                     Class function1Base Declaration
\*---------------------------------------------------------------------------*/

class function1Base
:
    public refCount
{
protected:

    // Protected Data

        //- Name of entry
        const word name_;


    // Protected Member Functions

        //- No copy assignment
        void operator=(const function1Base&) = delete;


public:

    // Constructors

        //- Construct from entry name
        explicit function1Base(const word& entryName);

        //- Construct from entry name and dictionary (unused)
        function1Base(const word& entryName, const dictionary& dict);

        //- Copy construct
        explicit function1Base(const function1Base& rhs);


    //- Destructor
    virtual ~function1Base() = default;


    // Member Functions

    // Access

        //- The name of the entry
        const word& name() const
        {
            return name_;
        }


    // Manipulation

        //- Convert time
        virtual void convertTimeBase(const Time& t);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
