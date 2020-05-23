/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::functionObjects::parProfiling

Group
    grpUtilitiesFunctionObjects

Description
    Simple (simplistic) mpi-profiling.

Usage
    Example of function object specification:
    \verbatim
    profiling
    {
        type  parProfiling;
        libs  (utilityFunctionObjects);

        // Report stats on exit only (instead of every time step)
        executeControl  onEnd;
        writeControl    none;
    }
    \endverbatim

SourceFiles
    parProfiling.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_parProfiling_H
#define functionObjects_parProfiling_H

#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class parProfiling Declaration
\*---------------------------------------------------------------------------*/

class parProfiling
:
    public functionObject
{
    // Private Member Functions

        //- No copy construct
        parProfiling(const parProfiling&) = delete;

        //- No copy assignment
        void operator=(const parProfiling&) = delete;


public:

    //- Runtime type information
    TypeName("parProfiling");


    // Constructors

        //- Construct from Time and dictionary
        parProfiling
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~parProfiling();


    // Member Functions

        //- Report the current profiling information
        void report();

        //- Report
        virtual bool execute();

        //- Do nothing
        virtual bool write();

        //- Report
        virtual bool end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
