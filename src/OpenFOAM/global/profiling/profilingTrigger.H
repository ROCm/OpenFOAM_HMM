/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2016 Bernhard Gschaider
     \\/     M anipulation  | Copyright (C) 2016-2107 OpenCFD Ltd.
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
    Foam::profilingTrigger

Description
    Triggers for starting/stopping code profiling.

SourceFiles
    profilingTrigger.C

\*---------------------------------------------------------------------------*/

#ifndef profilingTrigger_H
#define profilingTrigger_H

#include "clockTime.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class profilingInformation;

/*---------------------------------------------------------------------------*\
                      Class profilingTrigger Declaration
\*---------------------------------------------------------------------------*/

class profilingTrigger
{
    // Private Data Members

        //- The timer for the profiling information
        clockTime clock_;

        //- The profiling information
        profilingInformation *ptr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        profilingTrigger(const profilingTrigger&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const profilingTrigger&) = delete;


public:

    // Constructors

        //- Construct profiling with given description.
        //  Descriptions beginning with 'application::' are reserved for
        //  internal use.
        profilingTrigger(const char* name);

        //- Construct profiling with given description.
        //  Descriptions beginning with 'application::' are reserved for
        //  internal use.
        profilingTrigger(const string& name);


    //- Destructor
    ~profilingTrigger();


    // Member Functions

      // Access

        //- True if the triggered profiling is active
        bool running() const;


      // Edit

        //- Stop triggered profiling
        void stop();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macros

//- Define profiling with specified name and description string
//  This is required if the description contains space, colons etc.
//  \sa addProfiling0
//  \sa endProfiling
#define addProfiling(name,descr)                                               \
    ::Foam::profilingTrigger  profilingTriggerFor##name(descr)

//- Define profiling with specified name and description correspond to the name
//  \sa addProfiling
//  \sa endProfiling
#define addProfiling0(name)                                                    \
    ::Foam::profilingTrigger  profilingTriggerFor##name(#name)

//- Define profiling with specified name and description correspond to the
//  compiler-defined function name string:
//  \sa addProfiling
//  \sa endProfiling
#ifdef __GNUC__
    #define addProfilingInFunction(name)                                       \
    ::Foam::profilingTrigger  profilingTriggerFor##name(__PRETTY_FUNCTION__)
#else
    #define addProfilingInFunction(name)                                       \
    ::Foam::profilingTrigger  profilingTriggerFor##name(__func__)
#endif

//- Remove profiling with specified name
//  \sa addProfiling
//  \sa addProfiling0
#define endProfiling(name)    profilingTriggerFor##name.stop()


#endif

// ************************************************************************* //
