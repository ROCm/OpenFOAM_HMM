/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2016 Bernhard Gschaider
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class profilingInformation;

/*---------------------------------------------------------------------------*\
                      Class profilingTrigger Declaration
\*---------------------------------------------------------------------------*/

class profilingTrigger
{
    // Private Data Members

        //- The profiling information
        profilingInformation *ptr_;


    // Private Member Functions

        //- No copy construct
        profilingTrigger(const profilingTrigger&) = delete;

        //- No copy assignment
        void operator=(const profilingTrigger&) = delete;


public:

    // Constructors

        //- Default construct, no profiling trigger
        profilingTrigger();

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

        //- True if the triggered profiling is active
        bool running() const;

        //- Stop triggered profiling
        void stop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macros

//- Define profiling trigger with specified name and description string
//  \sa endProfiling
#define addProfiling(name,descr)                                               \
    ::Foam::profilingTrigger  profilingTriggerFor##name(descr)

//- Define profiling trigger with specified name and description
//- corresponding to the compiler-defined function name string
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
#define endProfiling(name)    profilingTriggerFor##name.stop()


#endif

// ************************************************************************* //
