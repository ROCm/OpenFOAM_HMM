/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::sigFpe

Description
    Set up trapping for floating point exceptions (signal FPE).

    Defined by controlDict InfoSwitch entries:
    - \par trapFpe
      Enable floating point exception trapping.

    - \par setNaN
      Initialization all malloced memory to NaN.
      Combined with \c trapFpe, this causes usage of uninitialized scalars
      to trigger an abort.

    Environment variables:
      - \par FOAM_SIGFPE (true|false)
        overrides \c trapFpe
      - \par FOAM_SETNAN  (true|false)
        overrides \c setNaN

    Note that trapping can be set/removed through the static member functions
    or through the scope of the object (constructor sets trapping; destructor
    restores original). The class behaves as a singleton.

SourceFiles
    sigFpe.C

\*---------------------------------------------------------------------------*/

#ifndef sigFpe_H
#define sigFpe_H

#include <signal.h>
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class T> class UList;

/*---------------------------------------------------------------------------*\
                           Class sigFpe Declaration
\*---------------------------------------------------------------------------*/

class sigFpe
{
    // Private data

        //- Flag that floating point trapping should be used.
        //  Can override with FOAM_SIGFPE env variable
        static bool switchFpe_;

        //- Flag that NaN initialisation should be used.
        //  Can override with FOAM_SETNAN env variable
        static bool switchNan_;

        //- Flag to indicate floating point trapping is currently active
        static bool sigActive_;

        //- Flag to indicate mallocNan is currently active
        static bool nanActive_;

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;


    // Static data members

        #if defined(__linux__) && defined(__GNUC__)
        //- Handler for caught signals
        static void sigHandler(int);
        #endif


public:

    // Constructors

        //- Construct null
        sigFpe();


    //- Destructor
    ~sigFpe();


    // Member functions

        //- Check if SIGFPE signals handler is to be enabled.
        //  This is controlled by the trapFpe entry or the FOAM_SIGFPE
        //  environment variable
        static bool requested();

        //- True if SIGFPE handling is currently active.
        static inline bool active()
        {
            return sigActive_;
        }

        //- True if NaN memory initialisation is currently active.
        static inline bool nanActive()
        {
            return nanActive_;
        }

        //- Activate SIGFPE signal handler when FOAM_SIGFPE is %set
        //  Fill memory with NaN when FOAM_SETNAN is %set
        static void set(bool verbose=false);

        //- Deactivate SIGFPE signal handler and NaN memory initialisation
        static void unset(bool verbose=false);

        #ifdef __linux__
        //- Malloc function which initializes to NaN
        static void* mallocNan(size_t size);
        #endif

        //- Fill data block with NaN values
        static void fillNan(UList<scalar>& list);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
