/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "sigFpe.H"
#include "error.H"
#include "JobInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"
#include "Switch.H"

#ifdef LINUX_GNUC
    #ifndef __USE_GNU
        #define __USE_GNU
    #endif
    #include <fenv.h>
    #include <malloc.h>
#elif defined(sgiN32) || defined(sgiN32Gcc)
    #include <sigfpe.h>
#endif

#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::sigFpe::switchFpe_(Foam::debug::optimisationSwitch("trapFpe", 0));
bool Foam::sigFpe::switchNan_(Foam::debug::optimisationSwitch("setNaN", 0));

bool Foam::sigFpe::sigActive_ = false;
bool Foam::sigFpe::mallocNanActive_ = false;

struct sigaction Foam::sigFpe::oldAction_;


// File-scope function.
// Controlled by env variable containing a bool (true|false|on|off ...)
// or by the specified flag
static bool isTrue(const char* envName, const bool flag)
{
    const std::string str = Foam::getEnv(envName);

    if (str.size())
    {
        Foam::Switch sw(str, true); // silently ignore bad input
        if (sw.valid())
        {
            return sw;
        }
    }

    // Env was not set or did not contain a valid bool value
    return flag;
}


void Foam::sigFpe::fillNan(UList<scalar>& lst)
{
    lst = std::numeric_limits<scalar>::signaling_NaN();
}


#ifdef LINUX
extern "C"
{
    extern void* __libc_malloc(size_t size);

    // Override the GLIBC malloc to support mallocNan
    void* malloc(size_t size)
    {
        if (Foam::sigFpe::mallocNanActive_)
        {
            return Foam::sigFpe::mallocNan(size);
        }
        else
        {
            return __libc_malloc(size);
        }
    }
}

void* Foam::sigFpe::mallocNan(size_t size)
{
    // Call the low-level GLIBC malloc function
    void * result = __libc_malloc(size);

    // Initialize to signalling NaN
    UList<scalar> lst(reinterpret_cast<scalar*>(result), size/sizeof(scalar));
    sigFpe::fillNan(lst);

    return result;
}
#endif


#ifdef LINUX_GNUC
void Foam::sigFpe::sigHandler(int)
{
    // Reset old handling
    if (sigaction(SIGFPE, &oldAction_, nullptr) < 0)
    {
        FatalErrorInFunction
            << "Cannot reset SIGFPE trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo.signalEnd();

    error::printStack(Perr);

    // Throw signal (to old handler)
    raise(SIGFPE);
}
#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    set(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    unset(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sigFpe::requested()
{
    return isTrue("FOAM_SIGFPE", switchFpe_);
}


void Foam::sigFpe::set(const bool verbose)
{
    if (!sigActive_ && requested())
    {
        bool supported = false;

        #ifdef LINUX_GNUC
        supported = true;

        feenableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        struct sigaction newAction;
        newAction.sa_handler = sigHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(SIGFPE, &newAction, &oldAction_) < 0)
        {
            FatalErrorInFunction
                << "Cannot set SIGFPE trapping"
                << abort(FatalError);
        }

        sigActive_ = true;

        #elif defined(sgiN32) || defined(sgiN32Gcc)
        supported = true;

        sigfpe_[_DIVZERO].abort=1;
        sigfpe_[_OVERFL].abort=1;
        sigfpe_[_INVALID].abort=1;

        sigfpe_[_DIVZERO].trace=1;
        sigfpe_[_OVERFL].trace=1;
        sigfpe_[_INVALID].trace=1;

        handle_sigfpes
        (
            _ON,
            _EN_DIVZERO
          | _EN_INVALID
          | _EN_OVERFL,
            0,
            _ABORT_ON_ERROR,
            nullptr
        );

        sigActive_ = true;

        #endif


        if (verbose)
        {
            Info<< "trapFpe: Floating point exception trapping ";

            if (supported)
            {
                Info<< "enabled (FOAM_SIGFPE)." << endl;
            }
            else
            {
                Info<< "- not supported on this platform" << endl;
            }
        }
    }


    if (isTrue("FOAM_SETNAN", switchNan_))
    {
        #ifdef LINUX
        mallocNanActive_ = true;
        #endif

        if (verbose)
        {
            Info<< "setNaN : Initialise allocated memory to NaN ";

            if (mallocNanActive_)
            {
                Info<< "enabled (FOAM_SETNAN)." << endl;
            }
            else
            {
                Info<< " - not supported on this platform" << endl;
            }
        }
    }
}


void Foam::sigFpe::unset(const bool verbose)
{
    #ifdef LINUX_GNUC
    // Reset signal
    if (sigActive_)
    {
        if (verbose)
        {
            Info<< "sigFpe : Disabling floating point exception trapping"
                << endl;
        }

        if (sigaction(SIGFPE, &oldAction_, nullptr) < 0)
        {
            FatalErrorInFunction
                << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }

        // Reset exception raising
        int oldExcept = fedisableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        if (oldExcept == -1)
        {
            FatalErrorInFunction
                << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }
        sigActive_ = false;
    }
    #endif

    #ifdef LINUX
    // Disable initialization to NaN
    mallocNanActive_ = false;
    #endif
}


// ************************************************************************* //
