/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2011 Symscape
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
#include "UList.H"

#include <float.h>  // For *fp functions
#include <limits>

// File-local functions
#include "signalMacros.C"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::sigFpe::switchFpe_(Foam::debug::optimisationSwitch("trapFpe", 0));
bool Foam::sigFpe::switchNan_(Foam::debug::optimisationSwitch("setNaN", 0));

bool Foam::sigFpe::sigActive_ = false;
bool Foam::sigFpe::nanActive_ = false;

// Saved old FPE signal trapping setting (file-local variable)
static unsigned int oldFpe_ = 0u;


static void clearFpe()
{
    #ifndef Foam_no_sigFpe
    _clearfp();
    _controlfp(oldFpe_, 0xFFFFFFFF);
    #endif
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Can turn on/off via env variable containing a bool (true|false|on|off ...)
// or by the specified flag
static bool isTrue(const char* envName, bool deflt)
{
    Foam::Switch sw(Foam::Switch::find(Foam::getEnv(envName)));

    if (sw.good())
    {
        return static_cast<bool>(sw);
    }

    // Env was not set or did not contain a valid bool value
    return deflt;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigFpe::sigHandler(int)
{
    resetHandler("SIGFPE", SIGFPE);

    JobInfo::shutdown();        // From running -> finished
    error::printStack(Perr);
    clearFpe();
    ::raise(SIGFPE);            // Throw signal (to old handler)
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    set(false);
}


Foam::sigFpe::ignore::ignore()
:
    wasActive_(sigFpe::active())
{
    if (wasActive_)
    {
        sigFpe::unset();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    unset(false);
}


Foam::sigFpe::ignore::~ignore()
{
    restore();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigFpe::ignore::restore()
{
    if (wasActive_)
    {
        sigFpe::set();
    }
    wasActive_ = false;
}


bool Foam::sigFpe::requested()
{
    return isTrue("FOAM_SIGFPE", switchFpe_);
}


void Foam::sigFpe::set(bool verbose)
{
    if (!sigActive_ && requested())
    {
        #ifdef Foam_no_sigFpe

        if (verbose)
        {
            Info<< "trapFpe: Floating point exception trapping ";
            Info<< "- disabled on this platform" << endl;
        }

        #else

        oldFpe_ = _controlfp(0, 0);

        const unsigned int newFpe =
        (
            oldFpe_ & ~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW)
        );

        _controlfp(newFpe, _MCW_EM);

        setHandler("SIGFPE", SIGFPE, sigHandler);

        sigActive_ = true;

        if (verbose)
        {
            Info<< "trapFpe: Floating point exception trapping ";

            if (sigActive_)
            {
                Info<< "enabled (FOAM_SIGFPE)." << endl;
            }
            else
            {
                Info<< "- not supported on this platform" << endl;
            }
        }
        #endif
    }


    nanActive_ = false;
    if (isTrue("FOAM_SETNAN", switchNan_))
    {
        if (verbose)
        {
            Info<< "setNaN : Initialise allocated memory to NaN "
                << "- not supported on this platform" << endl;
        }
    }
}


void Foam::sigFpe::unset(bool verbose)
{
    if (sigActive_)
    {
        if (verbose)
        {
            Info<< "sigFpe : Disabling floating point exception trapping"
                << endl;
        }

        sigActive_ = false;

        clearFpe();

        resetHandler("SIGFPE", SIGFPE);
    }

    nanActive_ = false;
}


void Foam::sigFpe::fillNan(UList<scalar>& list)
{
    list = std::numeric_limits<scalar>::signaling_NaN();
}


// ************************************************************************* //
