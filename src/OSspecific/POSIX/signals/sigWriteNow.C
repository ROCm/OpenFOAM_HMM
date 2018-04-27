/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "sigWriteNow.H"
#include "error.H"
#include "JobInfo.H"
#include "IOstreams.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Signal number to catch
int Foam::sigWriteNow::signal_
(
    Foam::debug::optimisationSwitch("writeNowSignal", -1)
);

Foam::Time* Foam::sigWriteNow::runTimePtr_ = nullptr;

struct sigaction Foam::sigWriteNow::oldAction_;


namespace Foam
{

// Register re-reader
class addwriteNowSignalToOpt
:
    public ::Foam::simpleRegIOobject
{

public:

    addwriteNowSignalToOpt(const char* name)
    :
        ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
    {}

    virtual ~addwriteNowSignalToOpt() = default;

    virtual void readData(Foam::Istream& is)
    {
        sigWriteNow::signal_ = readLabel(is);
        sigWriteNow::set(true);
    }

    virtual void writeData(Foam::Ostream& os) const
    {
        os << sigWriteNow::signal_;
    }
};

addwriteNowSignalToOpt addwriteNowSignalToOpt_("writeNowSignal");

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigWriteNow::sigHandler(int)
{
    if (runTimePtr_)
    {
        Info<< "sigWriteNow :"
            << " setting up write at end of the next iteration" << nl << endl;
        runTimePtr_->writeOnce();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigWriteNow::sigWriteNow()
{}


Foam::sigWriteNow::sigWriteNow(Time& runTime, bool verbose)
{
    runTimePtr_ = &runTime;     // Store runTime
    set(verbose);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigWriteNow::~sigWriteNow()
{
    // Reset old handling
    if (signal_ > 0)
    {
        if (sigaction(signal_, &oldAction_, nullptr) < 0)
        {
            FatalErrorInFunction
                << "Cannot reset " << signal_ << " trapping"
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigWriteNow::set(bool verbose)
{
    if (signal_ >= 0)
    {
        struct sigaction newAction;
        newAction.sa_handler = sigHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(signal_, &newAction, &oldAction_) < 0)
        {
            FatalErrorInFunction
                << "Cannot set " << signal_ << " trapping"
                << abort(FatalError);
        }

        if (verbose)
        {
            Info<< "sigWriteNow :"
                << " Enabling writing upon signal " << signal_
                << endl;
        }
    }
}


bool Foam::sigWriteNow::active() const
{
    return signal_ > 0;
}


// ************************************************************************* //
