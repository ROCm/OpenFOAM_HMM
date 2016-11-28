/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "noiseModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noiseModel, 0);
    defineRunTimeSelectionTable(noiseModel, dictionary);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::scalar Foam::noiseModel::checkUniformTimeStep
(
    const scalarList& times
) const
{
    scalar deltaT = -1.0;

    if (times.size() > 1)
    {
        for (label i = 1; i < times.size(); i++)
        {
            scalar dT = times[i] - times[i-1];

            if (deltaT < 0)
            {
                deltaT = dT;
            }

            if (mag(dT/deltaT - 1) > 1e-8)
            {
                FatalErrorInFunction
                    << "Unable to process data with a variable time step"
                    << exit(FatalError);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unable to create FFT with a single value"
            << exit(FatalError);
    }

    return deltaT;
}


Foam::label Foam::noiseModel::findStartTimeIndex
(
    const instantList& allTimes,
    const scalar startTime
) const
{
    forAll(allTimes, timeI)
    {
        const instant& i = allTimes[timeI];

        if (i.value() >= startTime)
        {
            return timeI;
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noiseModel::noiseModel(const dictionary& dict, const bool readFields)
:
    dict_(dict),
    rhoRef_(1),
    nSamples_(65536),
    fLower_(25),
    fUpper_(10000),
    startTime_(0),
    windowModelPtr_(),
    graphFormat_("raw")
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::noiseModel::~noiseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::noiseModel::read(const dictionary& dict)
{
    dict.readIfPresent("rhoRef", rhoRef_);
    dict.readIfPresent("N", nSamples_);
    dict.readIfPresent("fl", fLower_);
    dict.readIfPresent("fu", fUpper_);
    dict.readIfPresent("startTime", startTime_);
    dict.readIfPresent("graphFormat", graphFormat_);

    // Check number of samples  - must be a power of 2 for our FFT
    bool powerOf2 = ((nSamples_ != 0) && !(nSamples_ & (nSamples_ - 1)));
    if (!powerOf2)
    {
        FatalIOErrorInFunction(dict)
            << "N: Number of samples in sampling windows must be a "
            << "power of 2"
            << exit(FatalIOError);
    }

    if (fLower_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "fl: lower frequency bound must be greater than zero"
            << exit(FatalIOError);

    }

    if (fUpper_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "fu: upper frequency bound must be greater than zero"
            << exit(FatalIOError);

    }

    if (fUpper_ < fLower_)
    {
        FatalIOErrorInFunction(dict)
            << "fu: upper frequency bound must be greater than lower "
            << "frequency bound (fl)"
            << exit(FatalIOError);

    }

    windowModelPtr_ = windowModel::New(dict, nSamples_);

    return true;
}


// ************************************************************************* //
