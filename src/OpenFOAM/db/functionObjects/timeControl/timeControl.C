/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "timeControl.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::timeControl::timeControls
>
Foam::timeControl::timeControlNames_
{
    { timeControl::ocTimeStep, "timeStep" },
    { timeControl::ocWriteTime, "writeTime" },
    { timeControl::ocOutputTime, "outputTime" },
    { timeControl::ocAdjustableRunTime, "adjustableRunTime" },
    { timeControl::ocRunTime, "runTime" },
    { timeControl::ocClockTime, "clockTime" },
    { timeControl::ocCpuTime, "cpuTime" },
    { timeControl::ocOnEnd, "onEnd" },
    { timeControl::ocNone, "none" },
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeControl::timeControl
(
    const Time& t,
    const dictionary& dict,
    const word& prefix
)
:
    time_(t),
    prefix_(prefix),
    timeControl_(ocTimeStep),
    intervalSteps_(0),
    interval_(-1),
    executionIndex_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeControl::~timeControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::timeControl::entriesPresent
(
    const dictionary& dict,
    const word& prefix
)
{
    const word controlName(prefix + "Control");

    if (dict.found(controlName))
    {
        return true;
    }

    return false;
}


void Foam::timeControl::read(const dictionary& dict)
{
    word controlName(prefix_ + "Control");
    word intervalName(prefix_ + "Interval");

    // For backward compatibility support the deprecated 'outputControl' option
    // now superseded by 'writeControl' for compatibility with Time
    if (prefix_ == "write" && dict.found("outputControl"))
    {
        IOWarningInFunction(dict)
            << "Using deprecated 'outputControl'" << nl
            << "    Please use 'writeControl' with 'writeInterval'"
            << endl;

        // Change to the old names for this option
        controlName = "outputControl";
        intervalName = "outputInterval";
    }

    if (dict.found(controlName))
    {
        timeControl_ = timeControlNames_.lookup(controlName, dict);
    }
    else
    {
        timeControl_ = ocTimeStep;
    }

    switch (timeControl_)
    {
        case ocTimeStep:
        {
            intervalSteps_ = dict.lookupOrDefault<label>(intervalName, 0);
            break;
        }

        case ocWriteTime:
        case ocOutputTime:
        {
            intervalSteps_ = dict.lookupOrDefault<label>(intervalName, 1);
            break;
        }

        case ocClockTime:
        case ocRunTime:
        case ocCpuTime:
        case ocAdjustableRunTime:
        {
            const scalar userTime = readScalar(dict.lookup(intervalName));
            interval_ = time_.userTimeToTime(userTime);
            break;
        }

        case ocOnEnd:
        default:
        {
            break;
        }
    }
}


bool Foam::timeControl::execute()
{
    switch (timeControl_)
    {
        case ocTimeStep:
        {
            return
            (
                (intervalSteps_ <= 1)
             || !(time_.timeIndex() % intervalSteps_)
            );
            break;
        }

        case ocWriteTime:
        case ocOutputTime:
        {
            if (time_.writeTime())
            {
                executionIndex_++;
                return !(executionIndex_ % intervalSteps_);
            }
            break;
        }

        case ocRunTime:
        case ocAdjustableRunTime:
        {
            label executionIndex = label
            (
                (
                    (time_.value() - time_.startTime().value())
                  + 0.5*time_.deltaTValue()
                )
               /interval_
            );

            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case ocCpuTime:
        {
            label executionIndex = label
            (
                returnReduce(time_.elapsedCpuTime(), maxOp<double>())
               /interval_
            );
            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case ocClockTime:
        {
            label executionIndex = label
            (
                returnReduce(label(time_.elapsedClockTime()), maxOp<label>())
               /interval_
            );
            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case ocOnEnd:
        {
            scalar endTime = time_.endTime().value() - 0.5*time_.deltaTValue();
            return time_.value() > endTime;
            break;
        }

        case ocNone:
        {
            return false;
        }

        default:
        {
            FatalErrorInFunction
                << "Undefined time control: "
                << timeControlNames_[timeControl_] << nl
                << abort(FatalError);
            break;
        }
    }

    return false;
}


// ************************************************************************* //
