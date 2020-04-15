/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "Time.H"
#include "FIFOStack.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::runTimeControls::averageCondition::calc
(
    const label fieldi,
    bool& satisfied,
    bool& processed
)
{
    const word& fieldName = fieldNames_[fieldi];

    const word valueType =
        state_.objectResultType(functionObjectName_, fieldName);

    if (pTraits<Type>::typeName != valueType)
    {
        return;
    }

    const scalar dt = obr_.time().deltaTValue();

    const Type currentValue =
        state_.getObjectResult<Type>(functionObjectName_, fieldName);
    const word meanName(fieldName + "Mean");

    // Current mean value
    Type meanValue = state_.getResult<Type>(meanName);

    switch (windowType_)
    {
        case windowType::NONE:
        {
            const scalar Dt = totalTime_[fieldi];
            const scalar alpha = (Dt - dt)/Dt;
            const scalar beta = dt/Dt;
            meanValue = alpha*meanValue + beta*currentValue;

            break;
        }
        case windowType::APPROXIMATE:
        {
            const scalar Dt = totalTime_[fieldi];
            scalar alpha = (Dt - dt)/Dt;
            scalar beta = dt/Dt;
            if (Dt - dt >= window_)
            {
                alpha = (window_ - dt)/window_;
                beta = dt/window_;
            }
            else
            {
                // Ensure that averaging is performed over window time
                // before condition can be satisfied
                satisfied = false;
            }

            meanValue = alpha*meanValue + beta*currentValue;
            totalTime_[fieldi] += dt;
            break;
        }
        case windowType::EXACT:
        {
            FIFOStack<scalar> windowTimes;
            FIFOStack<Type> windowValues;
            dictionary& dict = this->conditionDict().subDict(fieldName);
            dict.readIfPresent("windowTimes", windowTimes);
            dict.readIfPresent("windowValues", windowValues);

            // Increment time for all existing values
            for (scalar& dti : windowTimes)
            {
                dti += dt;
            }

            // Remove any values outside the window
            bool removeValue = true;
            while (removeValue && windowTimes.size())
            {
                removeValue = windowTimes.first() > window_;

                if (removeValue)
                {
                    windowTimes.pop();
                    windowValues.pop();
                }
            }

            // Add the current value
            windowTimes.push(dt);
            windowValues.push(currentValue);

            // Calculate the window average
            auto timeIter = windowTimes.cbegin();
            auto valueIter = windowValues.cbegin();

            meanValue = pTraits<Type>::zero;
            Type valueOld(pTraits<Type>::zero);

            for
            (
                label i = 0;
                timeIter.good();
                ++i, ++timeIter, ++valueIter
            )
            {
                const Type& value = valueIter();
                const scalar dt = timeIter();

                meanValue += dt*value;

                if (i)
                {
                    meanValue -= dt*valueOld;
                }

                valueOld = value;
            }

            meanValue /= windowTimes.first();

            // Store the state information for the next step
            dict.set("windowTimes", windowTimes);
            dict.set("windowValues", windowValues);

            break;
        }
    }

    scalar delta = mag(meanValue - currentValue);

    Log << "        " << meanName << ": " << meanValue
        << ", delta: " << delta << nl;

    // Note: Writing result to owner function object and not the local run-time
    // condition
    state_.setResult(meanName, meanValue);

    if (delta > tolerance_)
    {
        satisfied = false;
    }

    processed = true;
}


// ************************************************************************* //
