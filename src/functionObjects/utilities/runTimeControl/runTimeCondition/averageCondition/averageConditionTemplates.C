/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::runTimeControls::averageCondition::calc
(
    const word& fieldName,
    const scalar alpha,
    const scalar beta,
    bool& satisfied,
    bool& processed
)
{
    const word valueType =
        state_.objectResultType(functionObjectName_, fieldName);

    if (pTraits<Type>::typeName != valueType)
    {
        return;
    }

    Type currentValue =
        state_.getObjectResult<Type>(functionObjectName_, fieldName);

    const word meanName(fieldName + "Mean");

    Type meanValue = state_.getResult<Type>(meanName);
    meanValue = alpha*meanValue + beta*currentValue;

    scalar delta = mag(meanValue - currentValue);

    if (log_)
    {
        Info<< "        " << meanName << ": " << meanValue
            << ", delta: " << delta << nl;
    }

    state_.setResult(meanName, meanValue);

    if (delta > tolerance_)
    {
        satisfied = false;
    }

    processed = true;
}


// ************************************************************************* //
