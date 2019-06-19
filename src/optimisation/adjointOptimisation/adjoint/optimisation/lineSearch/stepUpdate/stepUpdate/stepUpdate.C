/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "stepUpdate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(stepUpdate, 0);
defineRunTimeSelectionTable(stepUpdate, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stepUpdate::stepUpdate(const dictionary& dict)
:
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

autoPtr<stepUpdate> stepUpdate::New(const dictionary& dict)
{
    const word modelType =
        dict.lookupOrDefault<word>("stepUpdateType", "bisection");

    Info<< "stepUpdate type : " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown stepUpdate type " << modelType
            << nl << nl
            << "Valid stepUpdate types are : " << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<stepUpdate>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void stepUpdate::setDeriv(const scalar deriv)
{
    // Does nothing in base
}


void stepUpdate::setNewMeritValue(const scalar value)
{
    // Does nothing in base
}


void stepUpdate::setOldMeritValue(const scalar value)
{
    // Does nothing in base
}


void stepUpdate::setInitialStep(const scalar value)
{
    // Does nothing in base
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
