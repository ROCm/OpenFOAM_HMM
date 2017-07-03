/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModel> Foam::diameterModel::New
(
    const dictionary& dict,
    const phaseModel& phase
)
{
    const word modelType(dict.lookup("diameterModel"));

    Info << "Selecting diameterModel for phase "
        << phase.name()
        << ": "
        << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
           << "Unknown diameterModel type "
           << modelType << nl << nl
           << "Valid diameterModel types :" << endl
           << dictionaryConstructorTablePtr_->sortedToc()
           << exit(FatalError);
    }

    return cstrIter()
    (
        dict.optionalSubDict(modelType + "Coeffs"),
        phase
    );
}


// ************************************************************************* //
