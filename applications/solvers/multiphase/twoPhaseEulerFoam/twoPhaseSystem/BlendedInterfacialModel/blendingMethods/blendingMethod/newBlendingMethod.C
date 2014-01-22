/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::blendingMethod> Foam::blendingMethod::New
(
    const dictionary& dict,
    const phaseModel& phase1,
    const phaseModel& phase2
)
{
    word blendingMethodType(dict.lookup("type"));

    Info<< "Selecting blendingMethod for "
        << phase1.name() << " and " << phase2.name() << ": " 
        << blendingMethodType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(blendingMethodType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("blendingMethod::New")
            << "Unknown blendingMethodType type "
            << blendingMethodType << endl << endl
            << "Valid blendingMethod types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, phase1, phase2);
}


// ************************************************************************* //
