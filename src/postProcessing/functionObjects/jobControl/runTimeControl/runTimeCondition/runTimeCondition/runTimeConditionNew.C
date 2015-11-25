/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "runTimeCondition.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::runTimeCondition> Foam::runTimeCondition::New
(
    const word& conditionName,
    const objectRegistry& obr,
    const dictionary& dict,
    functionObjectState& state
)
{
    word conditionType(dict.lookup("type"));

    Info<< "Selecting runTimeCondition " << conditionType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(conditionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "runTimeCondition::New"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "functionObjectState&"
            ")"
        )   << "Unknown runTimeCondition type "
            << conditionType << nl << nl
            << "Valid runTimeCondition types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<runTimeCondition>
        (
            cstrIter()(conditionName, obr, dict, state)
        );
}


// ************************************************************************* //
