/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "decompositionConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionConstraint, 1);
    defineRunTimeSelectionTable(decompositionConstraint, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraint::decompositionConstraint
(
    const dictionary& constraintDict
)
:
    coeffDict_(constraintDict)
{}


Foam::decompositionConstraint::decompositionConstraint
(
    const dictionary& constraintDict,
    const word&
)
:
    coeffDict_(constraintDict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionConstraint>
Foam::decompositionConstraint::New
(
    const dictionary& dict
)
{
    return decompositionConstraint::New
    (
        dict,
        dict.get<word>("type")
    );
}


Foam::autoPtr<Foam::decompositionConstraint>
Foam::decompositionConstraint::New
(
    const dictionary& dict,
    const word& modelType
)
{
    Info<< "Selecting decompositionConstraint " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "decompositionConstraint",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<decompositionConstraint>(cstrIter()(dict));
}


// ************************************************************************* //
