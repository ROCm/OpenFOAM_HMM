/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "interpolationWeights.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(interpolationWeights, 0);
defineRunTimeSelectionTable(interpolationWeights, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolationWeights::interpolationWeights(const scalarField& samples)
:
    samples_(samples)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interpolationWeights> Foam::interpolationWeights::New
(
    const word& type,
    const scalarField& samples
)
{
    if (debug)
    {
        InfoInFunction
            << "Selecting interpolationWeights "
            << type << endl;
    }

    auto cstrIter = wordConstructorTablePtr_->cfind(type);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown interpolationWeights type "
            << type << nl << nl
            << "Valid interpolationWeights types :" << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<interpolationWeights>(cstrIter()(samples));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interpolationWeights::~interpolationWeights()
{}


// ************************************************************************* //
