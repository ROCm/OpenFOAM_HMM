/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "massRosinRammler.H"
#include "MathFunctions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(massRosinRammler, 0);
    addToRunTimeSelectionTable(distributionModel, massRosinRammler, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::massRosinRammler::massRosinRammler
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    lambda_(distributionModelDict_.getCompat<scalar>("lambda", {{"d", 2112}})),
    n_(distributionModelDict_.get<scalar>("n"))
{
    if (lambda_ < VSMALL || n_ < VSMALL)
    {
        FatalErrorInFunction
            << "Scale/Shape parameter cannot be equal to or less than zero:"
            << "    lambda = " << lambda_
            << "    n = " << n_
            << exit(FatalError);
    }

    check();
}


Foam::distributionModels::massRosinRammler::massRosinRammler
(
    const massRosinRammler& p
)
:
    distributionModel(p),
    lambda_(p.lambda_),
    n_(p.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::massRosinRammler::sample() const
{
    scalar d;

    // Re-sample if the calculated d is out of the physical range
    do
    {
        const scalar a = 3/n_ + 1;
        const scalar P = rndGen_.sample01<scalar>();
        const scalar x = Math::invIncGamma(a, P);
        d = lambda_*pow(x, 1/n_);
    } while (d < minValue_ || d > maxValue_);

    return d;
}


Foam::scalar Foam::distributionModels::massRosinRammler::meanValue() const
{
    return lambda_;
}


// ************************************************************************* //
