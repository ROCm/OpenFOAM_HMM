/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "normal.H"
#include "addToRunTimeSelectionTable.H"
#include "MathFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(normal, 0);
    addToRunTimeSelectionTable(distributionModel, normal, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::normal::normal
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    mu_
    (
        distributionModelDict_.getCompat<scalar>
        (
            "mu",
            {{"expectation", 2112}}
        )
    ),
    sigma_
    (
        distributionModelDict_.getCompat<scalar>
        (
            "sigma",
            {{"variance", 2112}}
        )
    )
{
    if (mag(sigma_) == 0)
    {
        FatalErrorInFunction
            << "Standard deviation cannot be zero." << nl
            << "    sigma = " << sigma_ << nl
            << exit(FatalError);
    }

    check();
}


Foam::distributionModels::normal::normal(const normal& p)
:
    distributionModel(p),
    mu_(p.mu_),
    sigma_(p.sigma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::normal::sample() const
{

    scalar a = erf((minValue_ - mu_)/sigma_);
    scalar b = erf((maxValue_ - mu_)/sigma_);

    scalar y = rndGen_.sample01<scalar>();
    scalar x = Math::erfInv(y*(b - a) + a)*sigma_ + mu_;

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    x = min(max(x, minValue_), maxValue_);

    return x;
}


Foam::scalar Foam::distributionModels::normal::meanValue() const
{
    return mu_;
}


// ************************************************************************* //
