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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable(distributionModel, exponential, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::exponential::exponential
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    lambda_(distributionModelDict_.get<scalar>("lambda"))
{
    if (lambda_ < VSMALL)
    {
        FatalErrorInFunction
            << "Rate parameter cannot be equal to or less than zero:" << nl
            << "    lambda = " << lambda_
            << exit(FatalError);
    }

    check();
}


Foam::distributionModels::exponential::exponential(const exponential& p)
:
    distributionModel(p),
    lambda_(p.lambda_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::exponential::~exponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::exponential::sample() const
{
    scalar y = rndGen_.sample01<scalar>();
    scalar K = exp(-lambda_*maxValue_) - exp(-lambda_*minValue_);
    return -(1.0/lambda_)*log(exp(-lambda_*minValue_) + y*K);
}


Foam::scalar Foam::distributionModels::exponential::meanValue() const
{
    return 1.0/lambda_;
}


// ************************************************************************* //
