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
    minValue_(distributionModelDict_.get<scalar>("minValue")),
    maxValue_(distributionModelDict_.get<scalar>("maxValue")),
    expectation_(distributionModelDict_.get<scalar>("expectation")),
    variance_(distributionModelDict_.get<scalar>("variance")),
    a_(0.147)
{
    if (minValue_ < 0)
    {
        FatalErrorInFunction
            << "Minimum value must be greater than zero. "
            << "Supplied minValue = " << minValue_
            << abort(FatalError);
    }

    if (maxValue_ < minValue_)
    {
        FatalErrorInFunction
            << "Maximum value is smaller than the minimum value:"
            << "    maxValue = " << maxValue_ << ", minValue = " << minValue_
            << abort(FatalError);
    }
}


Foam::distributionModels::normal::normal(const normal& p)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    expectation_(p.expectation_),
    variance_(p.variance_),
    a_(p.a_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::normal::~normal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::normal::sample() const
{

    scalar a = erf((minValue_ - expectation_)/variance_);
    scalar b = erf((maxValue_ - expectation_)/variance_);

    scalar y = rndGen_.sample01<scalar>();
    scalar x = Math::erfInv(y*(b - a) + a)*variance_ + expectation_;

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    x = min(max(x, minValue_), maxValue_);

    return x;
}


Foam::scalar Foam::distributionModels::normal::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::normal::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::normal::meanValue() const
{
    return expectation_;
}


// ************************************************************************* //
