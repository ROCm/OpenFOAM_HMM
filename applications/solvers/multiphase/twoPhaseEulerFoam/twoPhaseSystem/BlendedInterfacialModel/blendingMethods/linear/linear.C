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

#include "linear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(linear, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        linear,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::linear
(
    const dictionary& dict,
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    blendingMethod(dict, phase1, phase2),
    maxFullyDispersedAlpha1_
    (
        "maxFullyDispersedAlpha1",
        dimless,
        dict.lookup
        (
            IOobject::groupName("maxFullyDispersedAlpha", phase1.name())
        )
    ),
    maxPartlyDispersedAlpha1_
    (
        "maxPartlyDispersedAlpha1",
        dimless,
        dict.lookup
        (
            IOobject::groupName("maxPartlyDispersedAlpha", phase1.name())
        )
    ),
    maxFullyDispersedAlpha2_
    (
        "maxFullyDispersedAlpha2",
        dimless,
        dict.lookup
        (
            IOobject::groupName("maxFullyDispersedAlpha", phase2.name())
        )
    ),
    maxPartlyDispersedAlpha2_
    (
        "maxPartlyDispersedAlpha2",
        dimless,
        dict.lookup
        (
            IOobject::groupName("maxPartlyDispersedAlpha", phase2.name())
        )
    )
{
    if
    (
        maxFullyDispersedAlpha1_ > maxPartlyDispersedAlpha1_
     || maxFullyDispersedAlpha2_ > maxPartlyDispersedAlpha2_
    )
    {
        FatalErrorIn
        (
            "Foam::blendingMethods::linear::linear"
            "("
                "const dictionary& dict,"
                "const phaseModel& phase1,"
                "const phaseModel& phase2"
            ")"
        )   << "The supplied fully dispersed volume fraction is greater than "
            << "the partly dispersed value"
            << endl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::~linear()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linear::f1() const
{
    return
        min
        (
            max
            (
                (phase1_ - maxFullyDispersedAlpha1_)
               /(maxPartlyDispersedAlpha1_ - maxFullyDispersedAlpha1_ + SMALL),
                0.0
            ),
            1.0
        );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linear::f2() const
{
    return
        min
        (
            max
            (
                (maxPartlyDispersedAlpha2_ - phase2_)
               /(maxPartlyDispersedAlpha2_ - maxFullyDispersedAlpha2_ + SMALL),
                0.0
            ),
            1.0
        );
}


// ************************************************************************* //
