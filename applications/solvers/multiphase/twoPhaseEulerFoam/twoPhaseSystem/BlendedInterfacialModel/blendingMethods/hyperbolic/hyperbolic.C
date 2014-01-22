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

#include "hyperbolic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(hyperbolic, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        hyperbolic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::hyperbolic::hyperbolic
(
    const dictionary& dict,
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    blendingMethod(dict, phase1, phase2),
    maxDispersedAlpha1_
    (
        "maxDispersedAlpha1",
        dimless,
        dict.lookup
        (
            IOobject::groupName("maxDispersedAlpha", phase1.name())
        )
    ),
    maxDispersedAlpha2_
    (
        "maxDispersedAlpha2",
        dimless,
        dict.lookup
        (
            IOobject::groupName("maxDispersedAlpha", phase2.name())
        )
    ),
    transitionAlphaScale_
    (
        "transitionAlphaScale",
        dimless,
        dict.lookup("transitionAlphaScale")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::hyperbolic::~hyperbolic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::f1() const
{
    return
        (
            1
          + tanh
            (
                (4/transitionAlphaScale_)
               *(phase1_ - maxDispersedAlpha1_)
            )
        )/2;
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::f2() const
{
    return
        (
            1
          + tanh
            (
                (4/transitionAlphaScale_)
               *(maxDispersedAlpha2_ - phase2_)
            )
        )/2;
}


// ************************************************************************* //
