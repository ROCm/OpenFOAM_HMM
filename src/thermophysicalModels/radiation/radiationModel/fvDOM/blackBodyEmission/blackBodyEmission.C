/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "blackBodyEmission.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::blackBodyEmission::blackBodyEmission
(
    const fileName& name,
    const word& instance,
    const label nLambda,
    const volScalarField& T
)
:
    blackBodyEmissiveTable_(name, instance, T.mesh()),
    C1_("C1",dimensionSet(1, 4, 3, 0, 0, 0, 0), 3.7419e-16),
    C2_("C2",dimensionSet(0, 1, 0, 1, 0, 0, 0), 14.388e-6),
    bLambda_(nLambda),
    T_(T)
{
    forAll(bLambda_, lambdaI)
    {
        bLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "bLambda_" + Foam::name(lambdaI) ,
                    T.mesh().time().timeName(),
                    T.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                radiation::sigmaSB*pow4(T)
            )
        );

    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::blackBodyEmission::~blackBodyEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::blackBodyEmission::fLambdaT
(
    const scalar lambdaT
) const
{
    return  blackBodyEmissiveTable_.lookUp(lambdaT*1.0e6)[1];
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::blackBodyEmission::EbDeltaLambdaT
(
    const volScalarField& T,
    const Vector2D<scalar>& band
) const
{
    tmp<volScalarField> Eb
    (
        new volScalarField
        (
            IOobject
            (
                "Eb",
                T.mesh().time().timeName(),
                T.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            radiation::sigmaSB*pow4(T)
        )
    );


    if (band == Vector2D<scalar>::one)
    {
        return Eb;
    }
    else
    {
        forAll(T, i)
        {
            scalar T1 = fLambdaT(band[1]*T[i]);
            scalar T2 = fLambdaT(band[0]*T[i]);
            dimensionedScalar fLambdaDelta
            (
                "fLambdaDelta",
                dimless,
                T1 - T2
            );
            Eb()[i] = Eb()[i]*fLambdaDelta.value();
        }
        return Eb;
    }
}


void Foam::radiation::blackBodyEmission::correct
(
    const label lambdaI,
    const Vector2D<scalar>& band
)
{
    bLambda_[lambdaI] = EbDeltaLambdaT(T_, band);
}


// ************************************************************************* //
