/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "constrainHbyA.H"
#include "volFields.H"
#include "fixedFluxExtrapolatedPressureFvPatchScalarField.H"


#ifdef USE_ROCTX
#include <roctracer/roctx.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::constrainHbyA
(
    const tmp<volVectorField>& tHbyA,
    const volVectorField& U,
    const volScalarField& p
)
{

    #ifdef USE_ROCTX
    roctxRangePush("Foam::constrainHbyA");
    #endif

    tmp<volVectorField> tHbyANew;

    if (tHbyA.isTmp())
    {
        tHbyANew = tHbyA;
        tHbyANew.ref().rename("HbyA");
    }
    else
    {
        tHbyANew = new volVectorField("HbyA", tHbyA);
    }

    volVectorField& HbyA = tHbyANew.ref();
    volVectorField::Boundary& HbyAbf = HbyA.boundaryFieldRef();


    #ifdef USE_ROCTX
    roctxRangePush("Foam::constrainHbyA_loop");
    #endif

    forAll(U.boundaryField(), patchi)
    {
        if
        (
           !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            HbyAbf[patchi] = U.boundaryField()[patchi];
        }
    }

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return tHbyANew;
}


// ************************************************************************* //
