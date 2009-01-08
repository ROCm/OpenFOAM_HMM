/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "SpalartAllmarasIDDES.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmarasIDDES, 0);
addToRunTimeSelectionTable(LESModel, SpalartAllmarasIDDES, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmarasIDDES::alpha() const
{
    return
        0.25
      - wallDist(mesh_).y()
       /dimensionedScalar("hMax", dimLength, max(cmptMax(delta())));
}


tmp<volScalarField> SpalartAllmarasIDDES::ft(const volScalarField& S) const
{
    return tanh(pow3(sqr(ct_)*r(nuSgs_, S)));
}


tmp<volScalarField> SpalartAllmarasIDDES::fl(const volScalarField& S) const
{
    return tanh(pow(sqr(cl_)*r(transport().nu(), S), 10));
}


tmp<volScalarField> SpalartAllmarasIDDES::rd
(
    const volScalarField& visc,
    const volScalarField& S
) const
{
    volScalarField d = wallDist(mesh_).y();

    tmp<volScalarField> trd
    (
        new volScalarField
        (
            min
            (
                visc
               /(
                   max
                   (
                       S,
                       dimensionedScalar("SMALL", S.dimensions(), SMALL)
                   )*sqr(kappa_*d)
                 + dimensionedScalar
                   (
                       "ROOTVSMALL",
                       dimensionSet(0, 2 , -1, 0, 0),
                       ROOTVSMALL
                   )
               ), scalar(10.0)
            )
        )
    );

    return trd;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmarasIDDES::fd(const volScalarField& S) const
{
    return 1.0 - tanh(pow3(8.0*rd(nuSgs_+transport().nu(), S)));
}


void SpalartAllmarasIDDES::dTildaUpdate(const volScalarField& S)
{
    volScalarField alpha = this->alpha();

    volScalarField expTerm = exp(sqr(alpha));

    volScalarField fHill =
        2.0*(pos(alpha)*pow(expTerm, -11.09) + neg(alpha)*pow(expTerm, -9.0));


    volScalarField fStep = min(2.0*pow(expTerm, -9.0), scalar(1));
    volScalarField fHyb = max(1.0 - fd(S), fStep);

    volScalarField fAmp = 1.0 - max(ft(S), fl(S));

    volScalarField fRestore = max(fHill - 1.0, scalar(0))*fAmp;

    // volScalarField ft2 = IGNORING ft2 terms

    volScalarField Psi = sqrt
    (
        min
        (
            scalar(100),
            (1.0 - Cb1_/(Cw1_*sqr(kappa_)*fwStar_)*fv2())/max(SMALL, fv1())
        )
    );

    dTilda_ = max
    (
        dimensionedScalar("zero", dimLength, 0.0),
        fHyb*(1.0 + fRestore*Psi)*wallDist(mesh_).y()
      + (1.0 - fHyb)*CDES_*Psi*delta()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmarasIDDES::SpalartAllmarasIDDES
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    SpalartAllmaras(U, phi, transport, typeName),

    fwStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fwStar",
            coeffDict(),
            0.424
        )
    ),
    cl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cl",
            coeffDict(),
            3.55
        )
    ),
    ct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ct",
            coeffDict(),
            1.63
        )
    )

{}


bool SpalartAllmarasIDDES::read()
{
    if (SpalartAllmaras::read())
    {
        fwStar_.readIfPresent(coeffDict());
        cl_.readIfPresent(coeffDict());
        ct_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
