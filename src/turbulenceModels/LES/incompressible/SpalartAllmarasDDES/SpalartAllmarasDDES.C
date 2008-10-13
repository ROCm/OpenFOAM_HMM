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

#include "SpalartAllmarasDDES.H"
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

defineTypeNameAndDebug(SpalartAllmarasDDES, 0);
addToRunTimeSelectionTable(LESModel, SpalartAllmarasDDES, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmarasDDES::rd
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


tmp<volScalarField> SpalartAllmarasDDES::fd(const volScalarField& S)
{
    return 1.0 - tanh(pow3(8.0*rd(nuSgs_ + transport().nu(), S)));
}


void SpalartAllmarasDDES::dTildaUpdate(const volScalarField& S)
{
    dTilda_ =
        wallDist(mesh_).y()
      - fd(S)*max
        (
            dimensionedScalar("zero", dimLength, 0.0),
            wallDist(mesh_).y() - CDES_*delta()
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmarasDDES::SpalartAllmarasDDES
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    SpalartAllmaras(U, phi, transport, typeName)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
