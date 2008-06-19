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

#include "spectEddyVisc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LES
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(spectEddyVisc, 0);
addToRunTimeSelectionTable(LESmodel, spectEddyVisc, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
spectEddyVisc::spectEddyVisc
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESmodel(typeName, U, phi, transport),
    GenEddyVisc(U, phi, transport),

    cB_(LESmodelProperties().lookupOrAddDefault<scalar>("cB", 8.22)),
    cK1_(LESmodelProperties().lookupOrAddDefault<scalar>("cK1", 0.83)),
    cK2_(LESmodelProperties().lookupOrAddDefault<scalar>("cK2", 1.03)),
    cK3_(LESmodelProperties().lookupOrAddDefault<scalar>("cK3", 4.75)),
    cK4_(LESmodelProperties().lookupOrAddDefault<scalar>("cK4", 2.55))
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> spectEddyVisc::k() const
{
    volScalarField Eps = 2*nuEff()*magSqr(symm(fvc::grad(U())));

    return
        cK1_*pow(delta(), 2.0/3.0)*pow(Eps, 2.0/3.0)
        *exp(-cK2_*pow(delta(), -4.0/3.0)*nu()/pow(Eps, 1.0/3.0))
      - cK3_*pow(Eps*nu(), 1.0/2.0)
       *erfc(cK4_*pow(delta(), -2.0/3.0)*pow(Eps, -1.0/6.0));
}


void spectEddyVisc::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);

    volScalarField Re = sqr(delta())*mag(symm(gradU))/nu();

    for (label i=0; i<5; i++)
    {
        nuSgs_ =
            nu()
           /(
               scalar(1)
               - exp(-cB_*pow(nu()/(nuSgs_ + nu()), 1.0/3.0)*pow(Re, -2.0/3.0))
            );
    }

    nuSgs_.correctBoundaryConditions();
}


bool spectEddyVisc::read()
{
    if (GenEddyVisc::read())
    {
        LESmodelProperties().readIfPresent<scalar>("cB", cB_);
        LESmodelProperties().readIfPresent<scalar>("cK1", cK1_);
        LESmodelProperties().readIfPresent<scalar>("cK2", cK2_);
        LESmodelProperties().readIfPresent<scalar>("cK3", cK3_);
        LESmodelProperties().readIfPresent<scalar>("cK4", cK4_);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LES
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
