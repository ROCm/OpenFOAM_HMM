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

#include "PDRkEpsilon.H"
#include "wallFvPatch.H"
#include "PDRDragModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RAS
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PDRkEpsilon, 0);
addToRunTimeSelectionTable(turbulenceModel, PDRkEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PDRkEpsilon::PDRkEpsilon
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    basicThermo& thermophysicalModel
)
:
    turbulenceModel(typeName, rho, U, phi, thermophysicalModel),

    Cmu(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Cmu", 0.09)),
    C1(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C1", 1.44)),
    C2(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C2", 1.92)),
    alphak(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphak", 1.0)),
    alphaEps
    (
        turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphaEps", 0.76923)
    ),
    alphah(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphah", 1.0)),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Cmu*rho_*sqr(k_)/(epsilon_ + epsilonSmall_)
    )
{
#   include "wallViscosityI.H"

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> PDRkEpsilon::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - (mut_/rho_)*dev(twoSymm(fvc::grad(U_))),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> PDRkEpsilon::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -muEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> PDRkEpsilon::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U) - fvc::div(muEff()*dev2(fvc::grad(U)().T()))
    );
}


bool PDRkEpsilon::read()
{
    if (turbulenceModel::read())
    {
        turbulenceModelCoeffs_.readIfPresent<scalar>("Cmu", Cmu);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C1", C1);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C2", C2);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphak", alphak);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphaEps", alphaEps);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphah", alphah);

        return true;
    }
    else
    {
        return false;
    }
}


void PDRkEpsilon::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*Cmu*sqr(k_)/(epsilon_ + epsilonSmall_);
#       include "wallViscosityI.H"
        return;
    }

    turbulenceModel::correct();

    volScalarField divU = fvc::div(phi_/fvc::interpolate(rho_));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volScalarField G = 2*mut_*(tgradU() && dev(symm(tgradU())));
    tgradU.clear();

    // Add the blockage generation term so that it is included consistently
    // in both the k and epsilon equations
    const volScalarField& betav = U_.db().lookupObject<volScalarField>("betav");

    const PDRDragModel& drag =
        U_.db().lookupObject<PDRDragModel>("PDRDragModel");

    volScalarField GR = drag.Gk();

#   include "wallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        betav*fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1*(betav*G + GR)*epsilon_/k_
      - fvm::SuSp(((2.0/3.0)*C1)*betav*rho_*divU, epsilon_)
      - fvm::Sp(C2*betav*rho_*epsilon_/k_, epsilon_)
    );

#   include "wallDissipationI.H"

    epsEqn().relax();

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        betav*fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        betav*G + GR
      - fvm::SuSp((2.0/3.0)*betav*rho_*divU, k_)
      - fvm::Sp(betav*rho_*epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    mut_ = rho_*Cmu*sqr(k_)/epsilon_;

#   include "wallViscosityI.H"

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RAS
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
