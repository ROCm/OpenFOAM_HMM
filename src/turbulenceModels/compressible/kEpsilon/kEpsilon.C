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

#include "kEpsilon.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace turbulenceModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsilon, 0);
addToRunTimeSelectionTable(turbulenceModel, kEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilon::kEpsilon
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    basicThermo& thermophysicalModel
)
:
    turbulenceModel(typeName, rho, U, phi, thermophysicalModel),

    Cmu(turbulenceModelCoeffs_.lookupOrDefault<scalar>("Cmu", 0.09)),
    C1(turbulenceModelCoeffs_.lookupOrDefault<scalar>("C1", 1.44)),
    C2(turbulenceModelCoeffs_.lookupOrDefault<scalar>("C2", 1.92)),
    C3(turbulenceModelCoeffs_.lookupOrDefault<scalar>("C3", -0.33)),
    alphak(turbulenceModelCoeffs_.lookupOrDefault<scalar>("alphak", 1.0)),
    alphaEps
    (
        turbulenceModelCoeffs_.lookupOrDefault<scalar>("alphaEps", 0.76923)
    ),
    alphah(turbulenceModelCoeffs_.lookupOrDefault<scalar>("alphah", 1.0)),

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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kEpsilon::R() const
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


tmp<volSymmTensorField> kEpsilon::devRhoReff() const
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


tmp<fvVectorMatrix> kEpsilon::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U)
      - fvc::div(muEff()*dev2(fvc::grad(U)().T()))
    );
}


bool kEpsilon::read()
{
    if (turbulenceModel::read())
    {
        Cmu = turbulenceModelCoeffs_.lookupOrDefault<scalar>("Cmu", 0.09);
        C1 = turbulenceModelCoeffs_.lookupOrDefault<scalar>("C1", 1.44);
        C2 = turbulenceModelCoeffs_.lookupOrDefault<scalar>("C2", 1.92);
        C3 = turbulenceModelCoeffs_.lookupOrDefault<scalar>("C3", -0.33);
        alphak = turbulenceModelCoeffs_.lookupOrDefault<scalar>("alphak", 1.0);
        alphaEps = turbulenceModelCoeffs_.lookupOrDefault<scalar>
            (
                "alphaEps",
                0.76923
            );
        alphah = turbulenceModelCoeffs_.lookupOrDefault<scalar>("alphah", 1.0);

        return true;
    }
    else
    {
        return false;
    }
}


void kEpsilon::correct()
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
    volScalarField G = mut_*(tgradU() && dev(twoSymm(tgradU())));
    tgradU.clear();

#   include "wallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1*G*epsilon_/k_
      - fvm::SuSp(((2.0/3.0)*C1 + C3)*rho_*divU, epsilon_)
      - fvm::Sp(C2*rho_*epsilon_/k_, epsilon_)
    );

#   include "wallDissipationI.H"

    epsEqn().relax();

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::SuSp((2.0/3.0)*rho_*divU, k_)
      - fvm::Sp(rho_*epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    mut_ = rho_*Cmu*sqr(k_)/epsilon_;

#   include "wallViscosityI.H"

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
