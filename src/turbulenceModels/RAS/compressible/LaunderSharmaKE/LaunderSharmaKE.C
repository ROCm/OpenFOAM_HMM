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

#include "LaunderSharmaKE.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RAS
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaunderSharmaKE, 0);
addToRunTimeSelectionTable(turbulenceModel, LaunderSharmaKE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> LaunderSharmaKE::fMu() const
{
    return exp(-3.4/sqr(scalar(1) + rho_*sqr(k_)/(mu()*epsilon_)/50.0));
}


tmp<volScalarField> LaunderSharmaKE::f2() const
{
    return
        scalar(1)
      - 0.3*exp(-min(sqr(rho_*sqr(k_)/(mu()*epsilon_)), scalar(50.0)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
LaunderSharmaKE::LaunderSharmaKE
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
    C3(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C3", -0.33)),
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
        Cmu*fMu()*rho_*sqr(k_)/(epsilon_ + epsilonSmall_)
    )
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LaunderSharmaKE::R() const
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


tmp<volSymmTensorField> LaunderSharmaKE::devRhoReff() const
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


tmp<fvVectorMatrix> LaunderSharmaKE::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U) - fvc::div(muEff()*dev2(fvc::grad(U)().T()))
    );
}


bool LaunderSharmaKE::read()
{
    if (turbulenceModel::read())
    {
        turbulenceModelCoeffs_.readIfPresent<scalar>("Cmu", Cmu);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C1", C1);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C2", C2);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C3", C3);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphaEps", alphaEps);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphah", alphah);

        return true;
    }
    else
    {
        return false;
    }
}


void LaunderSharmaKE::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*Cmu*fMu()*sqr(k_)/(epsilon_ + epsilonSmall_);
        return;
    }

    turbulenceModel::correct();

    // Calculate parameters and coefficients for Launder-Sharma low-Reynolds
    // number model

    volScalarField E = 2.0*mu()*mut_*fvc::magSqrGradGrad(U_)/rho_;
    volScalarField D = 2.0*mu()*magSqr(fvc::grad(sqrt(k_)))/rho_;

    volScalarField divU = fvc::div(phi_/fvc::interpolate(rho_));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volScalarField G = mut_*(tgradU() && dev(twoSymm(tgradU())));
    tgradU.clear();


    // Dissipation equation

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1*G*epsilon_/k_ + fvm::SuSp((C3 - 2.0/3.0*C1)*rho_*divU, epsilon_)
      - fvm::Sp(C2*f2()*rho_*epsilon_/k_, epsilon_)
    //+ 0.75*1.5*flameKproduction*epsilon_/k_
      + E
    );

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
        G - fvm::SuSp(2.0/3.0*rho_*divU, k_)
      - fvm::Sp(rho_*(epsilon_ + D)/k_, k_)
    //+ flameKproduction
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    mut_ = Cmu*fMu()*rho_*sqr(k_)/epsilon_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RAS
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
