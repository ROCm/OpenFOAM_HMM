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

#include "RNGkEpsilon.H"
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

defineTypeNameAndDebug(RNGkEpsilon, 0);
addToRunTimeSelectionTable(turbulenceModel, RNGkEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
RNGkEpsilon::RNGkEpsilon
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    basicThermo& thermophysicalModel
)
:
    turbulenceModel(typeName, rho, U, phi, thermophysicalModel),

    Cmu(turbulenceModelCoeffs_.lookup("Cmu")),
    C1(turbulenceModelCoeffs_.lookup("C1")),
    C2(turbulenceModelCoeffs_.lookup("C2")),
    C3(turbulenceModelCoeffs_.lookup("C3")),
    alphak(turbulenceModelCoeffs_.lookup("alphak")),
    alphaEps(turbulenceModelCoeffs_.lookup("alphaEps")),
    alphah(turbulenceModelCoeffs_.lookup("alphah")),
    eta0(turbulenceModelCoeffs_.lookup("eta0")),
    beta(turbulenceModelCoeffs_.lookup("beta")),

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

tmp<volSymmTensorField> RNGkEpsilon::R() const
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


tmp<volSymmTensorField> RNGkEpsilon::devRhoReff() const
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


tmp<fvVectorMatrix> RNGkEpsilon::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U) - fvc::div(muEff()*dev2(fvc::grad(U)().T()))
    );
}


bool RNGkEpsilon::read()
{
    if (turbulenceModel::read())
    {
        turbulenceModelCoeffs_.lookup("Cmu") >> Cmu;
        turbulenceModelCoeffs_.lookup("C1") >> C1;
        turbulenceModelCoeffs_.lookup("C2") >> C2;
        turbulenceModelCoeffs_.lookup("C3") >> C3;
        turbulenceModelCoeffs_.lookup("alphak") >> alphak;
        turbulenceModelCoeffs_.lookup("alphaEps") >> alphaEps;
        turbulenceModelCoeffs_.lookup("alphah") >> alphah;
        turbulenceModelCoeffs_.lookup("eta0") >> eta0;
        turbulenceModelCoeffs_.lookup("beta") >> beta;

        return true;
    }
    else
    {
        return false;
    }
}


void RNGkEpsilon::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*Cmu*sqr(k_)/(epsilon_ + epsilonSmall_);
        return;
    }

    turbulenceModel::correct();

    volScalarField divU = fvc::div(phi_/fvc::interpolate(rho_));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volScalarField S2 = (tgradU() && dev(twoSymm(tgradU())));
    tgradU.clear();

    volScalarField G = mut_*S2;

    volScalarField eta = sqrt(mag(S2))*k_/epsilon_;
    volScalarField eta3 = eta*sqr(eta);

    volScalarField R = 
        ((eta*(-eta/eta0 + scalar(1)))/(beta*eta3 + scalar(1)));

#   include "wallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        (C1 - R)*G*epsilon_/k_
      - fvm::SuSp(((2.0/3.0)*C1 + C3)*rho_*divU, epsilon_)
      - fvm::Sp(C2*rho_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "wallDissipationI.H"

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
      - fvm::Sp(rho_*(epsilon_)/k_, k_)
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
