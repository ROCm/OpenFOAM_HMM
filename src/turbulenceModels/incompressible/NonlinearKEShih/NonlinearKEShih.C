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

#include "NonlinearKEShih.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NonlinearKEShih, 0);
addToRunTimeSelectionTable(turbulenceModel, NonlinearKEShih, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
NonlinearKEShih::NonlinearKEShih
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    turbulenceModel(typeName, U, phi, lamTransportModel),

    C1(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C1", 1.44)),
    C2(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C2", 1.92)),
    alphak(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphak", 1.0)),
    alphaEps
    (
        turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphaEps", 0.76923)
    ),
    A1(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("A1", 1.25)),
    A2(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("A2", 1000.0)),
    Ctau1(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Ctau1", -4.0)),
    Ctau2(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Ctau2", 13.0)),
    Ctau3(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Ctau3", -2.0)),
    alphaKsi
    (
        turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphaKsi", 0.9)
    ),

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

    gradU(fvc::grad(U)),
    eta(k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU + gradU.T())))),
    ksi(k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU - gradU.T())))),
    Cmu(2.0/(3.0*(A1 + eta + alphaKsi*ksi))),
    fEta(A2 + pow(eta, 3.0)),

    nut_(Cmu*sqr(k_)/(epsilon_ + epsilonSmall_)),

    nonlinearStress
    (
        "nonlinearStress",
        symm
        (
        pow(k_, 3.0)/sqr(epsilon_)*
        (
            Ctau1/fEta*
            (
                (gradU & gradU)
              + (gradU & gradU)().T()
            )
          + Ctau2/fEta*(gradU & gradU.T())
          + Ctau3/fEta*(gradU.T() & gradU)
        )
        )
    )
{
#   include "wallNonlinearViscosityI.H"

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> NonlinearKEShih::R() const
{
    volTensorField gradU = fvc::grad(U_);

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
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(gradU) + nonlinearStress,
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> NonlinearKEShih::devReff() const
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
           -nuEff()*dev(twoSymm(fvc::grad(U_))) + nonlinearStress
        )
    );
}


tmp<fvVectorMatrix> NonlinearKEShih::divDevReff(volVectorField& U) const
{
    return
    (
        fvc::div(nonlinearStress)
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool NonlinearKEShih::read()
{
    if (turbulenceModel::read())
    {
        turbulenceModelCoeffs_.readIfPresent<scalar>("C1", C1);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C2", C2);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphak", alphak);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphaEps", alphaEps);
        turbulenceModelCoeffs_.readIfPresent<scalar>("A1", A1);
        turbulenceModelCoeffs_.readIfPresent<scalar>("A2", A2);
        turbulenceModelCoeffs_.readIfPresent<scalar>("Ctau1", Ctau1);
        turbulenceModelCoeffs_.readIfPresent<scalar>("Ctau2", Ctau2);
        turbulenceModelCoeffs_.readIfPresent<scalar>("Ctau3", Ctau3);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphaKsi", alphaKsi);

        return true;
    }
    else
    {
        return false;
    }
}


void NonlinearKEShih::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    turbulenceModel::correct();

    gradU = fvc::grad(U_);

    // generation term
    volScalarField S2 = symm(gradU) && gradU;

    volScalarField G =
        Cmu*sqr(k_)/epsilon_*S2
      - (nonlinearStress && gradU);

#   include "nonLinearWallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1*G*epsilon_/k_
      - fvm::Sp(C2*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "wallDissipationI.H"

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity

    eta = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU + gradU.T())));
    ksi = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU - gradU.T())));
    Cmu = 2.0/(3.0*(A1 + eta + alphaKsi*ksi));
    fEta = A2 + pow(eta, 3.0);

    nut_ = Cmu*sqr(k_)/epsilon_;

#   include "wallNonlinearViscosityI.H"

    nonlinearStress = symm
    (
        pow(k_, 3.0)/sqr(epsilon_)*
        (
            Ctau1/fEta*
            (
                (gradU & gradU)
              + (gradU & gradU)().T()
            )
          + Ctau2/fEta*(gradU & gradU.T())
          + Ctau3/fEta*(gradU.T() & gradU)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceModels
} // End namespace Foam

// ************************************************************************* //
