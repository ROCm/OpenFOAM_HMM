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

#include "LienCubicKELowRe.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RAS
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LienCubicKELowRe, 0);
addToRunTimeSelectionTable(RASmodel, LienCubicKELowRe, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
LienCubicKELowRe::LienCubicKELowRe
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASmodel(typeName, U, phi, lamTransportModel),

    C1(RASmodelCoeffs_.lookupOrAddDefault<scalar>("C1", 1.44)),
    C2(RASmodelCoeffs_.lookupOrAddDefault<scalar>("C2", 1.92)),
    alphak(RASmodelCoeffs_.lookupOrAddDefault<scalar>("alphak", 1.0)),
    alphaEps
    (
        RASmodelCoeffs_.lookupOrAddDefault<scalar>("alphaEps", 0.76923)
    ),
    A1(RASmodelCoeffs_.lookupOrAddDefault<scalar>("A1", 1.25)),
    A2(RASmodelCoeffs_.lookupOrAddDefault<scalar>("A2", 1000.0)),
    Ctau1(RASmodelCoeffs_.lookupOrAddDefault<scalar>("Ctau1", -4.0)),
    Ctau2(RASmodelCoeffs_.lookupOrAddDefault<scalar>("Ctau2", 13.0)),
    Ctau3(RASmodelCoeffs_.lookupOrAddDefault<scalar>("Ctau3", -2.0)),
    alphaKsi(RASmodelCoeffs_.lookupOrAddDefault<scalar>("alphaKsi", 0.9)),
    CmuWall(RASmodelCoeffs_.lookupOrAddDefault<scalar>("Cmu", 0.09)),
    Am(RASmodelCoeffs_.lookupOrAddDefault<scalar>("Am", 0.016)),
    Aepsilon
    (
        RASmodelCoeffs_.lookupOrAddDefault<scalar>("Aepsilon", 0.263)
    ),
    Amu(RASmodelCoeffs_.lookupOrAddDefault<scalar>("Amu", 0.00222)),

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

    y_(mesh_),

    gradU(fvc::grad(U)),
    eta(k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU + gradU.T())))),
    ksi(k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU - gradU.T())))),
    Cmu(2.0/(3.0*(A1 + eta + alphaKsi*ksi))),
    fEta(A2 + pow(eta, 3.0)),

    C5viscosity
    (
        -2.0*pow(Cmu, 3.0)*pow(k_, 4.0)/pow(epsilon_, 3.0)*
        (magSqr(gradU + gradU.T()) - magSqr(gradU - gradU.T()))
    ),

    yStar(sqrt(k_)*y_/nu() + SMALL),

    nut_
    (
        Cmu*
        (
            scalar(1) - exp(-Am*yStar))
           /(scalar(1) - exp(-Aepsilon*yStar) + SMALL
        )
        *sqr(k_)/(epsilon_ + epsilonSmall_)
        // cubic term C5, implicit part
      + max
        (
            C5viscosity,
            dimensionedScalar("0", C5viscosity.dimensions(), 0.0)
        )
    ),
    // turbulent viscosity, with implicit part of C5

    nonlinearStress
    (
        "nonlinearStress",
        symm
        (
        // quadratic terms
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
        // cubic term C4
      - 20.0*pow(k_, 4.0)/pow(epsilon_, 3.0)*
        pow(Cmu, 3.0)*
        (
            ((gradU & gradU) & gradU.T())
          + ((gradU & gradU.T()) & gradU.T())
          - ((gradU.T() & gradU) & gradU)
          - ((gradU.T() & gradU.T()) & gradU)
        )
        // cubic term C5, explicit part
      + min
        (
            C5viscosity,
            dimensionedScalar("0", C5viscosity.dimensions(), 0.0)
        )*gradU
        )
    )
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LienCubicKELowRe::R() const
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
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(gradU) + nonlinearStress,
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> LienCubicKELowRe::devReff() const
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


tmp<fvVectorMatrix> LienCubicKELowRe::divDevReff(volVectorField& U) const
{
    return
    (
        fvc::div(nonlinearStress)
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool LienCubicKELowRe::read()
{
    if (RASmodel::read())
    {
        RASmodelCoeffs_.readIfPresent<scalar>("C1", C1);
        RASmodelCoeffs_.readIfPresent<scalar>("C2", C2);
        RASmodelCoeffs_.readIfPresent<scalar>("alphak", alphak);
        RASmodelCoeffs_.readIfPresent<scalar>("alphaEps", alphaEps);
        RASmodelCoeffs_.readIfPresent<scalar>("A1", A1);
        RASmodelCoeffs_.readIfPresent<scalar>("A2", A2);
        RASmodelCoeffs_.readIfPresent<scalar>("Ctau1", Ctau1);
        RASmodelCoeffs_.readIfPresent<scalar>("Ctau2", Ctau2);
        RASmodelCoeffs_.readIfPresent<scalar>("Ctau3", Ctau3);
        RASmodelCoeffs_.readIfPresent<scalar>("alphaKsi", alphaKsi);
        RASmodelCoeffs_.readIfPresent<scalar>("Cmu", CmuWall);
        RASmodelCoeffs_.readIfPresent<scalar>("Am", Am);
        RASmodelCoeffs_.readIfPresent<scalar>("Aepsilon", Aepsilon);
        RASmodelCoeffs_.readIfPresent<scalar>("Amu", Amu);

        return true;
    }
    else
    {
        return false;
    }
}


void LienCubicKELowRe::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    RASmodel::correct();

    if (mesh_.changing())
    {
        y_.correct();
    }

    gradU = fvc::grad(U_);

    // generation term
    volScalarField S2 = symm(gradU) && gradU;

    yStar = sqrt(k_)*y_/nu() + SMALL;
    volScalarField Rt = sqr(k_)/(nu()*epsilon_);

    volScalarField fMu =
        (scalar(1) - exp(-Am*yStar))
       /(scalar(1) - exp(-Aepsilon*yStar) + SMALL);

    volScalarField f2 = scalar(1) - 0.3*exp(-sqr(Rt));

    volScalarField G =
        Cmu*fMu*sqr(k_)/epsilon_*S2 - (nonlinearStress && gradU);

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1*G*epsilon_/k_
        // E-term
      + C2*f2*pow(Cmu, 0.75)*pow(k_, scalar(0.5))
       /(kappa_*y_*(scalar(1) - exp(-Aepsilon*yStar)))
       *exp(-Amu*sqr(yStar))*epsilon_
      - fvm::Sp(C2*f2*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "LienCubicKELowReSetWallDissipation.H"
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

    C5viscosity =
        - 2.0*pow(Cmu, 3.0)*pow(k_, 4.0)/pow(epsilon_, 3.0)*
        (magSqr(gradU + gradU.T()) - magSqr(gradU - gradU.T()));

    nut_ =
        Cmu*fMu*sqr(k_)/epsilon_
        // C5 term, implicit
        + max
        (
            C5viscosity,
            dimensionedScalar("0", C5viscosity.dimensions(), 0.0)
        );

    nonlinearStress = symm
    (
        // quadratic terms
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
        // cubic term C4
      - 20.0*pow(k_, 4.0)/pow(epsilon_, 3.0)*
        pow(Cmu, 3.0)*
        (
            ((gradU & gradU) & gradU.T())
          + ((gradU & gradU.T()) & gradU.T())
          - ((gradU.T() & gradU) & gradU)
          - ((gradU.T() & gradU.T()) & gradU)
        )
        // cubic term C5, explicit part
      + min
        (
            C5viscosity,
            dimensionedScalar("0", C5viscosity.dimensions(), 0.0)
        )*gradU
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RAS
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
