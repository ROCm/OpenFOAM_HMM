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

#include "LienLeschzinerLowRe.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LienLeschzinerLowRe, 0);
addToRunTimeSelectionTable(turbulenceModel, LienLeschzinerLowRe, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
LienLeschzinerLowRe::LienLeschzinerLowRe
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    turbulenceModel(typeName, U, phi, lamTransportModel),

    C1(turbulenceModelCoeffs_.lookupOrDefault<scalar>("C1", 1.44)),
    C2(turbulenceModelCoeffs_.lookupOrDefault<scalar>("C2", 1.92)),
    alphak(turbulenceModelCoeffs_.lookupOrDefault<scalar>("alphak", 1.0)),
    alphaEps
    (
        turbulenceModelCoeffs_.lookupOrDefault<scalar>("alphaEps", 0.76923)
    ),
    Cmu(turbulenceModelCoeffs_.lookupOrDefault<scalar>("Cmu", 0.09)),
    Am(turbulenceModelCoeffs_.lookupOrDefault<scalar>("Am", 0.016)),
    Aepsilon(turbulenceModelCoeffs_.lookupOrDefault<scalar>("Aepsilon", 0.263)),
    Amu(turbulenceModelCoeffs_.lookupOrDefault<scalar>("Amu", 0.00222)),

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

    yStar_(sqrt(k_)*y_/nu() + SMALL),

    nut_
    (
        Cmu*(scalar(1) - exp(-Am*yStar_))
       /(scalar(1) - exp(-Aepsilon*yStar_) + SMALL)*sqr(k_)
       /(epsilon_ + epsilonSmall_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LienLeschzinerLowRe::R() const
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
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(gradU),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> LienLeschzinerLowRe::devReff() const
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
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> LienLeschzinerLowRe::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
    //- (fvc::grad(U) & fvc::grad(nuEff()))
      - fvc::div(nuEff()*fvc::grad(U)().T())
    );
}


bool LienLeschzinerLowRe::read()
{
    if (turbulenceModel::read())
    {
        C1 = turbulenceModelCoeffs_.lookupOrDefault<scalar>("C1", 1.44);
        C2 = turbulenceModelCoeffs_.lookupOrDefault<scalar>("C2", 1.92);
        alphak = turbulenceModelCoeffs_.lookupOrDefault<scalar>("alphak", 1.0);
        alphaEps = turbulenceModelCoeffs_.lookupOrDefault<scalar>
            (
                "alphaEps",
                0.76923
            );
        Cmu = turbulenceModelCoeffs_.lookupOrDefault<scalar>("Cmu", 0.09);
        Am = turbulenceModelCoeffs_.lookupOrDefault<scalar>("Am", 0.016);
        Aepsilon = turbulenceModelCoeffs_.lookupOrDefault<scalar>
            (
                "Aepsilon",
                0.263
            );
        Amu = turbulenceModelCoeffs_.lookupOrDefault<scalar>("Amu", 0.00222);

        return true;
    }
    else
    {
        return false;
    }
}


void LienLeschzinerLowRe::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    turbulenceModel::correct();

    if (mesh_.changing())
    {
        y_.correct();
    }

    scalar Cmu75 = pow(Cmu.value(), 0.75);

    volTensorField gradU = fvc::grad(U_);

    // generation term
    volScalarField S2 = symm(gradU) && gradU;

    yStar_ = sqrt(k_)*y_/nu() + SMALL;
    volScalarField Rt = sqr(k_)/(nu()*epsilon_);

    volScalarField fMu =
        (scalar(1) - exp(-Am*yStar_))
       /(scalar(1) - exp(-Aepsilon*yStar_) + SMALL);

    volScalarField f2 = scalar(1) - 0.3*exp(-sqr(Rt));

    volScalarField G = Cmu*fMu*sqr(k_)/epsilon_*S2;


    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1*G*epsilon_/k_
        // E-term
        + C2*f2*Cmu75*pow(k_, scalar(0.5))
        /(kappa_*y_*(scalar(1) - exp(-Aepsilon*yStar_)))
       *exp(-Amu*sqr(yStar_))*epsilon_
      - fvm::Sp(C2*f2*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "LienLeschzinerLowReSetWallDissipation.H"
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
    nut_ = Cmu*fMu*sqr(k_)/epsilon_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceModels
} // End namespace Foam

// ************************************************************************* //
