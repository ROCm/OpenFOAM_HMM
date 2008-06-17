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
namespace incompressible
{
namespace RAS
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RNGkEpsilon, 0);
addToRunTimeSelectionTable(turbulenceModel, RNGkEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RNGkEpsilon::RNGkEpsilon
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    turbulenceModel(typeName, U, phi, lamTransportModel),

    Cmu(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Cmu", 0.0845)),
    C1(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C1", 1.42)),
    C2(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C2", 1.68)),
    alphak(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphak", 1.39)),
    alphaEps
    (
        turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphaEps", 1.39)
    ),
    eta0(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("eta0", 4.38)),
    beta(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("beta", 0.012)),

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

    nut_(Cmu*sqr(k_)/(epsilon_ + epsilonSmall_))
{
#   include "wallViscosityI.H"

    printCoeffs();
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
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> RNGkEpsilon::devReff() const
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


tmp<fvVectorMatrix> RNGkEpsilon::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool RNGkEpsilon::read()
{
    if (turbulenceModel::read())
    {
        turbulenceModelCoeffs_.readIfPresent<scalar>("Cmu", Cmu);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C1", C1);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C2", C2);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphak", alphak);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphaEps", alphaEps);
        turbulenceModelCoeffs_.readIfPresent<scalar>("eta0", eta0);
        turbulenceModelCoeffs_.readIfPresent<scalar>("beta", beta);

        return true;
    }
    else
    {
        return false;
    }
}


void RNGkEpsilon::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    turbulenceModel::correct();

    volScalarField S2 = 2*magSqr(symm(fvc::grad(U_)));

    volScalarField G = nut_*S2;

    volScalarField eta = sqrt(S2)*k_/epsilon_;
    volScalarField R =
        ((eta*(scalar(1) - eta/eta0))/(scalar(1) + beta*eta*sqr(eta)));

#   include "wallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        (C1 - R)*G*epsilon_/k_
    //- fvm::SuSp(R*G/k_, epsilon_)
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
        G - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    nut_ = Cmu*sqr(k_)/epsilon_;

#   include "wallViscosityI.H"

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RAS
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
