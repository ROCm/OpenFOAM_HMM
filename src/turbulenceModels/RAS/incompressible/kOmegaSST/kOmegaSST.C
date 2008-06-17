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

#include "kOmegaSST.H"
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

defineTypeNameAndDebug(kOmegaSST, 0);
addToRunTimeSelectionTable(RASmodel, kOmegaSST, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> kOmegaSST::F1(const volScalarField& CDkOmega) const
{
    volScalarField CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    volScalarField arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

tmp<volScalarField> kOmegaSST::F2() const
{
    volScalarField arg2 = min
    (
        max
        (
            (scalar(2)/betaStar)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSST::kOmegaSST
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASmodel(typeName, U, phi, lamTransportModel),

    alphaK1
    (
        RASmodelCoeffs_.lookupOrAddDefault<scalar>("alphaK1", 0.85034)
    ),
    alphaK2(RASmodelCoeffs_.lookupOrAddDefault<scalar>("alphaK2", 1.0)),
    alphaOmega1
    (
        RASmodelCoeffs_.lookupOrAddDefault<scalar>("alphaOmega1", 0.5)
    ),
    alphaOmega2
    (
        RASmodelCoeffs_.lookupOrAddDefault<scalar>("alphaOmega2", 0.85616)
    ),
    gamma1(RASmodelCoeffs_.lookupOrAddDefault<scalar>("gamma1", 0.5532)),
    gamma2(RASmodelCoeffs_.lookupOrAddDefault<scalar>("gamma2", 0.4403)),
    beta1(RASmodelCoeffs_.lookupOrAddDefault<scalar>("beta1", 0.075)),
    beta2(RASmodelCoeffs_.lookupOrAddDefault<scalar>("beta2", 0.0828)),
    betaStar
    (
        RASmodelCoeffs_.lookupOrAddDefault<scalar>("betaStar", 0.09)
    ),
    a1(RASmodelCoeffs_.lookupOrAddDefault<scalar>("a1", 0.31)),
    c1(RASmodelCoeffs_.lookupOrAddDefault<scalar>("c1", 10.0)),

    omega0_("omega0", dimless/dimTime, SMALL),
    omegaSmall_("omegaSmall", dimless/dimTime, SMALL),

    Cmu(RASmodelCoeffs_.lookupOrAddDefault<scalar>("Cmu", 0.09)),

    y_(mesh_),

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

    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nut_(a1*k_/max(a1*(omega_ + omegaSmall_), F2()*mag(symm(fvc::grad(U_)))))
{
#   include "kOmegaWallViscosityI.H"

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaSST::R() const
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


tmp<volSymmTensorField> kOmegaSST::devReff() const
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


tmp<fvVectorMatrix> kOmegaSST::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool kOmegaSST::read()
{
    if (RASmodel::read())
    {
        RASmodelCoeffs_.readIfPresent<scalar>("alphaK1", alphaK1);
        RASmodelCoeffs_.readIfPresent<scalar>("alphaK2", alphaK2);
        RASmodelCoeffs_.readIfPresent<scalar>
        (
            "alphaOmega1",
            alphaOmega1
        );
        RASmodelCoeffs_.readIfPresent<scalar>
        (
            "alphaOmega2",
            alphaOmega2
        );
        RASmodelCoeffs_.readIfPresent<scalar>("gamma1", gamma1);
        RASmodelCoeffs_.readIfPresent<scalar>("gamma2", gamma2);
        RASmodelCoeffs_.readIfPresent<scalar>("beta1", beta1);
        RASmodelCoeffs_.readIfPresent<scalar>("beta2", beta2);
        RASmodelCoeffs_.readIfPresent<scalar>("betaStar", betaStar);
        RASmodelCoeffs_.readIfPresent<scalar>("a1", a1);
        RASmodelCoeffs_.readIfPresent<scalar>("c1", c1);
        RASmodelCoeffs_.readIfPresent<scalar>("Cmu", Cmu);

        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaSST::correct()
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

    volScalarField S2 = magSqr(symm(fvc::grad(U_)));
    volScalarField G = nut_*2*S2;

#   include "kOmegaWallFunctionsI.H"

    volScalarField CDkOmega =
        (2*alphaOmega2)*(fvc::grad(k_) & fvc::grad(omega_))/omega_;

    volScalarField F1 = this->F1(CDkOmega);

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)*2*S2
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();

#   include "wallOmegaI.H"

    solve(omegaEqn);
    bound(omega_, omega0_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G, c1*betaStar*k_*omega_)
      - fvm::Sp(betaStar*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    nut_ = a1*k_/max(a1*omega_, F2()*sqrt(S2));

#   include "kOmegaWallViscosityI.H"

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RAS
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
