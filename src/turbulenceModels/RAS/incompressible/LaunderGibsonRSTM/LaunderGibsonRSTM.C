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

#include "LaunderGibsonRSTM.H"
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

defineTypeNameAndDebug(LaunderGibsonRSTM, 0);
addToRunTimeSelectionTable(turbulenceModel, LaunderGibsonRSTM, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
LaunderGibsonRSTM::LaunderGibsonRSTM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    turbulenceModel(typeName, U, phi, lamTransportModel),

    Cmu(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Cmu", 0.09)),
    Clg1(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Clg1", 1.8)),
    Clg2(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Clg2", 0.6)),
    C1(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C1", 1.44)),
    C2(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C2", 1.92)),
    Cs(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Cs", 0.25)),
    Ceps(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("Ceps", 0.15)),
    alphaR(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphaR", 1.22)),
    alphaEps
    (
        turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("alphaEps", 0.76923)
    ),
    C1Ref(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C1Ref", 0.5)),
    C2Ref(turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("C2Ref", 0.3)),
    couplingFactor_
    (
        turbulenceModelCoeffs_.lookupOrAddDefault<scalar>("couplingFactor", 0.0)
    ),

    yr_(mesh_),

    R_
    (
        IOobject
        (
            "R",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
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

    nut_(Cmu*sqr(k_)/(epsilon_ + epsilonSmall_))
{
#   include "wallViscosityI.H"

    if (couplingFactor_ < 0.0 || couplingFactor_ > 1.0)
    {
        FatalErrorIn
        (
            "LaunderGibsonRSTM::LaunderGibsonRSTM"
            "(const volVectorField& U, const surfaceScalarField& phi,"
            "transportModel& lamTransportModel)"
        )   << "couplingFactor = " << couplingFactor_
            << " is not in range 0 - 1" << nl
            << exit(FatalError);
    }

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LaunderGibsonRSTM::devReff() const
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
            R_ - nu()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> LaunderGibsonRSTM::divDevReff(volVectorField& U) const
{
    if (couplingFactor_ > 0.0)
    {
        return
        (
            fvc::div(R_ + couplingFactor_*nut_*fvc::grad(U), "div(R)")
          + fvc::laplacian((1.0-couplingFactor_)*nut_, U, "laplacian(nuEff,U)")
          - fvm::laplacian(nuEff(), U)
        );
    }
    else
    {
        return
        (
            fvc::div(R_)
          + fvc::laplacian(nut_, U, "laplacian(nuEff,U)")
          - fvm::laplacian(nuEff(), U)
        );
    }
}


bool LaunderGibsonRSTM::read()
{
    if (turbulenceModel::read())
    {
        turbulenceModelCoeffs_.readIfPresent<scalar>("Cmu", Cmu);
        turbulenceModelCoeffs_.readIfPresent<scalar>("Clg1", Clg1);
        turbulenceModelCoeffs_.readIfPresent<scalar>("Clg2", Clg2);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C1", C1);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C2", C2);
        turbulenceModelCoeffs_.readIfPresent<scalar>("Cs", Cs);
        turbulenceModelCoeffs_.readIfPresent<scalar>("Ceps", Ceps);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphaR", alphaR);
        turbulenceModelCoeffs_.readIfPresent<scalar>("alphaEps", alphaEps);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C1Ref", C1Ref);
        turbulenceModelCoeffs_.readIfPresent<scalar>("C2Ref", C2Ref);

        turbulenceModelCoeffs_.readIfPresent<scalar>
        (
            "couplingFactor",
            couplingFactor_
        );

        if (couplingFactor_ < 0.0 || couplingFactor_ > 1.0)
        {
            FatalErrorIn("LaunderGibsonRSTM::read()")
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1"
                << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


void LaunderGibsonRSTM::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    turbulenceModel::correct();

    if (mesh_.changing())
    {
        yr_.correct();
    }

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G = 0.5*tr(P);

#   include "wallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
    //- fvm::laplacian(Ceps*(k_/epsilon_)*R_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1*G*epsilon_/k_
      - fvm::Sp(C2*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "wallDissipationI.H"

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Reynolds stress equation

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (typeid(curPatch) == typeid(wallFvPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                P[faceCelli] *=
                    min(G[faceCelli]/(0.5*tr(P[faceCelli]) + SMALL), 1.0);
            }
        }
    }

    volSymmTensorField reflect = C1Ref*epsilon_/k_*R_ - C2Ref*Clg2*dev(P);

    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(R_)
      + fvm::div(phi_, R_)
    //- fvm::laplacian(Cs*(k_/epsilon_)*R_, R_)
      - fvm::laplacian(DREff(), R_)
      + fvm::Sp(Clg1*epsilon_/k_, R_)
      ==
        P
      + (2.0/3.0*(Clg1 - 1)*I)*epsilon_
      - Clg2*dev(P)

        // wall reflection terms
      + symm
        (
            I*((yr_.n() & reflect) & yr_.n())
          - 1.5*(yr_.n()*(reflect & yr_.n())
          + (yr_.n() & reflect)*yr_.n())
        )*pow(Cmu, 0.75)*pow(k_, 1.5)/(kappa_*yr_*epsilon_)
    );

    REqn().relax();
    solve(REqn);

    R_.max
    (
        dimensionedSymmTensor
        (
            "zero",
            R_.dimensions(),
            symmTensor
            (
                k0_.value(), -GREAT, -GREAT,
                             k0_.value(), -GREAT,
                                          k0_.value()
            )
        )
    );

    k_ == 0.5*tr(R_);
    bound(k_, k0_);


    // Re-calculate turbulent viscosity
    nut_ = Cmu*sqr(k_)/epsilon_;


#   include "wallViscosityI.H"


    // Correct wall shear stresses

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (typeid(curPatch) == typeid(wallFvPatch))
        {
            symmTensorField& Rw = R_.boundaryField()[patchi];

            const scalarField& nutw = nut_.boundaryField()[patchi];

            vectorField snGradU = U_.boundaryField()[patchi].snGrad();

            const vectorField& faceAreas
                = mesh_.Sf().boundaryField()[patchi];

            const scalarField& magFaceAreas
                = mesh_.magSf().boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                // Calculate near-wall velocity gradient
                tensor gradUw
                    = (faceAreas[facei]/magFaceAreas[facei])*snGradU[facei];

                // Calculate near-wall shear-stress tensor
                tensor tauw = -nutw[facei]*2*symm(gradUw);

                // Reset the shear components of the stress tensor
                Rw[facei].xy() = tauw.xy();
                Rw[facei].xz() = tauw.xz();
                Rw[facei].yz() = tauw.yz();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RAS
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
