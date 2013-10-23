/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "twoPhaseSystem.H"
#include "fvMatrix.H"
#include "PhaseIncompressibleTurbulenceModel.H"
#include "surfaceInterpolate.H"
#include "MULES.H"
#include "subCycle.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcCurl.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::twoPhaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phase1_
    (
        *this,
        *this,
        wordList(lookup("phases"))[0]
    ),

    phase2_
    (
        *this,
        *this,
        wordList(lookup("phases"))[1]
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->calcPhi()
    ),

    dgdt_
    (
        IOobject
        (
            "dgdt",
            mesh.time().timeName(),
            mesh
        ),
        pos(phase2_)*fvc::div(phi_)/max(phase2_, scalar(0.0001))
    ),

    sigma_
    (
        "sigma",
        dimensionSet(1, 0, -2, 0, 0),
        lookup("sigma")
    ),

    Cvm_
    (
        "Cvm",
        dimless,
        lookup("Cvm")
    ),

    Cl_
    (
        "Cl",
        dimless,
        lookup("Cl")
    ),

    drag1_
    (
        dragModel::New
        (
            subDict("drag"),
            phase1_,
            phase1_,
            phase2_
        )
    ),

    drag2_
    (
        dragModel::New
        (
            subDict("drag"),
            phase2_,
            phase2_,
            phase1_
        )
    ),

    heatTransfer1_
    (
        heatTransferModel::New
        (
            subDict("heatTransfer"),
            phase1_,
            phase1_,
            phase2_
        )
    ),

    heatTransfer2_
    (
        heatTransferModel::New
        (
            subDict("heatTransfer"),
            phase2_,
            phase2_,
            phase1_
        )
    ),

    dispersedPhase_(lookup("dispersedPhase")),

    residualPhaseFraction_
    (
        readScalar(lookup("residualPhaseFraction"))
    ),

    residualSlip_
    (
        "residualSlip",
        dimVelocity,
        lookup("residualSlip")
    )
{
    if
    (
        !(
            dispersedPhase_ == phase1_.name()
         || dispersedPhase_ == phase2_.name()
         || dispersedPhase_ == "both"
        )
    )
    {
        FatalErrorIn("twoPhaseSystem::twoPhaseSystem(const fvMesh& mesh)")
            << "invalid dispersedPhase " << dispersedPhase_
            << exit(FatalError);
    }

    Info << "dispersedPhase is " << dispersedPhase_ << endl;

    // Ensure the phase-fractions sum to 1
    phase2_.volScalarField::operator=(scalar(1) - phase1_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::rho() const
{
    return phase1_*phase1_.thermo().rho() + phase2_*phase2_.thermo().rho();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::U() const
{
    return phase1_*phase1_.U() + phase2_*phase2_.U();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::calcPhi() const
{
    return
        fvc::interpolate(phase1_)*phase1_.phi()
      + fvc::interpolate(phase2_)*phase2_.phi();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::dragCoeff() const
{
    tmp<volScalarField> tdragCoeff
    (
        new volScalarField
        (
            IOobject
            (
                "dragCoeff",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("dragCoeff", dimensionSet(1, -3, -1, 0, 0), 0)
        )
    );
    volScalarField& dragCoeff = tdragCoeff();

    volVectorField Ur(phase1_.U() - phase2_.U());
    volScalarField magUr(mag(Ur) + residualSlip_);

    if (dispersedPhase_ == phase1_.name())
    {
        dragCoeff = drag1().K(magUr);
    }
    else if (dispersedPhase_ == phase2_.name())
    {
        dragCoeff = drag2().K(magUr);
    }
    else if (dispersedPhase_ == "both")
    {
        dragCoeff =
        (
            phase2_*drag1().K(magUr)
          + phase1_*drag2().K(magUr)
        );
    }
    else
    {
        FatalErrorIn("twoPhaseSystem::dragCoeff()")
            << "dispersedPhase: " << dispersedPhase_ << " is incorrect"
            << exit(FatalError);
    }

    volScalarField alphaCoeff(max(phase1_*phase2_, residualPhaseFraction_));
    dragCoeff *= alphaCoeff;

    // Remove drag at fixed-flux boundaries
    forAll(phase1_.phi().boundaryField(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                phase1_.phi().boundaryField()[patchi]
            )
        )
        {
            dragCoeff.boundaryField()[patchi] = 0.0;
        }
    }

    return tdragCoeff;
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::liftForce
(
    const volVectorField& U
) const
{
    tmp<volVectorField> tliftForce
    (
        new volVectorField
        (
            IOobject
            (
                "liftForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector
            (
                "liftForce",
                dimensionSet(1, -2, -2, 0, 0),
                vector::zero
            )
        )
    );
    volVectorField& liftForce = tliftForce();

    volVectorField Ur(phase1_.U() - phase2_.U());

    liftForce =
        Cl_*(phase1_*phase1_.rho() + phase2_*phase2_.rho())
       *(Ur ^ fvc::curl(U));

    // Remove lift at fixed-flux boundaries
    forAll(phase1_.phi().boundaryField(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                phase1_.phi().boundaryField()[patchi]
            )
        )
        {
            liftForce.boundaryField()[patchi] = vector::zero;
        }
    }

    return tliftForce;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::heatTransferCoeff() const
{
    tmp<volScalarField> theatTransferCoeff
    (
        new volScalarField
        (
            IOobject
            (
                "heatTransferCoeff",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "heatTransferCoeff",
                dimensionSet(1, -1, -3, -1, 0),
                0
            )
        )
    );
    volScalarField& heatTransferCoeff = theatTransferCoeff();

    volVectorField Ur(phase1_.U() - phase2_.U());
    volScalarField magUr(mag(Ur) + residualSlip_);

    if (dispersedPhase_ == phase1_.name())
    {
        heatTransferCoeff = heatTransfer1().K(magUr);
    }
    else if (dispersedPhase_ == phase2_.name())
    {
        heatTransferCoeff = heatTransfer2().K(magUr);
    }
    else if (dispersedPhase_ == "both")
    {
        heatTransferCoeff =
        (
            phase2_*heatTransfer1().K(magUr)
          + phase1_*heatTransfer2().K(magUr)
        );
    }
    else
    {
        FatalErrorIn("twoPhaseSystem::heatTransferCoeff()")
            << "dispersedPhase: " << dispersedPhase_ << " is incorrect"
            << exit(FatalError);
    }

    volScalarField alphaCoeff(max(phase1_*phase2_, residualPhaseFraction_));
    heatTransferCoeff *= alphaCoeff;

    // Remove heatTransfer at fixed-flux boundaries
    forAll(phase1_.phi().boundaryField(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                phase1_.phi().boundaryField()[patchi]
            )
        )
        {
            heatTransferCoeff.boundaryField()[patchi] = 0.0;
        }
    }

    return theatTransferCoeff;
}


void Foam::twoPhaseSystem::solve()
{
    const Time& runTime = mesh_.time();

    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    const surfaceScalarField& phi1 = phase1_.phi();
    const surfaceScalarField& phi2 = phase2_.phi();

    const dictionary& alphaControls = mesh_.solverDict
    (
        alpha1.name()
    );

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));
    Switch implicitPhasePressure
    (
        alphaControls.lookupOrDefault<Switch>("implicitPhasePressure", false)
    );

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    alpha1.correctBoundaryConditions();


    surfaceScalarField phic("phic", phi_);
    surfaceScalarField phir("phir", phi1 - phi2);

    surfaceScalarField alpha1f(fvc::interpolate(max(alpha1, scalar(0))));

    tmp<surfaceScalarField> pPrimeByA;

    if (implicitPhasePressure)
    {
        const volScalarField& rAU1 = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName("rAU", phase1_.name())
        );
        const volScalarField& rAU2 = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName("rAU", phase2_.name())
        );

        pPrimeByA =
            fvc::interpolate((1.0/phase1_.rho())
           *rAU1*phase1_.turbulence().pPrime())
          + fvc::interpolate((1.0/phase2_.rho())
           *rAU2*phase2_.turbulence().pPrime());

        surfaceScalarField phiP
        (
            pPrimeByA()*fvc::snGrad(alpha1, "bounded")*mesh_.magSf()
        );

        phic += alpha1f*phiP;
        phir += phiP;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", dgdt_.dimensions(), 0.0)
        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi_)*min(alpha1, scalar(1))
        );

        forAll(dgdt_, celli)
        {
            if (dgdt_[celli] > 0.0 && alpha1[celli] > 0.0)
            {
                Sp[celli] -= dgdt_[celli]*alpha1[celli];
                Su[celli] += dgdt_[celli]*alpha1[celli];
            }
            else if (dgdt_[celli] < 0.0 && alpha1[celli] < 1.0)
            {
                Sp[celli] += dgdt_[celli]*(1.0 - alpha1[celli]);
            }
        }

        dimensionedScalar totalDeltaT = runTime.deltaT();
        if (nAlphaSubCycles > 1)
        {
            phase1_.phiAlpha() =
                dimensionedScalar("0", phase1_.phiAlpha().dimensions(), 0);
        }

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            surfaceScalarField alphaPhic1
            (
                fvc::flux
                (
                    phic,
                    alpha1,
                    alphaScheme
                )
              + fvc::flux
                (
                    -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                    alpha1,
                    alpharScheme
                )
            );

            // Ensure that the flux at inflow BCs is preserved
            forAll(alphaPhic1.boundaryField(), patchi)
            {
                fvsPatchScalarField& alphaPhic1p =
                    alphaPhic1.boundaryField()[patchi];

                if (!alphaPhic1p.coupled())
                {
                    const scalarField& phi1p = phi1.boundaryField()[patchi];
                    const scalarField& alpha1p = alpha1.boundaryField()[patchi];

                    forAll(alphaPhic1p, facei)
                    {
                        if (phi1p[facei] < 0)
                        {
                            alphaPhic1p[facei] = alpha1p[facei]*phi1p[facei];
                        }
                    }
                }
            }

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhic1,
                Sp,
                Su,
                1,
                0
            );

            if (nAlphaSubCycles > 1)
            {
                phase1_.phiAlpha() += (runTime.deltaT()/totalDeltaT)*alphaPhic1;
            }
            else
            {
                phase1_.phiAlpha() = alphaPhic1;
            }
        }

        if (implicitPhasePressure)
        {
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian(alpha1f*pPrimeByA, alpha1, "bounded")
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            phase1_.phiAlpha() += alpha1Eqn.flux();
        }

        phase2_.phiAlpha() = phi_ - phase1_.phiAlpha();
        alpha2 = scalar(1) - alpha1;

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
    }
}


void Foam::twoPhaseSystem::correct()
{
    phase1_.correct();
    phase2_.correct();
}


void Foam::twoPhaseSystem::correctTurbulence()
{
    phase1_.turbulence().correct();
    phase2_.turbulence().correct();
}


bool Foam::twoPhaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        readOK &= phase1_.read(*this);
        readOK &= phase2_.read(*this);

        lookup("sigma") >> sigma_;
        lookup("Cvm") >> Cvm_;
        lookup("Cl") >> Cl_;

        // drag1_->read(*this);
        // drag2_->read(*this);

        // heatTransfer1_->read(*this);
        // heatTransfer2_->read(*this);

        lookup("dispersedPhase") >> dispersedPhase_;
        lookup("residualPhaseFraction") >> residualPhaseFraction_;
        lookup("residualSlip") >> residualSlip_;

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
