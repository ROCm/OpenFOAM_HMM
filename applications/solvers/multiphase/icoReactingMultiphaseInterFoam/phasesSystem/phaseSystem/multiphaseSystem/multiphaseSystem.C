/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "multiphaseSystem.H"

#include "fixedValueFvsPatchFields.H"
#include "Time.H"
#include "subCycle.H"
#include "fvcMeshPhi.H"

#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcDdt.H"
#include "fvcFlux.H"
#include "fvmDdt.H"
#include "fvcAverage.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "CMULES.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseSystem, 0);
    defineRunTimeSelectionTable(multiphaseSystem, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseSystem::multiphaseSystem
(
     const fvMesh& mesh
)
:
    phaseSystem(mesh),
    cAlphas_(mesh.solverDict("alpha").lookup("cAlphas")),
    ddtAlphaMax_(0.0),
    limitedPhiAlphas_(phaseModels_.size()),
    Su_(phaseModels_.size()),
    Sp_(phaseModels_.size())
{
    label phasei = 0;
    phases_.setSize(phaseModels_.size());
    forAllIters(phaseModels_, iter)
    {
        phaseModel& pm = iter()();
        phases_.set(phasei++, &pm);
    }

    // Initiate Su and Sp
    forAllConstIters(phaseModels_, iter)
    {
        const phaseModel& pm = iter()();

        Su_.insert
        (
            pm.name(),
            volScalarField::Internal
            (
                IOobject
                (
                    "Su" + pm.name(),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(dimless/dimTime, Zero)
            )
        );

        Sp_.insert
        (
            pm.name(),
            volScalarField::Internal
            (
                IOobject
                (
                    "Sp" + pm.name(),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(dimless/dimTime, Zero)
            )
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiphaseSystem::calculateSuSp()
{
    forAllConstIters(totalPhasePairs_, iter)
    {
        const phasePair& pair = iter()();

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& alpha1 = pair.phase1();
        const volScalarField& alpha2 = pair.phase2();

        tmp<volScalarField> tCoeffs1 = this->coeffs(phase1.name());
        const volScalarField&  coeffs1 = tCoeffs1();

        tmp<volScalarField> tCoeffs2 = this->coeffs(phase2.name());
        const volScalarField&  coeffs2 = tCoeffs2();

        // Phase 1 to phase 2
        const phasePairKey key12
        (
            phase1.name(),
            phase2.name(),
            true
        );


        tmp<volScalarField> tdmdt12(this->dmdt(key12));
        const volScalarField& dmdt12 = tdmdt12();

        // Phase 2 to phase 1
        const phasePairKey key21
        (
            phase2.name(),
            phase1.name(),
            true
        );

        tmp<volScalarField> tdmdt21(this->dmdt(key21));
        const volScalarField& dmdt21 = tdmdt21();

        volScalarField::Internal& SpPhase1 = Sp_[phase1.name()];

        volScalarField::Internal& SuPhase1 = Su_[phase1.name()];

        volScalarField::Internal& SpPhase2 = Sp_[phase2.name()];

        volScalarField::Internal& SuPhase2 = Su_[phase2.name()];

        const volScalarField dmdtNet(dmdt21 - dmdt12);

        const volScalarField coeffs12(coeffs1 - coeffs2);

        // NOTE: dmdtNet is distributed in terms =
        //  Source for phase 1 =
        //      dmdtNet/rho1
        //    - alpha1*dmdtNet(1/rho1 - 1/rho2)

        forAll(dmdtNet, celli)
        {
            scalar dmdt21 = dmdtNet[celli];
            scalar coeffs12Cell = coeffs12[celli];

            scalar alpha1Limited = max(min(alpha1[celli], 1.0), 0.0);

            // exp.
            SuPhase1[celli] += coeffs1[celli]*dmdt21;

            if (dmdt21 > 0)
            {
                if (coeffs12Cell > 0)
                {
                    // imp
                    SpPhase1[celli] -= dmdt21*coeffs12Cell;
                }
                else if (coeffs12Cell < 0)
                {
                    // exp
                    SuPhase1[celli] -=
                        dmdt21*coeffs12Cell*alpha1Limited;
                }
            }
            else if (dmdt21 < 0)
            {
                if (coeffs12Cell > 0)
                {
                    // exp
                    SuPhase1[celli] -=
                        dmdt21*coeffs12Cell*alpha1Limited;
                }
                else if (coeffs12Cell < 0)
                {
                    // imp
                    SpPhase1[celli] -= dmdt21*coeffs12Cell;
                }
            }
        }

        forAll(dmdtNet, celli)
        {
            scalar dmdt12 = -dmdtNet[celli];
            scalar coeffs21Cell = -coeffs12[celli];

            scalar alpha2Limited = max(min(alpha2[celli], 1.0), 0.0);

            // exp
            SuPhase2[celli] += coeffs2[celli]*dmdt12;

            if (dmdt12 > 0)
            {
                if (coeffs21Cell > 0)
                {
                    // imp
                    SpPhase2[celli] -= dmdt12*coeffs21Cell;
                }
                else if (coeffs21Cell < 0)
                {
                    // exp
                    SuPhase2[celli] -=
                        dmdt12*coeffs21Cell*alpha2Limited;
                }
            }
            else if (dmdt12 < 0)
            {
                if (coeffs21Cell > 0)
                {
                    // exp
                    SuPhase2[celli] -=
                        coeffs21Cell*dmdt12*alpha2Limited;
                }
                else if (coeffs21Cell < 0)
                {
                    // imp
                    SpPhase2[celli] -= dmdt12*coeffs21Cell;
                }
            }
        }

        // Update ddtAlphaMax
        ddtAlphaMax_ =
            max
            (
                ddtAlphaMax_.value(),
                max(gMax((dmdt21*coeffs1)()), gMax((dmdt12*coeffs2)()))
            );
    }
}


void Foam::multiphaseSystem::solve()
{
    const fvMesh& mesh = this->mesh();

    const dictionary& alphaControls = mesh.solverDict("alpha");
    label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));
    label nAlphaCorr(alphaControls.get<label>("nAlphaCorr"));
    mesh.solverDict("alpha").readEntry("cAlphas", cAlphas_);

    // Reset ddtAlphaMax
    ddtAlphaMax_ = dimensionedScalar(dimless, Zero);

    PtrList<surfaceScalarField> phiAlphaCorrs(phases_.size());

    const surfaceScalarField& phi = this->phi();

    surfaceScalarField phic(mag((phi)/mesh_.magSf()));

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    surfaceScalarField::Boundary& phicBf = phic.boundaryFieldRef();
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        label phasei = 0;
        for (phaseModel& phase1 : phases_)
        {
            const volScalarField& alpha1 = phase1;

            phiAlphaCorrs.set
            (
                phasei,
                new surfaceScalarField
                (
                    "phi" + alpha1.name() + "Corr",
                    fvc::flux
                    (
                        phi,
                        alpha1,
                        "div(phi," + alpha1.name() + ')'
                    )
                )
            );

            surfaceScalarField& phiAlphaCorr = phiAlphaCorrs[phasei];

            for (phaseModel& phase2 : phases_)
            {
                const volScalarField& alpha2 = phase2;

                if (&phase2 == &phase1) continue;

                const phasePairKey key12(phase1.name(), phase2.name());

                if (!cAlphas_.found(key12))
                {
                    FatalErrorInFunction
                        << "Phase compression factor (cAlpha) not found for : "
                        << key12
                        << exit(FatalError);
                }
                scalar cAlpha = cAlphas_.find(key12)();

                phic = min(cAlpha*phic, max(phic));

                surfaceScalarField phir(phic*nHatf(alpha1, alpha2));

                word phirScheme
                (
                    "div(phir," + alpha2.name() + ',' + alpha1.name() + ')'
                );

                phiAlphaCorr += fvc::flux
                (
                   -fvc::flux(-phir, alpha2, phirScheme),
                    alpha1,
                    phirScheme
                );
            }

            // Ensure that the flux at inflow BCs is preserved
            forAll(phiAlphaCorr.boundaryField(), patchi)
            {
                fvsPatchScalarField& phiAlphaCorrp =
                    phiAlphaCorr.boundaryFieldRef()[patchi];

                if (!phiAlphaCorrp.coupled())
                {
                    const scalarField& phi1p = phi.boundaryField()[patchi];
                    const scalarField& alpha1p =
                        alpha1.boundaryField()[patchi];

                    forAll(phiAlphaCorrp, facei)
                    {
                        if (phi1p[facei] < 0)
                        {
                            phiAlphaCorrp[facei] = alpha1p[facei]*phi1p[facei];
                        }
                    }
                }
            }

            ++phasei;
        }

        // Set Su and Sp to zero
        for (const phaseModel& phase : phases_)
        {
            Su_[phase.name()] = dimensionedScalar("Su", dimless/dimTime, Zero);
            Sp_[phase.name()] = dimensionedScalar("Sp", dimless/dimTime, Zero);

            // Add alpha*div(U)
            const volScalarField& alpha = phase;
            Su_[phase.name()] +=
                fvc::div(phi)*min(max(alpha, scalar(0)), scalar(1));
        }


        // Fill Su and Sp
        calculateSuSp();

        // Limit phiAlphaCorr on each phase
        phasei = 0;
        for (phaseModel& phase : phases_)
        {
            volScalarField& alpha1 = phase;

            surfaceScalarField& phiAlphaCorr = phiAlphaCorrs[phasei];

            volScalarField::Internal& Su = Su_[phase.name()];
            volScalarField::Internal& Sp = Sp_[phase.name()];

            MULES::limit
            (
                1.0/mesh_.time().deltaT().value(),
                geometricOneField(),
                alpha1,
                phi,
                phiAlphaCorr,
                Sp,
                Su,
                oneField(),
                zeroField(),
                true
            );
            ++phasei;
        }

        MULES::limitSum(phiAlphaCorrs);

        volScalarField sumAlpha
        (
            IOobject
            (
                "sumAlpha",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );

        phasei = 0;
        for (phaseModel& phase : phases_)
        {
            volScalarField& alpha1 = phase;

            const volScalarField::Internal& Su = Su_[phase.name()];

            const volScalarField::Internal& Sp = Sp_[phase.name()];

            surfaceScalarField& phiAlpha = phiAlphaCorrs[phasei];

            // Add a bounded upwind U-mean flux
            //phiAlpha += upwind<scalar>(mesh_, phi).flux(alpha1);
            fvScalarMatrix alpha1Eqn
            (
                fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              + fv::gaussConvectionScheme<scalar>
                (
                    mesh,
                    phi,
                    upwind<scalar>(mesh, phi)
                ).fvmDiv(phi, alpha1)
              ==
                 Su + fvm::Sp(Sp, alpha1)
            );

            alpha1Eqn.solve();

            phiAlpha += alpha1Eqn.flux();

            if (nAlphaSubCycles > 1)
            {
                for
                (
                    subCycle<volScalarField> alphaSubCycle
                    (
                        alpha1,
                        nAlphaSubCycles
                    );
                    !(++alphaSubCycle).end();
                )
                {
                    MULES::explicitSolve
                    (
                        geometricOneField(),
                        alpha1,
                        phi,
                        phiAlpha,
                        (alphaSubCycle.index()*Sp)(),
                        (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                        oneField(),
                        zeroField()
                    );

                    if (alphaSubCycle.index() == 1)
                    {
                        phase.alphaPhi() = phiAlpha;
                    }
                    else
                    {
                        phase.alphaPhi() += phiAlpha;
                    }
                }

                phase.alphaPhi() /= nAlphaSubCycles;
            }
            else
            {
                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi,
                    phiAlpha,
                    Sp,
                    Su,
                    oneField(),
                    zeroField()
                );

                phase.alphaPhi() = phiAlpha;
            }

            ++phasei;
        }

        if (acorr == nAlphaCorr - 1)
        {
            volScalarField sumAlpha
            (
                IOobject
                (
                    "sumAlpha",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            );

            // Reset rhoPhi
            rhoPhi_ = dimensionedScalar("rhoPhi", dimMass/dimTime, Zero);

            for (phaseModel& phase : phases_)
            {
                volScalarField& alpha1 = phase;
                sumAlpha += alpha1;

                // Update rhoPhi
                rhoPhi_ += fvc::interpolate(phase.rho()) * phase.alphaPhi();

            }

            Info<< "Phase-sum volume fraction, min, max = "
                << sumAlpha.weightedAverage(mesh_.V()).value()
                << ' ' << min(sumAlpha).value()
                << ' ' << max(sumAlpha).value()
                << endl;

            volScalarField sumCorr(1.0 - sumAlpha);

            for (phaseModel& phase : phases_)
            {
                volScalarField& alpha = phase;
                alpha += alpha*sumCorr;

                Info<< alpha.name() << " volume fraction = "
                    << alpha.weightedAverage(mesh.V()).value()
                    << "  Min(alpha) = " << min(alpha).value()
                    << "  Max(alpha) = " << max(alpha).value()
                    << endl;
            }
        }
    }
}


const Foam::UPtrList<Foam::phaseModel>& Foam::multiphaseSystem::phases() const
{
    return phases_;
}


Foam::UPtrList<Foam::phaseModel>& Foam::multiphaseSystem::phases()
{
    return phases_;
}


const Foam::phaseModel& Foam::multiphaseSystem::phase(const label i) const
{
    return phases_[i];
}


Foam::phaseModel& Foam::multiphaseSystem::phase(const label i)
{
    return phases_[i];
}


Foam::dimensionedScalar Foam::multiphaseSystem::ddtAlphaMax() const
{
    return ddtAlphaMax_;
}


Foam::scalar Foam::multiphaseSystem::maxDiffNo() const
{
    auto iter = phaseModels_.cbegin();

    scalar maxVal = max(iter()->diffNo()).value();

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        maxVal = max(maxVal, max(iter()->diffNo()).value());
    }

    return maxVal * mesh_.time().deltaT().value();
}


const Foam::multiphaseSystem::compressionFluxTable&
Foam::multiphaseSystem::limitedPhiAlphas() const
{
    return limitedPhiAlphas_;
}


Foam::multiphaseSystem::SuSpTable& Foam::multiphaseSystem::Su()
{
    return Su_;
}


Foam::multiphaseSystem::SuSpTable& Foam::multiphaseSystem::Sp()
{
    return Sp_;
}


bool Foam::multiphaseSystem::read()
{
    return true;
}


// ************************************************************************* //
