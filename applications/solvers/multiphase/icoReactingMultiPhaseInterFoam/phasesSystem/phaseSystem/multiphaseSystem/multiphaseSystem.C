/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C)  2017 OpenCFD Ltd.
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
#include "MULES.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcDdt.H"
#include "fvcFlux.H"
#include "fvcAverage.H"
#include "fvMatrix.H"

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
    dmdtContinuityError_
    (
        IOobject
        (
            "dmdtContinuityError",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", inv(dimTime), 0)
    )
{
    label phaseI = 0;
    phases_.setSize(phaseModels_.size());
    forAllConstIter(HashTable<autoPtr<phaseModel> >, phaseModels_, iter)
    {
        phaseModel& pm = const_cast<phaseModel&>(iter()());
        phases_.set(phaseI++, &pm);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseSystem::~multiphaseSystem()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiphaseSystem::solve()
{
    const fvMesh& mesh = this->mesh();

    const dictionary& alphaControls = mesh.solverDict("alpha");
    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));
    scalar cAlpha(readScalar(alphaControls.lookup("cAlpha")));
    scalar icAlpha(readScalar(alphaControls.lookup("icAlpha")));

    PtrList<surfaceScalarField> phiAlphaCorrs(phases_.size());

    const surfaceScalarField& phi = this->phi();

    surfaceScalarField phic(mag((phi)/mesh_.magSf()));

    phic = min(cAlpha*phic, max(phic));

     // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (cAlpha*icAlpha)*fvc::interpolate(mag(this->U()));
    }

    forAllIter(UPtrList<phaseModel>, phases_, iter)
    {
        phaseModel& phase1 = iter();
        volScalarField& alpha1 = phase1;
        alpha1.correctBoundaryConditions();
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        int phasei = 0;
        forAllIter(UPtrList<phaseModel>, phases_, iter)
        {
            phaseModel& phase1 = iter();
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

            forAllIter(UPtrList<phaseModel>, phases_, iter2)
            {
                phaseModel& phase2 = iter2();
                const volScalarField& alpha2 = phase2;

                if (&phase2 == &phase1)
                {
                    continue;
                }

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

            phasei++;
        }

        // Set Su and Sp tp zero
        forAllIter(UPtrList<phaseModel>, phases_, iter)
        {
            phaseModel& phase = iter();
            Su_[phase.name()] = dimensionedScalar("Su", dimless/dimTime, 0.0);
            Sp_[phase.name()] = dimensionedScalar("Sp", dimless/dimTime, 0.0);

            // Add alpha*div(U)
            const volScalarField& alpha = phase;
            Su_[phase.name()] +=
                fvc::div(phi)*min(max(alpha, scalar(0)), scalar(1));
        }


        // Fill Su and Sp
        forAllIter(phaseSystem::phasePairTable, totalPhasePairs_, iter)
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

            volScalarField::Internal& SpPhase1 =
                Sp_[phase1.name()];

            volScalarField::Internal& SuPhase1 =
                Su_[phase1.name()];

            volScalarField::Internal& SpPhase2 =
                Sp_[phase2.name()];

            volScalarField::Internal& SuPhase2 =
                Su_[phase2.name()];

            const volScalarField dmdtNet(dmdt21 - dmdt12);

            const volScalarField coeffs12(coeffs1 - coeffs2);

            //DebugVar(max(coeffs1.internalField()));
            //DebugVar(min(coeffs1.internalField()));

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

            DebugVar(key12);
            DebugVar(key21);

            DebugVar(max(dmdt21.internalField()));
            DebugVar(max(dmdt12.internalField()));

            DebugVar(min(dmdt21.internalField()));
            DebugVar(min(dmdt12.internalField()));

            //DebugVar(max(SpPhase2));
            //DebugVar(max(SpPhase1));

            //DebugVar(min(SpPhase2));
            //DebugVar(min(SpPhase1));

            //DebugVar(max(SuPhase2));
            //DebugVar(max(SuPhase1));

            //DebugVar(min(SuPhase2));
            //DebugVar(min(SuPhase1));
        }


        // Limit phiAlphaCorr on each phase
        phasei = 0;
        forAllIter(UPtrList<phaseModel>, phases_, iter)
        {
            phaseModel& phase = iter();
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
                1,
                0,
                true
            );
            phasei ++;
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
            dimensionedScalar("sumAlpha", dimless, 0)
        );

        phasei = 0;
        forAllIter(UPtrList<phaseModel>, phases_, iter)
        {
            phaseModel& phase = iter();
            volScalarField& alpha1 = phase;


            const volScalarField::Internal& Su =
                Su_[phase.name()];

            const volScalarField::Internal& Sp =
                Sp_[phase.name()];

            surfaceScalarField& phiAlpha = phiAlphaCorrs[phasei];

            // Add a bounded upwind U-mean flux
            phiAlpha += upwind<scalar>(mesh_, phi).flux(alpha1);

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
                    //surfaceScalarField phiAlphaCorrs0(phiAlphaCorrs[phasei]);

                    MULES::explicitSolve
                    (
                        geometricOneField(),
                        alpha1,
                        phi,
                        phiAlpha,
                        (alphaSubCycle.index()*Sp)(),
                        (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                        phase.alphaMax(),
                        0
                    );

                    if (alphaSubCycle.index() == 1)
                    {
                        phase.alphaPhi() = phiAlpha;//phiAlphaCorrs0;
                    }
                    else
                    {
                        phase.alphaPhi() += phiAlpha;//phiAlphaCorrs0;
                    }
                }

                phase.alphaPhi() /= nAlphaSubCycles;

            }
            else
            {
                phaseModel& phase = iter();
                volScalarField& alpha1 = phase;

                //surfaceScalarField& phiAlpha = phiAlphaCorrs[phasei];

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi,
                    phiAlpha,
                    Sp,
                    Su,
                    phase.alphaMax(),
                    0
                );


                phase.alphaPhi() = phiAlpha;
            }

            phasei++;
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
                dimensionedScalar("sumAlpha", dimless, 0)
            );

            // Reset rhoPhi
            rhoPhi_ = dimensionedScalar("rhoPhi", dimMass/dimTime, 0.0);

            forAllIter(UPtrList<phaseModel>, phases_, iter)
            {
                phaseModel& phase = iter();
                volScalarField& alpha1 = phase;
                sumAlpha += alpha1;

                // Update rhoPhi
                rhoPhi_ +=
                    fvc::interpolate(phase.rho())*phase.alphaPhi();


                Info<< alpha1.name() << " volume fraction = "
                    << alpha1.weightedAverage(mesh.V()).value()
                    << "  Min(alpha) = " << min(alpha1).value()
                    << "  Max(alpha) = " << max(alpha1).value()
                << endl;

                //DebugVar(max(alpha1.internalField()));
                //DebugVar(min(alpha1.internalField()));

                volScalarField dAlphadt("dAlphadt", fvc::ddt(alpha1));

                //DebugVar(max(dAlphadt.internalField()));
                //DebugVar(min(dAlphadt.internalField()));
                if (mesh.time().outputTime())
                {
                    //volScalarField divrhoPhi("divrhoPhi", fvc::div(rhoPhi_));
                    //divrhoPhi.write();
                    //dAlphadt.write();
                }

            }

            Info<< "Phase-sum volume fraction, min, max = "
                << sumAlpha.weightedAverage(mesh_.V()).value()
                << ' ' << min(sumAlpha).value()
                << ' ' << max(sumAlpha).value()
                << endl;
        }
    }
}


const Foam::volScalarField& Foam::multiphaseSystem::dmdtContinuityError() const
{
    return dmdtContinuityError_;
}


Foam::volScalarField& Foam::multiphaseSystem::dmdtContinuityError()
{
    return dmdtContinuityError_;
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


// ************************************************************************* //
