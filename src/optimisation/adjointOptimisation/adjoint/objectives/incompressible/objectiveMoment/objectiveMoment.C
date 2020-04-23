/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "objectiveMoment.H"
#include "createZeroField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveMoment, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveMoment,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveMoment::objectiveMoment
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    momentPatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        ).sortedToc()
    ),
    momentDirection_(dict.get<vector>("direction")),
    rotationCentre_(dict.get<vector>("rotationCenter")),
    Aref_(dict.get<scalar>("Aref")),
    lRef_(dict.get<scalar>("lRef")),
    rhoInf_(dict.get<scalar>("rhoInf")),
    UInf_(dict.get<scalar>("UInf")),
    invDenom_(2./(rhoInf_*UInf_*UInf_*Aref_*lRef_)),
    stressXPtr_
    (
        Foam::createZeroFieldPtr<vector>
        (
            mesh_, "stressX", dimLength/sqr(dimTime)
        )
    ),
    stressYPtr_
    (
        Foam::createZeroFieldPtr<vector>
        (
            mesh_, "stressY", dimLength/sqr(dimTime)
        )
    ),
    stressZPtr_
    (
        Foam::createZeroFieldPtr<vector>
        (
            mesh_, "stressZ", dimLength/sqr(dimTime)
        )
    ),
    devReff_(vars_.turbulence()->devReff()())
{
    // Sanity check and print info
    if (momentPatches_.empty())
    {
        FatalErrorInFunction
            << "No valid patch name on which to minimize " << type() << endl
            << exit(FatalError);
    }
    if (debug)
    {
        Info<< "Minimizing " << type() << " in patches:" << endl;
        for (const label patchI : momentPatches_)
        {
            Info<< "\t " << mesh_.boundary()[patchI].name() << endl;
        }
    }

    // Allocate boundary field pointers
    bdJdpPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdxdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdxdbDirectMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveMoment::J()
{
    vector pressureMoment(Zero);
    vector viscousMoment(Zero);
    vector cumulativeMoment(Zero);

    // Update field here and use the same value for all functions
    const volScalarField& p = vars_.pInst();
    devReff_ = vars_.turbulence()->devReff()();

    for (const label patchI : momentPatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const vectorField& Sf = patch.Sf();
        vectorField dx(patch.Cf() - rotationCentre_);
        pressureMoment += gSum
        (
            rhoInf_*(dx ^ Sf)*p.boundaryField()[patchI]
        );

        // Viscous term calculated using the full tensor derivative
        viscousMoment += gSum
        (
            rhoInf_*(dx^(devReff_.boundaryField()[patchI] & Sf))
        );
    }

    cumulativeMoment = pressureMoment + viscousMoment;

    scalar moment = cumulativeMoment & momentDirection_;
    scalar Cm = moment*invDenom_;
    DebugInfo<<
        "Moment|Coeff " << moment << "|" << Cm << endl;
    J_ = Cm;
    return Cm;
}


void objectiveMoment::update_meanValues()
{
    if (computeMeanFields_)
    {
        const volVectorField& U = vars_.U();
        const autoPtr<incompressible::RASModelVariables>&
           turbVars = vars_.RASModelVariables();
        const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();

        devReff_ = turbVars->devReff(lamTransp, U)();
    }
}


void objectiveMoment::update_boundarydJdp()
{
    for (const label patchI : momentPatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        vectorField dx(patch.Cf() - rotationCentre_);
        bdJdpPtr_()[patchI] = (momentDirection_ ^ dx)*invDenom_*rhoInf_;
    }
}


void objectiveMoment::update_dSdbMultiplier()
{
    const volScalarField& p = vars_.p();

    for (const label patchI : momentPatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const vectorField dx(patch.Cf() - rotationCentre_);
        bdSdbMultPtr_()[patchI] =
            (
                (
                    rhoInf_*
                    (
                        (momentDirection_^dx) &
                        (
                            devReff_.boundaryField()[patchI]
                        )
                    )
                )
              + rhoInf_ * (momentDirection_^dx) * p.boundaryField()[patchI]
            )
           *invDenom_;
    }
}


void objectiveMoment::update_dxdbMultiplier()
{
    const volScalarField& p = vars_.p();
    const volVectorField& U = vars_.U();

    const autoPtr<incompressible::RASModelVariables>&
       turbVars = vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();
    volScalarField nuEff(lamTransp.nu() + turbVars->nutRef());
    volTensorField gradU(fvc::grad(U));

    // Explicitly correct the boundary gradient to get rid of the
    // tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
            const vectorField& nf = tnf();
            gradU.boundaryFieldRef()[patchI] =
                nf * U.boundaryField()[patchI].snGrad();
        }
    }

    volTensorField stress(nuEff*(gradU + T(gradU)));

    stressXPtr_().replace(0, stress.component(0));
    stressXPtr_().replace(1, stress.component(1));
    stressXPtr_().replace(2, stress.component(2));

    stressYPtr_().replace(0, stress.component(3));
    stressYPtr_().replace(1, stress.component(4));
    stressYPtr_().replace(2, stress.component(5));

    stressZPtr_().replace(0, stress.component(6));
    stressZPtr_().replace(1, stress.component(7));
    stressZPtr_().replace(2, stress.component(8));

    volTensorField gradStressX(fvc::grad(stressXPtr_()));
    volTensorField gradStressY(fvc::grad(stressYPtr_()));
    volTensorField gradStressZ(fvc::grad(stressZPtr_()));

    volVectorField gradp(fvc::grad(p));

    for (const label patchI : momentPatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();
        vectorField dx(patch.Cf() - rotationCentre_);
        vectorField aux(momentDirection_^dx);
        bdxdbMultPtr_()[patchI] =
        (
            (
                (
                   -(aux.component(0) * gradStressX.boundaryField()[patchI])
                   -(aux.component(1) * gradStressY.boundaryField()[patchI])
                   -(aux.component(2) * gradStressZ.boundaryField()[patchI])
                ) & nf
            )
            + (momentDirection_ & (dx^nf))*gradp.boundaryField()[patchI]
        )
        *invDenom_*rhoInf_;
    }
}


void objectiveMoment::update_dxdbDirectMultiplier()
{
    const volScalarField& p = vars_.p();

    for (const label patchI : momentPatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();
        const vectorField dx(patch.Cf() - rotationCentre_);
        const vectorField force
        (
            rhoInf_
           *(
                ((p.boundaryField()[patchI]*nf)
              + (devReff_.boundaryField()[patchI] & nf))
            )
        );
        bdxdbDirectMultPtr_()[patchI] =
            (force^momentDirection_)*invDenom_*rhoInf_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
