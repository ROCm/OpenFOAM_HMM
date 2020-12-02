/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 PCOpt/NTUA
    Copyright (C) 2020 FOSS GP
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

#include "shapeSensitivitiesIncompressible.H"
#include "adjointBoundaryConditions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(shapeSensitivities, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void shapeSensitivities::accumulateDirectSensitivityIntegrand(const scalar dt)
{
    // Accumulate direct sensitivities
    PtrList<objective>& functions(objectiveManager_.getObjectiveFunctions());
    for (const label patchI : sensitivityPatchIDs_)
    {
        const scalarField magSfDt(mesh_.boundary()[patchI].magSf()*dt);
        for (objective& func : functions)
        {
            const scalar wei(func.weight());
            dSfdbMult_()[patchI] += wei*func.dSdbMultiplier(patchI)*dt;
            dnfdbMult_()[patchI] += wei*func.dndbMultiplier(patchI)*magSfDt;
            dxdbDirectMult_()[patchI] +=
                wei*func.dxdbDirectMultiplier(patchI)*magSfDt;
        }
    }
}


void shapeSensitivities::accumulateBCSensitivityIntegrand(const scalar dt)
{
    auto& UaBoundary = adjointVars_.Ua().boundaryFieldRef();
    tmp<boundaryVectorField> DvDbMult(dvdbMult());

    // Accumulate sensitivities due to boundary conditions
    for (const label patchI : sensitivityPatchIDs_)
    {
        const scalarField magSfDt(mesh_.boundary()[patchI].magSf()*dt);
        fvPatchVectorField& Uab = UaBoundary[patchI];
        if (isA<adjointVectorBoundaryCondition>(Uab))
        {
            bcDxDbMult_()[patchI] +=
            (
                DvDbMult()[patchI]
              & refCast<adjointVectorBoundaryCondition>(Uab).dxdbMult()
            )*magSfDt;
        }
    }
}


tmp<boundaryVectorField> shapeSensitivities::dvdbMult() const
{
    tmp<boundaryVectorField>
        tres(createZeroBoundaryPtr<vector>(meshShape_).ptr());
    boundaryVectorField& res = tres.ref();

    // Grab references
    const volScalarField& pa = adjointVars_.pa();
    const volVectorField& Ua = adjointVars_.Ua();
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars_.adjointTurbulence();

    // Fields needed to calculate adjoint sensitivities
    const autoPtr<incompressible::RASModelVariables>&
       turbVars = primalVars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = primalVars_.laminarTransport();
    volScalarField nuEff(lamTransp.nu() + turbVars->nutRef());
    volTensorField gradUa(fvc::grad(Ua));

    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = meshShape_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();

        res[patchI] =
            (
                nuEff.boundaryField()[patchI]
              * (
                    Ua.boundaryField()[patchI].snGrad()
                  + (gradUa.boundaryField()[patchI] & nf)
                )
            )
          - (nf*pa.boundaryField()[patchI])
          + adjointTurbulence().adjointMomentumBCSource()[patchI];
    }

    return tres;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

shapeSensitivities::shapeSensitivities
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager
)
:
    adjointSensitivity
    (
        mesh,
        dict,
        primalVars,
        adjointVars,
        objectiveManager
    ),
    shapeSensitivitiesBase(mesh, dict),
    dSfdbMult_(createZeroBoundaryPtr<vector>(mesh_)),
    dnfdbMult_(createZeroBoundaryPtr<vector>(mesh_)),
    dxdbDirectMult_(createZeroBoundaryPtr<vector>(mesh_)),
    bcDxDbMult_(createZeroBoundaryPtr<vector>(mesh_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void shapeSensitivities::clearSensitivities()
{
    dSfdbMult_() = vector::zero;
    dnfdbMult_() = vector::zero;
    dxdbDirectMult_() = vector::zero;
    bcDxDbMult_() = vector::zero;

    adjointSensitivity::clearSensitivities();
    shapeSensitivitiesBase::clearSensitivities();
}


void shapeSensitivities::write(const word& baseName)
{
    adjointSensitivity::write(baseName);
    shapeSensitivitiesBase::write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
