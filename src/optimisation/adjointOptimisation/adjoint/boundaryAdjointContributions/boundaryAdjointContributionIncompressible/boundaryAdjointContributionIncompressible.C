/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "boundaryAdjointContributionIncompressible.H"
#include "adjointRASModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundaryAdjointContributionIncompressible, 0);
addToRunTimeSelectionTable
(
    boundaryAdjointContribution,
    boundaryAdjointContributionIncompressible,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundaryAdjointContributionIncompressible::
boundaryAdjointContributionIncompressible
(
    const word& managerName,
    const word& adjointSolverName,
    const word& simulationType,
    const fvPatch& patch
)
:
    boundaryAdjointContribution
    (
        managerName,
        adjointSolverName,
        simulationType,
        patch
    ),
    objectiveManager_
    (
        patch_.patch().boundaryMesh().mesh().
            lookupObjectRef<objectiveManager>(managerName)
    ),
    primalVars_
    (
        patch_.patch().boundaryMesh().mesh().
            lookupObject<incompressibleAdjointSolver>(adjointSolverName).
                getPrimalVars()
    ),
    adjointSolver_
    (
        patch_.patch().boundaryMesh().mesh().
            lookupObject<incompressibleAdjointSolver>(adjointSolverName)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> boundaryAdjointContributionIncompressible::velocitySource()
{
    // Objective function contribution
    PtrList<objective>& objectives = objectiveManager_.getObjectiveFunctions();
    tmp<vectorField> tsource =
        sumContributions
        (
            objectives,
            &objectiveIncompressible::boundarydJdv
        );
    vectorField& source = tsource.ref();

    // Turbulence model differentiation contribution.
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
        adjointVars().adjointTurbulence();
    source += adjointRAS().adjointMomentumBCSource()[patch_.index()];

    return tsource;
}


tmp<scalarField> boundaryAdjointContributionIncompressible::pressureSource()
{
    // Objective function contribution
    PtrList<objective>& objectives = objectiveManager_.getObjectiveFunctions();
    tmp<scalarField> tsource =
        sumContributions
        (
            objectives,
            &objectiveIncompressible::boundarydJdvn
        );

    scalarField& source = tsource.ref();

    // Turbulence model differentiation contribution.
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
        adjointVars().adjointTurbulence();
    const vectorField& adjointTurbulenceContr =
        adjointRAS().adjointMomentumBCSource()[patch_.index()];

    tmp<vectorField> tnf = patch_.nf();
    const vectorField& nf = tnf();

    source += adjointTurbulenceContr & nf;

    return (tsource);
}


tmp<vectorField>
boundaryAdjointContributionIncompressible::tangentVelocitySource()
{
    // Objective function contribution
    PtrList<objective>& objectives = objectiveManager_.getObjectiveFunctions();
    tmp<vectorField> tsource =
        sumContributions
        (
            objectives,
            &objectiveIncompressible::boundarydJdvt
        );

    vectorField& source = tsource.ref();

    // Turbulence model differentiation contribution.
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
        adjointVars().adjointTurbulence();
    const vectorField& adjointTurbulenceContr =
        adjointRAS().adjointMomentumBCSource()[patch_.index()];

    tmp<vectorField> tnf = patch_.nf();
    const vectorField& nf = tnf();

    source += adjointTurbulenceContr - (adjointTurbulenceContr & nf)*nf;

    return (tsource);
}


tmp<vectorField>
boundaryAdjointContributionIncompressible::normalVelocitySource()
{
    PtrList<objective>& objectives = objectiveManager_.getObjectiveFunctions();
    tmp<vectorField> tsource =
        sumContributions
        (
            objectives,
            &objectiveIncompressible::boundarydJdp
        );

    return (tsource);
}


tmp<scalarField> boundaryAdjointContributionIncompressible::energySource()
{
    PtrList<objective>& objectives = objectiveManager_.getObjectiveFunctions();
    tmp<scalarField> tsource =
        sumContributions
        (
            objectives,
            &objectiveIncompressible::boundarydJdT
        );

    return (tsource);
}


tmp<scalarField>
boundaryAdjointContributionIncompressible::adjointTMVariable1Source()
{
    PtrList<objective>& objectives = objectiveManager_.getObjectiveFunctions();
    tmp<scalarField> tsource =
        sumContributions
        (
            objectives,
            &objectiveIncompressible::boundarydJdTMvar1
        );

    return (tsource);
}


tmp<scalarField>
boundaryAdjointContributionIncompressible::adjointTMVariable2Source()
{
    PtrList<objective>& objectives = objectiveManager_.getObjectiveFunctions();
    tmp<scalarField> tsource =
        sumContributions
        (
            objectives,
            &objectiveIncompressible::boundarydJdTMvar2
        );

    return (tsource);
}


tmp<scalarField> boundaryAdjointContributionIncompressible::momentumDiffusion()
{
    tmp<scalarField> tnuEff(new scalarField(patch_.size(), Zero));
    scalarField& nuEff = tnuEff.ref();

    const autoPtr<incompressibleAdjoint::adjointRASModel>&
        adjointTurbulenceModel = adjointVars().adjointTurbulence();

    nuEff = adjointTurbulenceModel().nuEff()().boundaryField()[patch_.index()];

    return tnuEff;
}


tmp<scalarField> boundaryAdjointContributionIncompressible::laminarDiffusivity()
{
    tmp<scalarField> tnu(new scalarField(patch_.size(), Zero));
    scalarField& nu = tnu.ref();

    const autoPtr<incompressible::turbulenceModel>& turbulenceModel =
        primalVars_.turbulence();

    nu = turbulenceModel().nu()().boundaryField()[patch_.index()];

    return tnu;
}


tmp<scalarField> boundaryAdjointContributionIncompressible::thermalDiffusion()
{
    /*
    const polyMesh& mesh = patch_.patch().boundaryMesh().mesh();
    const compressible::turbulenceModel& turbulenceModel =
        mesh.lookupObject<compressible::turbulenceModel>("turbulenceModel");
    tmp<scalarField> talphaEff = turbulenceModel.alphaEff(patch_.index());
    */

    tmp<scalarField> talphaEff(new scalarField(patch_.size(), Zero));

    WarningInFunction
        << "no abstract thermalDiffusion is implemented. Returning zero field";


    return talphaEff;
}


tmp<scalarField> boundaryAdjointContributionIncompressible::wallDistance()
{
    tmp<scalarField> twallDist(new scalarField(patch_.size(), Zero));
    scalarField& wallDist = twallDist.ref();

    wallDist = primalVars_.turbulence()->y()[patch_.index()];

    return twallDist;
}


tmp<scalarField>
boundaryAdjointContributionIncompressible::TMVariable1Diffusion()
{
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
        adjointVars().adjointTurbulence();

    tmp<scalarField> tdiffCoeff =
        adjointRAS().diffusionCoeffVar1(patch_.index());

    return tdiffCoeff;
}


tmp<scalarField>
boundaryAdjointContributionIncompressible::TMVariable2Diffusion()
{
    const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
        adjointVars().adjointTurbulence();

    tmp<scalarField> tdiffCoeff =
        adjointRAS().diffusionCoeffVar2(patch_.index());

    return tdiffCoeff;
}


tmp<scalarField> boundaryAdjointContributionIncompressible::TMVariable1()
{
    const autoPtr<incompressible::RASModelVariables>& RASVariables =
        primalVars_.RASModelVariables();
    tmp<scalarField> tboundField(new scalarField(patch_.size(), Zero));
    scalarField& boundField = tboundField.ref();

    boundField = RASVariables().TMVar1().boundaryField()[patch_.index()];

    return tboundField;
}


tmp<scalarField> boundaryAdjointContributionIncompressible::TMVariable2()
{
    const autoPtr<incompressible::RASModelVariables>& RASVariables =
        primalVars_.RASModelVariables();
    tmp<scalarField> tboundField(new scalarField(patch_.size(), Zero));
    scalarField& boundField = tboundField.ref();

    boundField = RASVariables().TMVar2().boundaryField()[patch_.index()];

    return tboundField;
}


const fvPatchVectorField& boundaryAdjointContributionIncompressible::Ub() const
{
    return primalVars_.U().boundaryField()[patch_.index()];
}


const fvPatchScalarField& boundaryAdjointContributionIncompressible::pb() const
{
    return primalVars_.p().boundaryField()[patch_.index()];
}


const fvsPatchScalarField&
boundaryAdjointContributionIncompressible::phib() const
{
    return primalVars_.phi().boundaryField()[patch_.index()];
}


const fvPatchScalarField&
boundaryAdjointContributionIncompressible::turbulentDiffusivity() const
{
    return
        primalVars_.RASModelVariables()().nutRef().boundaryField()
        [
            patch_.index()
        ];
}


const fvPatchVectorField& boundaryAdjointContributionIncompressible::Uab() const
{
    return adjointVars().UaInst().boundaryField()[patch_.index()];
}


const fvPatchScalarField& boundaryAdjointContributionIncompressible::pab() const
{
    return adjointVars().paInst().boundaryField()[patch_.index()];
}


const fvsPatchScalarField&
boundaryAdjointContributionIncompressible::phiab() const
{
    return adjointVars().phiaInst().boundaryField()[patch_.index()];
}


const word boundaryAdjointContributionIncompressible::primalSolverName() const
{
    return primalVars_.solverName();
}


const word boundaryAdjointContributionIncompressible::adjointSolverName() const
{
    return adjointVars().solverName();
}


const incompressibleVars&
boundaryAdjointContributionIncompressible::primalVars() const
{
    return primalVars_;
}


const incompressibleAdjointVars&
boundaryAdjointContributionIncompressible::adjointVars() const
{
    return adjointSolver_.getAdjointVars();
}


objectiveManager&
boundaryAdjointContributionIncompressible::getObjectiveManager()
{
    return objectiveManager_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
