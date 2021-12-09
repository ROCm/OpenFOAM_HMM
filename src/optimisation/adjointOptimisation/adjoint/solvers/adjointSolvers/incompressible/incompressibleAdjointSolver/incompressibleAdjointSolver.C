/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "incompressibleAdjointSolver.H"
#include "incompressiblePrimalSolver.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleAdjointSolver, 0);
    defineRunTimeSelectionTable(incompressibleAdjointSolver, dictionary);
    addToRunTimeSelectionTable
    (
        adjointSolver,
        incompressibleAdjointSolver,
        adjointSolver
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleAdjointSolver::incompressibleAdjointSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName
)
:
    adjointSolver(mesh, managerType, dict, primalSolverName),
    primalVars_
    (
        mesh.lookupObjectRef<incompressiblePrimalSolver>(primalSolverName).
            getIncoVars()
    ),
    ATCModel_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::incompressibleAdjointSolver>
Foam::incompressibleAdjointSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName
)
{
    const word solverType(dict.get<word>("solver"));
    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "incompressibleAdjointSolver",
            solverType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return
        autoPtr<incompressibleAdjointSolver>
        (
            ctorPtr(mesh, managerType, dict, primalSolverName)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::incompressibleAdjointSolver::readDict(const dictionary& dict)
{
    if (adjointSolver::readDict(dict))
    {
        return true;
    }

    return false;
}

bool Foam::incompressibleAdjointSolver::useSolverNameForFields() const
{
    return getAdjointVars().useSolverNameForFields();
}


const Foam::incompressibleVars&
Foam::incompressibleAdjointSolver::getPrimalVars() const
{
    return primalVars_;
}



const Foam::incompressibleAdjointVars&
Foam::incompressibleAdjointSolver::getAdjointVars() const
{
    const incompressibleAdjointVars& adjointVars =
        refCast<incompressibleAdjointVars>(const_cast<variablesSet&>(vars_()));
    return adjointVars;
}


Foam::incompressibleAdjointVars&
Foam::incompressibleAdjointSolver::getAdjointVars()
{
    incompressibleAdjointVars& adjointVars =
        refCast<incompressibleAdjointVars>(const_cast<variablesSet&>(vars_()));
    return adjointVars;
}



const Foam::autoPtr<Foam::ATCModel>&
Foam::incompressibleAdjointSolver::getATCModel() const
{
    return ATCModel_;
}


Foam::autoPtr<Foam::ATCModel>& Foam::incompressibleAdjointSolver::getATCModel()
{
    return ATCModel_;
}


void Foam::incompressibleAdjointSolver::updatePrimalBasedQuantities()
{
    if (vars_)
    {
        getAdjointVars().adjointTurbulence()->setChangedPrimalSolution();
    }
}


Foam::tmp<Foam::volTensorField>
Foam::incompressibleAdjointSolver::computeGradDxDbMultiplier()
{
    // Term depending on the adjoint turbulence model
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS
    (
        getAdjointVars().adjointTurbulence()
    );
    tmp<volTensorField> tturbulenceTerm(adjointRAS->FISensitivityTerm());
    volTensorField& turbulenceTerm = tturbulenceTerm.ref();

    // nu effective
    tmp<volScalarField> tnuEff(adjointRAS->nuEff());
    const volScalarField& nuEff = tnuEff();

    tmp<volTensorField> tflowTerm
    (
        new volTensorField
        (
            IOobject
            (
               "flowTerm",
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(sqr(dimLength)/pow3(dimTime), Zero)
        )
    );
    volTensorField& flowTerm = tflowTerm.ref();

    const volScalarField& p = primalVars_.p();
    const volVectorField& U = primalVars_.U();
    const volScalarField& pa = getAdjointVars().pa();
    const volVectorField& Ua = getAdjointVars().Ua();
    volTensorField gradU(fvc::grad(U));
    volTensorField gradUa(fvc::grad(Ua));

    // Explicitly correct the boundary gradient to get rid of
    // the tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
            const vectorField& nf = tnf();
            gradU.boundaryFieldRef()[patchI] =
                nf*U.boundaryField()[patchI].snGrad();
            //gradUa.boundaryField()[patchI] =
            //    nf*Ua.boundaryField()[patchI].snGrad();
        }
    }

    volTensorField stress(nuEff*(gradU + T(gradU)));
    autoPtr<volVectorField> stressXPtr
    (
        createZeroFieldPtr<vector>(mesh_, "stressX", stress.dimensions())
    );
    autoPtr<volVectorField> stressYPtr
    (
        createZeroFieldPtr<vector>(mesh_, "stressY", stress.dimensions())
    );
    autoPtr<volVectorField> stressZPtr
    (
        createZeroFieldPtr<vector>(mesh_, "stressZ", stress.dimensions())
    );

    stressXPtr().replace(0, stress.component(0));
    stressXPtr().replace(1, stress.component(1));
    stressXPtr().replace(2, stress.component(2));

    stressYPtr().replace(0, stress.component(3));
    stressYPtr().replace(1, stress.component(4));
    stressYPtr().replace(2, stress.component(5));

    stressZPtr().replace(0, stress.component(6));
    stressZPtr().replace(1, stress.component(7));
    stressZPtr().replace(2, stress.component(8));

    volTensorField gradStressX(fvc::grad(stressXPtr()));
    volTensorField gradStressY(fvc::grad(stressYPtr()));
    volTensorField gradStressZ(fvc::grad(stressZPtr()));

    // Contribution from objective functions and constraints
    volTensorField objectiveContributions
    (
        IOobject
        (
            "objectiveContributions",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimLength)/pow3(dimTime), Zero)
    );
    PtrList<objective>& functions
        (objectiveManagerPtr_->getObjectiveFunctions());
    forAll(functions, funcI)
    {
        objectiveContributions +=
            functions[funcI].weight()
           *functions[funcI].gradDxDbMultiplier();
    }

    // Note:
    // term4 (Ua & grad(stress)) is numerically tricky.  Its div leads to third
    // order spatial derivs in E-SI based computations Applying the product
    // derivative rule (putting Ua inside the grad) gives better results in
    // NACA0012, SA, WF.  However, the original formulation should be kept at
    // the boundary in order to respect the Ua boundary conditions (necessary
    // for E-SI to give the same sens as FI).  A mixed approach is hence
    // followed
    volTensorField term4
    (
      - nuEff*(gradUa & (gradU + T(gradU)))
      + fvc::grad(nuEff * Ua & (gradU + T(gradU)))
    );

    forAll(mesh_.boundary(), pI)
    {
        if (!isA<coupledFvPatch>(mesh_.boundary()[pI]))
        {
            term4.boundaryFieldRef()[pI] =
                Ua.component(0)().boundaryField()[pI]
               *gradStressX.boundaryField()[pI]
              + Ua.component(1)().boundaryField()[pI]
               *gradStressY.boundaryField()[pI]
              + Ua.component(2)().boundaryField()[pI]
               *gradStressZ.boundaryField()[pI];
        }
    }

    // Compute dxdb multiplier
    flowTerm =
        // Term 1, ATC
        ATCModel_->getFISensitivityTerm()
        // Term 2
      - fvc::grad(p) * Ua
        // Term 3
      - nuEff*(gradU & (gradUa + T(gradUa)))
        // Term 4
      + term4
        // Term 5
      + (pa * gradU)
        // Term 6, from the adjoint turbulence model
      + turbulenceTerm.T()
        // Term 7, term from objective functions
      + objectiveContributions;

    flowTerm.correctBoundaryConditions();

    return (tflowTerm);
}


void Foam::incompressibleAdjointSolver::additionalSensitivityMapTerms
(
    boundaryVectorField& sensitivityMap,
    const labelHashSet& patchIDs,
    const scalar dt
)
{
    // Does nothing in base
}


// ************************************************************************* //
