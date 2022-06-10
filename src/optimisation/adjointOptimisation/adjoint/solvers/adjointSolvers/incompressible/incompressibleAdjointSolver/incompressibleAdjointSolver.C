/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
        ATCModel_().updatePrimalBasedQuantities();
        getAdjointVars().updatePrimalBasedQuantities();
    }
}


Foam::tmp<Foam::volTensorField>
Foam::incompressibleAdjointSolver::computeGradDxDbMultiplier()
{
    /*
    addProfiling
    (
        incompressibleAdjointSolver,
        "incompressibleAdjointSolver::computeGradDxDbMultiplier"
    );
    */
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS
    (
        getAdjointVars().adjointTurbulence()
    );

    const volScalarField& p = primalVars_.p();
    const volVectorField& U = primalVars_.U();
    const volScalarField& pa = getAdjointVars().pa();
    const volVectorField& Ua = getAdjointVars().Ua();

    // We only need to modify the boundaryField of gradU locally.
    // If grad(U) is cached then
    // a. The .ref() call fails since the tmp is initialised from a
    //    const ref
    // b. we would be changing grad(U) for all other places in the code
    //    that need it
    // So, always allocate new memory and avoid registering the new field
    tmp<volTensorField> tgradU =
        volTensorField::New("gradULocal", fvc::grad(U));
    volTensorField& gradU = tgradU.ref();
    volTensorField::Boundary& gradUbf = gradU.boundaryFieldRef();

    // Explicitly correct the boundary gradient to get rid of
    // the tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
            gradUbf[patchI] = tnf*U.boundaryField()[patchI].snGrad();
        }
    }

    tmp<volScalarField> tnuEff = adjointRAS->nuEff();
    tmp<volSymmTensorField> stress = tnuEff()*twoSymm(gradU);
    // Note:
    // term4 (Ua & grad(stress)) is numerically tricky.  Its div leads to third
    // order spatial derivs in E-SI based computations. Applying the product
    // derivative rule (putting Ua inside the grad) gives better results in
    // NACA0012, SA, WF.  However, the original formulation should be kept at
    // the boundary in order to respect the Ua boundary conditions (necessary
    // for E-SI to give the same sens as FI).  A mixed approach is hence
    // followed

    // Term 3, used also to allocated the return field
    tmp<volTensorField> tgradUa = fvc::grad(Ua);
    auto tflowTerm =
        tmp<volTensorField>::New
        (
            "flowTerm",
          - tnuEff*(gradU & twoSymm(tgradUa()))
        );
    volTensorField& flowTerm = tflowTerm.ref();
    // Term 4, only for the internal field
    flowTerm.ref() +=
        (
          - (tgradUa & stress())
          + fvc::grad(Ua & stress())
        )().internalField();

    // Boundary conditions from term 4
    for (label idir = 0; idir < pTraits<vector>::nComponents; ++idir)
    {
        autoPtr<volVectorField> stressDirPtr
        (
            createZeroFieldPtr<vector>
                (mesh_, "stressDir", stress().dimensions())
        );
        // Components need to be in the [0-5] range since stress is a
        // volSymmTensorField
        unzipRow(stress(), idir, stressDirPtr());
        volTensorField gradStressDir(fvc::grad(stressDirPtr()));
        forAll(mesh_.boundary(), pI)
        {
            if (!isA<coupledFvPatch>(mesh_.boundary()[pI]))
            {
                flowTerm.boundaryFieldRef()[pI] +=
                    Ua.component(idir)().boundaryField()[pI]
                   *gradStressDir.boundaryField()[pI];
            }
        }
    }
    // Release memory
    stress.clear();

    // Compute dxdb multiplier
    flowTerm +=
        // Term 1, ATC
        ATCModel_->getFISensitivityTerm()
        // Term 2
      - fvc::grad(p)*Ua;

    // Term 5
    flowTerm += pa*tgradU;

    // Term 6, from the adjoint turbulence model
    flowTerm += T(adjointRAS->FISensitivityTerm());

    // Term 7, term from objective functions
    PtrList<objective>& functions
        (objectiveManagerPtr_->getObjectiveFunctions());

    for (objective& objI : functions)
    {
        if (objI.hasGradDxDbMult())
        {
            flowTerm += objI.weight()*objI.gradDxDbMultiplier();
        }
    }

    flowTerm.correctBoundaryConditions();

  //profiling::writeNow();

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
