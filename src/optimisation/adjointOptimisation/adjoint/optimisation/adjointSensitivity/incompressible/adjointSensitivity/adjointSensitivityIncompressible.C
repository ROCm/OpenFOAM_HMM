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

#include "runTimeSelectionTables.H"
#include "adjointSensitivityIncompressible.H"
#include "boundaryAdjointContribution.H"
#include "incompressibleAdjointSolver.H"
#include "wallFvPatch.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointSensitivity, 0);
defineRunTimeSelectionTable(adjointSensitivity, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointSensitivity::adjointSensitivity
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager
)
:
    sensitivity(mesh, dict),
    derivatives_(0),
    primalVars_(primalVars),
    adjointVars_(adjointVars),
    objectiveManager_(objectiveManager)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<adjointSensitivity> adjointSensitivity::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "adjointSensitivity type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "adjointSensitivity",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<adjointSensitivity>
    (
        ctorPtr
        (
            mesh,
            dict,
            primalVars,
            adjointVars,
            objectiveManager
        )
    );
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

const scalarField& adjointSensitivity::calculateSensitivities()
{
    assembleSensitivities();
    write(type());
    return derivatives_;
}


const scalarField& adjointSensitivity::getSensitivities() const
{
    return derivatives_;
}


void adjointSensitivity::clearSensitivities()
{
    derivatives_ = scalar(0);
    if (fieldSensPtr_)
    {
        fieldSensPtr_().primitiveFieldRef() = scalar(0);
    }
}


void adjointSensitivity::write(const word& baseName)
{
    sensitivity::write(baseName);
}


tmp<volTensorField> adjointSensitivity::computeGradDxDbMultiplier()
{
    // Term depending on the adjoint turbulence model
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS
    (
        adjointVars_.adjointTurbulence()
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
    const volScalarField& pa = adjointVars_.pa();
    const volVectorField& Ua = adjointVars_.Ua();
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
    PtrList<objective>& functions(objectiveManager_.getObjectiveFunctions());
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

    const autoPtr<ATCModel>& ATCModel =
        mesh_.lookupObject<incompressibleAdjointSolver>
        (
            objectiveManager_.adjointSolverName()
        ).getATCModel();

    // Compute dxdb multiplier
    flowTerm =
        // Term 1, ATC
        ATCModel->getFISensitivityTerm()
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


tmp<volVectorField> adjointSensitivity::adjointMeshMovementSource()
{
    tmp<volTensorField> tgradDxDbMult = computeGradDxDbMultiplier();
    volTensorField& gradDxDbMult = tgradDxDbMult.ref();

    tmp<volVectorField> tadjointMeshMovementSource
    (
        new volVectorField
        (
            IOobject
            (
               "adjointMeshMovementSource",
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(gradDxDbMult.dimensions()/dimLength, Zero)
        )
    );

    volVectorField& source = tadjointMeshMovementSource.ref();

    source -= fvc::div(gradDxDbMult.T());

    // Terms from fvOptions
    fv::options::New(this->mesh_).postProcessSens
    (
        source.primitiveFieldRef(), adjointVars_.solverName()
    );

    return (tadjointMeshMovementSource);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
