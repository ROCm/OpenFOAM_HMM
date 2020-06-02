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

#include "adjointSolverManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointSolverManager, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointSolverManager::adjointSolverManager
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict
)
:
    regIOobject
    (
        IOobject
        (
            "adjointSolverManager" + dict.dictName(),
            mesh.time().system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true //register object
        )
    ),
    mesh_(mesh),
    dict_(dict),
    managerName_(dict.dictName()),
    primalSolverName_(dict.get<word>("primalSolver")),
    adjointSolvers_(0),
    objectiveSolverIDs_(0),
    constraintSolverIDs_(0),
    operatingPointWeight_
    (
        dict.getOrDefault<scalar>("operatingPointWeight", 1)
    )
{
    const dictionary& adjointSolversDict = dict.subDict("adjointSolvers");

    const wordList adjSolverNames = adjointSolversDict.toc();
    adjointSolvers_.setSize(adjSolverNames.size());
    objectiveSolverIDs_.setSize(adjSolverNames.size());
    constraintSolverIDs_.setSize(adjSolverNames.size());
    label nObjectives(0);
    label nConstraints(0);
    forAll(adjSolverNames, namei)
    {
        adjointSolvers_.set
        (
            namei,
            adjointSolver::New
            (
                mesh_,
                managerType,
                adjointSolversDict.subDict(adjSolverNames[namei]),
                primalSolverName_
            )
        );

        if (adjointSolvers_[namei].isConstraint())
        {
            constraintSolverIDs_[nConstraints++] = namei;
        }
        else
        {
            objectiveSolverIDs_[nObjectives++] = namei;
        }
    }
    objectiveSolverIDs_.setSize(nObjectives);
    constraintSolverIDs_.setSize(nConstraints);

    Info<< "Found " << nConstraints
        << " adjoint solvers acting as constraints" << endl;

    // Having more than one non-aggregated objectives per operating point
    // is needlessly expensive. Issue a warning
    if (objectiveSolverIDs_.size() > 1)
    {
        WarningInFunction
            << "Number of adjoint solvers corresponding to objectives "
            << "is greater than 1 (" << objectiveSolverIDs_.size() << ")" << nl
            << "Consider aggregating your objectives to one" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adjointSolverManager::readDict(const dictionary& dict)
{
    dict_ = dict;

    const dictionary& adjointSolversDict = dict.subDict("adjointSolvers");

    // Note: only updating existing solvers
    for (adjointSolver& solver : adjointSolvers_)
    {
        solver.readDict(adjointSolversDict.subDict(solver.name()));
    }

    return true;
}


const Foam::word& Foam::adjointSolverManager::managerName() const
{
    return managerName_;
}


const Foam::word& Foam::adjointSolverManager::primalSolverName() const
{
    return primalSolverName_;
}


const Foam::dictionary& Foam::adjointSolverManager::dict() const
{
    return dict_;
}


const Foam::PtrList<Foam::adjointSolver>&
Foam::adjointSolverManager::adjointSolvers() const
{
    return adjointSolvers_;
}


Foam::PtrList<Foam::adjointSolver>&
Foam::adjointSolverManager::adjointSolvers()
{
    return adjointSolvers_;
}


Foam::scalar Foam::adjointSolverManager::operatingPointWeight() const
{
    return operatingPointWeight_;
}


Foam::label Foam::adjointSolverManager::nConstraints() const
{
    return constraintSolverIDs_.size();
}


Foam::label Foam::adjointSolverManager::nObjectives() const
{
    return objectiveSolverIDs_.size();
}


Foam::label Foam::adjointSolverManager::nAdjointSolvers() const
{
    return nConstraints() + nObjectives();
}


void Foam::adjointSolverManager::solveAdjointEquations()
{
    for (adjointSolver& solver : adjointSolvers_)
    {
        // Solve the adjoint equations taking into consideration the weighted
        // contribution of possibly multiple objectives
        solver.solve();
    }
}


Foam::tmp<Foam::scalarField>
Foam::adjointSolverManager::aggregateSensitivities()
{
    tmp<scalarField> tsens(new scalarField(0));
    scalarField& sens = tsens.ref();

    // Sum sensitivities from all objectives expect the constraints
    for (const label solveri : objectiveSolverIDs_)
    {
        // Sum contributions
        const scalarField& solverSens =
            adjointSolvers_[solveri].getObjectiveSensitivities();

        if (sens.empty())
        {
            sens = scalarField(solverSens.size(), Zero);
        }
        sens += solverSens;
    }

    return tsens;
}


Foam::PtrList<Foam::scalarField>
Foam::adjointSolverManager::constraintSensitivities()
{
    PtrList<scalarField> constraintSens(constraintSolverIDs_.size());
    forAll(constraintSens, cI)
    {
        label consI = constraintSolverIDs_[cI];
        constraintSens.set
        (
            cI,
            new scalarField(adjointSolvers_[consI].getObjectiveSensitivities())
        );
    }

    return constraintSens;
}


void Foam::adjointSolverManager::computeAllSensitivities()
{
    for (adjointSolver& adjSolver : adjointSolvers_)
    {
        adjSolver.computeObjectiveSensitivities();
    }
}


void Foam::adjointSolverManager::clearSensitivities()
{
    for (adjointSolver& adjSolver : adjointSolvers_)
    {
        adjSolver.clearSensitivities();
    }
}


Foam::scalar Foam::adjointSolverManager::objectiveValue()
{
    scalar objValue(Zero);
    for (const label solveri : objectiveSolverIDs_)
    {
        objectiveManager& objManager =
            adjointSolvers_[solveri].getObjectiveManager();
        objValue += objManager.print();
    }

    return objValue;
}


Foam::tmp<Foam::scalarField> Foam::adjointSolverManager::constraintValues()
{
    tmp<scalarField> tconstraintValues
    (
        new scalarField(constraintSolverIDs_.size(), Zero)
    );
    scalarField& constraintValues = tconstraintValues.ref();
    forAll(constraintValues, cI)
    {
        objectiveManager& objManager =
            adjointSolvers_[constraintSolverIDs_[cI]].getObjectiveManager();
        constraintValues[cI] = objManager.print();
    }

    return tconstraintValues;
}


void Foam::adjointSolverManager::updatePrimalBasedQuantities(const word& name)
{
    if (primalSolverName_ == name)
    {
        for (adjointSolver& solver : adjointSolvers_)
        {
            solver.updatePrimalBasedQuantities();
        }
    }
}


// ************************************************************************* //
