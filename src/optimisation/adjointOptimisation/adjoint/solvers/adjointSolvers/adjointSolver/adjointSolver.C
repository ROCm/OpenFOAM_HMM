/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "adjointSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointSolver, 0);
    defineRunTimeSelectionTable(adjointSolver, adjointSolver);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointSolver::adjointSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName
)
:
    solver(mesh, managerType, dict),
    primalSolverName_(primalSolverName),
    objectiveManagerPtr_
    (
        objectiveManager::New
        (
            mesh,
            dict.subDict("objectives"),
            solverName_,
            primalSolverName
        )
    ),
    sensitivities_(nullptr),
    computeSensitivities_
    (
        dict.getOrDefault<bool>("computeSensitivities", true)
    ),
    isConstraint_(dict.getOrDefault<bool>("isConstraint", false))
{
    // Update objective-related quantities to get correct derivatives
    // in case of continuation
    objectiveManagerPtr_().update();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::adjointSolver> Foam::adjointSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName
)
{
    const word solverType(dict.get<word>("type"));

    auto* ctorPtr = adjointSolverConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "adjointSolver",
            solverType,
            *adjointSolverConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<adjointSolver>
    (
        ctorPtr(mesh, managerType, dict, primalSolverName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::primalSolver& Foam::adjointSolver::getPrimalSolver() const
{
    return mesh_.lookupObject<primalSolver>(primalSolverName_);
}


Foam::primalSolver& Foam::adjointSolver::getPrimalSolver()
{
    return
        const_cast<primalSolver&>
        (
            mesh_.lookupObject<primalSolver>(primalSolverName_)
        );
}


bool Foam::adjointSolver::readDict(const dictionary& dict)
{
    if (solver::readDict(dict))
    {
        computeSensitivities_ =
            dict.getOrDefault<bool>("computeSensitivities", true);

        objectiveManagerPtr_->readDict(dict.subDict("objectives"));

        return true;
    }

    return false;
}


const Foam::objectiveManager& Foam::adjointSolver::getObjectiveManager() const
{
    return objectiveManagerPtr_();
}


Foam::objectiveManager& Foam::adjointSolver::getObjectiveManager()
{
    return objectiveManagerPtr_();
}


bool Foam::adjointSolver::isConstraint()
{
    return isConstraint_;
}


void Foam::adjointSolver::clearSensitivities()
{
    sensitivities_.clear();
}


void Foam::adjointSolver::updatePrimalBasedQuantities()
{
    // Does nothing in base
}


// ************************************************************************* //
