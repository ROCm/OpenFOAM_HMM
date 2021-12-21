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

#include "optimisationManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(optimisationManager, 0);
    defineRunTimeSelectionTable(optimisationManager, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optimisationManager::optimisationManager(fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "optimisationDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            true
        )
    ),
    mesh_(mesh),
    time_(const_cast<Time&>(mesh.time())),
    primalSolvers_(),
    adjointSolverManagers_(),
    managerType_(get<word>("optimisationManager")),
    optType_(nullptr)
{
    const dictionary& primalSolversDict = subDict("primalSolvers");
    const wordList& primalSolverNames = primalSolversDict.toc();

    // Construct primal solvers
    primalSolvers_.setSize(primalSolverNames.size());
    forAll(primalSolvers_, solveri)
    {
        primalSolvers_.set
        (
            solveri,
            primalSolver::New
            (
                mesh,
                managerType_,
                primalSolversDict.subDict(primalSolverNames[solveri])
            )
        );
    }

    // Construct adjointSolverManagers
    const dictionary& adjointManagersDict = subDict("adjointManagers");
    const wordList& adjointManagerNames = adjointManagersDict.toc();
    adjointSolverManagers_.setSize(adjointManagerNames.size());

    label nAdjointSolvers(0);
    forAll(adjointSolverManagers_, manageri)
    {
        adjointSolverManagers_.set
        (
            manageri,
            new adjointSolverManager
            (
                mesh,
                managerType_,
                adjointManagersDict.subDict(adjointManagerNames[manageri])
            )
        );
        nAdjointSolvers += adjointSolverManagers_[manageri].nAdjointSolvers();
    }

    // Sanity checks on the naming convention
    if (primalSolvers_.size() > 1)
    {
        for (const primalSolver& solveri : primalSolvers_)
        {
            if (!solveri.useSolverNameForFields())
            {
                FatalErrorInFunction
                    << "Multiple primal solvers are present but "
                    << "useSolverNameForFields is set to false in "
                    << "primal solver " << solveri.solverName() << nl
                    << "This is considered fatal."
                    << exit(FatalError);
            }
        }
    }

    if (nAdjointSolvers > 1)
    {
        for (const adjointSolverManager& amI : adjointSolverManagers_)
        {
            const PtrList<adjointSolver>& adjointSolvers = amI.adjointSolvers();
            for (const adjointSolver& asI : adjointSolvers)
            {
                if (!asI.useSolverNameForFields())
                {
                    FatalErrorInFunction
                        << "Multiple adjoint solvers are present but "
                        << "useSolverNameForFields is set to false in "
                        << "adjoint solver " << asI.solverName() << nl
                        << "This is considered fatal."
                        << exit(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::optimisationManager> Foam::optimisationManager::New
(
    fvMesh& mesh
)
{
    const IOdictionary dict
    (
        IOobject
        (
            "optimisationDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false // Do not register
        )
    );

    const word modelType(dict.get<word>("optimisationManager"));

    Info<< "optimisationManager type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "optimisationManager",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<optimisationManager>(ctorPtr(mesh));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

bool Foam::optimisationManager::read()
{
    if (regIOobject::read())
    {
        // Note: Only changing existing solvers - not adding any new
        const dictionary& primalSolversDict = subDict("primalSolvers");
        for (primalSolver& sol : primalSolvers_)
        {
            sol.readDict(primalSolversDict.subDict(sol.solverName()));
        }

        const dictionary& adjointManagersDict = subDict("adjointManagers");
        for (adjointSolverManager& man : adjointSolverManagers_)
        {
            man.readDict(adjointManagersDict.subDict(man.managerName()));
        }

        return true;
    }

    return false;
}


Foam::PtrList<Foam::primalSolver>& Foam::optimisationManager::primalSolvers()
{
    return primalSolvers_;
}


Foam::PtrList<Foam::adjointSolverManager>&
Foam::optimisationManager::adjointSolverManagers()
{
    return adjointSolverManagers_;
}


void Foam::optimisationManager::solvePrimalEquations()
{
    // Solve all primal equations
    forAll(primalSolvers_, psI)
    {
        primalSolvers_[psI].solve();
    }
}


void Foam::optimisationManager::solveAdjointEquations()
{
    // Solve all adjoint solver equations
    forAll(adjointSolverManagers_, amI)
    {
        adjointSolverManagers_[amI].solveAdjointEquations();
    }
}


void Foam::optimisationManager::computeSensitivities()
{
    // Compute senstivities from all adjoint solvers
    forAll(adjointSolverManagers_, amI)
    {
        adjointSolverManagers_[amI].computeAllSensitivities();
    }
}


void Foam::optimisationManager::updatePrimalBasedQuantities()
{
    forAll(adjointSolverManagers_, amI)
    {
        PtrList<adjointSolver>& adjointSolvers =
            adjointSolverManagers_[amI].adjointSolvers();

        forAll(adjointSolvers, asI)
        {
            adjointSolvers[asI].updatePrimalBasedQuantities();
        }
    }
}


// ************************************************************************* //
