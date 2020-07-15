/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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

#include "steadyOptimisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyOptimisation, 0);
    addToRunTimeSelectionTable
    (
        optimisationManager,
        steadyOptimisation,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::steadyOptimisation::updateOptTypeSource()
{
    forAll(primalSolvers_, pI)
    {
        primalSolvers_[pI].updateOptTypeSource(optType_->sourcePtr());
    }

    forAll(adjointSolverManagers_, asmI)
    {
        PtrList<adjointSolver>& adjointSolvers =
            adjointSolverManagers_[asmI].adjointSolvers();

        forAll(adjointSolvers, aI)
        {
            adjointSolvers[aI].updateOptTypeSource(optType_->sourcePtr());
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::steadyOptimisation::lineSearchUpdate()
{
    // Compute direction of update
    tmp<scalarField> tdirection = optType_->computeDirection();
    scalarField& direction = tdirection.ref();

    // Grab reference to line search
    autoPtr<lineSearch>& lineSrch = optType_->getLineSearch();

    // Store starting point
    optType_->storeDesignVariables();

    // Compute merit function before update
    scalar meritFunction = optType_->computeMeritFunction();
    lineSrch->setOldMeritValue(meritFunction);

    // Get merit function derivative
    const scalar dirDerivative =
        optType_->meritFunctionDirectionalDerivative();
    lineSrch->setDeriv(dirDerivative);
    lineSrch->setDirection(direction);

    // Reset initial step.
    // Might be interpolated from previous optimisation cycles
    lineSrch->reset();

    // Perform line search
    for (label iter = 0; iter < lineSrch->maxIters(); ++iter)
    {
        Info<< "\n- - - - - - - - - - - - - - -"  << endl;
        Info<< "Line search iteration "   << iter << endl;
        Info<< "- - - - - - - - - - - - - - -\n"  << endl;

        // Update design variables. Multiplication with line search step
        // happens inside the update(direction) function
        optType_->update(direction);

        // Solve all primal equations
        solvePrimalEquations();

        // Compute and set new merit function
        meritFunction = optType_->computeMeritFunction();
        lineSrch->setNewMeritValue(meritFunction);

        if (lineSrch->converged())
        {
            // If line search criteria have been met, proceed
            Info<< "Line search converged in " << iter + 1
                << " iterations." << endl;
            scalarField scaledCorrection(lineSrch->step()*direction);
            optType_->updateOldCorrection(scaledCorrection);
            optType_->write();
            lineSrch()++;
            break;
        }
        else
        {
            // If maximum number of iteration has been reached, continue
            if (iter == lineSrch->maxIters() - 1)
            {
                Info<< "Line search reached max. number of iterations.\n"
                    << "Proceeding to the next optimisation cycle" << endl;
                scalarField scaledCorrection(lineSrch->step()*direction);
                optType_->updateOldCorrection(scaledCorrection);
                optType_->write();
                lineSrch()++;
            }
            // Reset to initial design variables and update step
            else
            {
                optType_->resetDesignVariables();
                lineSrch->updateStep();
            }
        }
    }
}


void Foam::steadyOptimisation::fixedStepUpdate()
{
    // Update design variables
    optType_->update();

    // Solve primal equations
    solvePrimalEquations();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::steadyOptimisation::steadyOptimisation(fvMesh& mesh)
:
    optimisationManager(mesh)
{
    optType_.reset
    (
        incompressible::optimisationType::New
        (
            mesh,
            subDict("optimisation"),
            adjointSolverManagers_
        ).ptr()
    );

    // Update source ptrs in all solvers to look at the source held in optType
    // Possible problem if mesh is adapted
    updateOptTypeSource();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::optimisationManager& Foam::steadyOptimisation::operator++()
{
    time_++;
    if (!end())
    {
        Info<< "\n* * * * * * * * * * * * * * * * *" << endl;
        Info<< "Optimisation cycle " << time_.value() << endl;
        Info<< "* * * * * * * * * * * * * * * * *\n" << endl;
    }
    return *this;
}


Foam::optimisationManager& Foam::steadyOptimisation::operator++(int)
{
    return operator++();
}


bool Foam::steadyOptimisation::checkEndOfLoopAndUpdate()
{
    if (update())
    {
        optType_->update();
    }
    return end();
}


bool Foam::steadyOptimisation::end()
{
    return time_.end();
}


bool Foam::steadyOptimisation::update()
{
    return (time_.timeIndex() != 1 && !end());
}


void Foam::steadyOptimisation::updateDesignVariables()
{
    // Update design variables using either a line-search scheme or
    // a fixed-step update
    if (optType_->getLineSearch())
    {
        lineSearchUpdate();
    }
    else
    {
        fixedStepUpdate();
    }

    // Reset adjoint sensitivities in all adjoint solver managers
    for (adjointSolverManager& adjSolverManager : adjointSolverManagers_)
    {
        adjSolverManager.clearSensitivities();
    }
}


// ************************************************************************* //
