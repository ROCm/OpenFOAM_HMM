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


Class
    Foam::optimisationManager

Description
    Abstract base class for optimisation methods

SourceFiles
    optimisationManager.C

\*---------------------------------------------------------------------------*/

#ifndef optimisationManager_H
#define optimisationManager_H

#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "optimisationTypeIncompressible.H"
#include "primalSolver.H"
#include "adjointSolverManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class optimisationManager Declaration
\*---------------------------------------------------------------------------*/

class optimisationManager
:
    public IOdictionary
{
protected:

    // Protected data

        fvMesh& mesh_;
        Time& time_;
        PtrList<primalSolver> primalSolvers_;
        PtrList<adjointSolverManager> adjointSolverManagers_;
        const word managerType_;
        autoPtr<incompressible::optimisationType> optType_;


private:

    // Private Member Functions

        //- No copy construct
        optimisationManager(const optimisationManager&) = delete;

        //- No copy assignment
        void operator=(const optimisationManager&) = delete;


public:

    //- Runtime type information
    TypeName("optimisationManager");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            optimisationManager,
            dictionary,
            (
                fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        optimisationManager(fvMesh& mesh);


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<optimisationManager> New(fvMesh& mesh);


    //- Destructor
    virtual ~optimisationManager() = default;


    // Member Functions

        virtual PtrList<primalSolver>& primalSolvers();

        virtual PtrList<adjointSolverManager>& adjointSolverManagers();

        virtual bool read();

        //- Prefix increment,
        virtual optimisationManager& operator++() = 0;

        //- Postfix increment, this is identical to the prefix increment
        virtual optimisationManager& operator++(int) = 0;

        //- Return true if end of optimisation run.
        //  Also, updates the design variables if needed
        virtual bool checkEndOfLoopAndUpdate() = 0;

        //- Return true if end of optimisation run
        virtual bool end() = 0;

        //- Whether to update the design variables
        virtual bool update() = 0;

        //- Update design variables.
        //  Might employ a line search to find a correction satisfying the step
        //  convergence criteria
        virtual void updateDesignVariables() = 0;

        //- Solve all primal equations
        virtual void solvePrimalEquations();

        //- Solve all adjoint equations
        virtual void solveAdjointEquations();

        //- Compute all adjoint sensitivities
        virtual void computeSensitivities();

        //- Solve all primal equations
        virtual void updatePrimalBasedQuantities();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
