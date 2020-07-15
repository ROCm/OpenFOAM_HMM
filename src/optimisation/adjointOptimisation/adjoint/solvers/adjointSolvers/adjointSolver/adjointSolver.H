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
    Foam::adjointSolver

Description
    Base class for adjoint solvers

\*---------------------------------------------------------------------------*/

#ifndef adjointSolver_H
#define adjointSolver_H

#include "fvMesh.H"
#include "Time.H"
#include "IOdictionary.H"
#include "solver.H"
#include "objectiveManager.H"
#include "sensitivity.H"
#include "primalSolver.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class adjointSolver Declaration
\*---------------------------------------------------------------------------*/

class adjointSolver
:
    public solver
{
private:

    // Private Member Functions

        //- No copy construct
        adjointSolver(const adjointSolver&) = delete;

        //- No copy assignment
        void operator=(const adjointSolver&) = delete;


protected:

    // Protected data

        //- Name of primal solver
        const word primalSolverName_;

        //- Object to manage objective functions
        autoPtr<objectiveManager> objectiveManagerPtr_;

        //- Sensitivities field
        autoPtr<scalarField> sensitivities_;

        //- Are sensitivities computed
        bool computeSensitivities_;

        //- Is the adjoint solver used to tackle a constraint
        bool isConstraint_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("adjointSolver");


    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            adjointSolver,
            adjointSolver,
            (
                fvMesh& mesh,
                const word& managerType,
                const dictionary& dict,
                const word& primalSolverName
            ),
            (mesh, managerType, dict, primalSolverName)
        );


    // Constructors

        //- Construct from mesh, dictionary, and primal solver name
        adjointSolver
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<adjointSolver> New
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~adjointSolver() = default;


    // Member Functions

        // Access

            virtual bool readDict(const dictionary& dict);

            //- Return the primal solver name
            const word& primalSolverName() const
            {
                return primalSolverName_;
            }

            //- Return a const-reference to the primal solver
            //- corresponding to this adjoint solver
            const primalSolver& getPrimalSolver() const;

            //- Return a non const-reference to the primal solver
            //- corresponding to this adjoint solver
            primalSolver& getPrimalSolver();

            //- Return a const reference to the objective manager
            const objectiveManager& getObjectiveManager() const;

            //- Return a reference to the objective manager
            objectiveManager& getObjectiveManager();


        // Evolution

            //- Compute sensitivities of the underlaying objectives
            virtual void computeObjectiveSensitivities() = 0;

            //- Is the solving referring to a constraint
            virtual bool isConstraint();

            //- Grab a reference to the computed sensitivities
            virtual const scalarField& getObjectiveSensitivities() = 0;

            //- Clears the sensitivity field known by the adjoint solver
            virtual void clearSensitivities();

            //- Return the base sensitivity object
            virtual sensitivity& getSensitivityBase() = 0;

            //- Update primal based quantities, e.g. the primal fields
            //- in adjoint turbulence models
            //  Does nothing in the base
            virtual void updatePrimalBasedQuantities();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
