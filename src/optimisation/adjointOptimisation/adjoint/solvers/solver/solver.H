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
    Foam::solver

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef solver_H
#define solver_H

#include "fvMesh.H"
#include "fvMatrix.H"
#include "localIOdictionary.H"
#include "variablesSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                            Class solver Declaration
\*---------------------------------------------------------------------------*/

class solver
:
    public localIOdictionary
{
protected:

    // Protected Data

        //- Reference to the mesh database
        fvMesh& mesh_;

        //- The optimisation type
        const word managerType_;

        //- Dictionary holding the solver info
        dictionary dict_;

        //- Solver name
        const word solverName_;

        //- Solve equations?
        bool active_;

        //- Pointer to a source term coming from the optimisationType
        //- (e.g. porosity from topologyOptimisation)
        //  Will never allocate new memory, so no need to be deleted
        //  in the destructor
        const volScalarField* optTypeSource_;

        //- Base variableSet pointer.
        //  To be allocated in derived classes
        autoPtr<variablesSet> vars_;


    // Protected Member Functions

        //- No copy construct
        solver(const solver&) = delete;

        //- No copy assignment
        void operator=(const solver&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solver");


    // Constructors

        //- Construct from mesh and dictionary
        solver
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict
        );


    //- Destructor
    virtual ~solver();


    // Member Functions

        virtual bool readDict(const dictionary& dict);


        // Access

            //- Return the solver mesh
            const fvMesh& mesh() const;

            //- Return the solver name
            const word& solverName() const;

            //- Use solver name as a suffix to the involved fields
            virtual bool useSolverNameForFields() const = 0;

            //- Return state of solver
            virtual bool active();

            //- Return the solver dictionary
            virtual const dictionary& dict() const;

            //- Return constant reference to variableSet used by the solver
            const variablesSet& getVariablesSet() const;

            //- Return non-constant reference to variableSet used by the solver
            variablesSet& getVariablesSet();


        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter() = 0;

            //- Main control loop
            virtual void solve() = 0;

            //- Looper (advances iters, time step)
            virtual bool loop() = 0;

            //- Restore initial field values if necessary
            virtual void restoreInitValues();

            //- Functions to be called before loop
            virtual void preLoop();

            //- Functions to be called after loop
            virtual void postLoop();

            //- Update source term related to optimisationType
            void updateOptTypeSource
            (
                const autoPtr<volScalarField>& optSourcePtr
            );

            //- Main control loop.
            //  Gets a list of function pointers to be called at the end of
            //  each solver iteration
            template<class Type>
            void solveWithArgs
            (
                Type& type,
                List<void (Type::*)()>& funcs
            );

            //- Add source from optimisationType to underlaying equation
            template<class Type>
            void addOptimisationTypeSource(fvMatrix<Type>& matrix) const;


        // IO

            //- Required by regIOobject
            virtual bool writeData(Ostream&) const
            {
                return true;
            }

            //- Workaround for turbulent fields on multi-point runs
            virtual bool write(const bool valid = true) const
            {
                return false;
            }

            //- Workaround for turbulent fields on multi-point runs
            virtual bool writeNow() const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
