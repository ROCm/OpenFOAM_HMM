/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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
    Foam::incompressible::optimisationType

Description
    Abstract base class for optimisation methods

SourceFiles
    optimisationType.C

\*---------------------------------------------------------------------------*/

#ifndef optimisationTypeIncompressible_H
#define optimisationTypeIncompressible_H

#include "adjointSolverManager.H"
#include "updateMethod.H"
#include "lineSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                      Class optimisationType Declaration
\*---------------------------------------------------------------------------*/

class optimisationType
{
protected:

    // Protected data

        fvMesh& mesh_;
        const dictionary dict_;
        PtrList<adjointSolverManager>& adjointSolvManagers_;
        autoPtr<updateMethod> updateMethod_;
        autoPtr<volScalarField> sourcePtr_;
        autoPtr<lineSearch> lineSearch_;

        //- Update the design variables given their correction
        virtual void updateDesignVariables(scalarField& correction) = 0;

        //- Compute eta if not set in the first step
        virtual void computeEta(scalarField& correction) = 0;


private:

    // Private Member Functions

        //- No copy construct
        optimisationType(const optimisationType&) = delete;

        //- No copy assignment
        void operator=(const optimisationType&) = delete;


public:

    //- Runtime type information
    TypeName("optimisationType");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            optimisationType,
            dictionary,
            (
                fvMesh& mesh,
                const dictionary& dict,
                PtrList<adjointSolverManager>& adjointSolverManagers
            ),
            (mesh, dict, adjointSolverManagers)
        );



    // Constructors

        //- Construct from components
        optimisationType
        (
            fvMesh& mesh,
            const dictionary& dict,
            PtrList<adjointSolverManager>& adjointSolverManagers
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<optimisationType> New
        (
            fvMesh& mesh,
            const dictionary& dict,
            PtrList<adjointSolverManager>& adjointSolverManagers
        );


    // Destructor

        virtual ~optimisationType() = default;

        //- Update design variables
        virtual void update();

        //- Update design variables based on a given correction
        virtual void update(scalarField& correction);

        //- Store design variables, as the starting point for line search
        virtual void storeDesignVariables() = 0;

        //- Reset to starting point of line search
        virtual void resetDesignVariables() = 0;

        //- Compute update direction
        virtual tmp<scalarField> computeDirection();

        //- Compute cumulative objective and constraint gradients
        virtual void updateGradientsAndValues
        (
            scalarField& objectiveSens,
            PtrList<scalarField>& constraintSens,
            scalar& objectiveValue,
            scalarField& constraintValues
        );

        //- Compute the merit function of the optimisation problem.
        //  Could be different than the objective function in case of
        //  constraint optimisation
        virtual scalar computeMeritFunction();

        //- Derivative of the merit function
        virtual scalar meritFunctionDirectionalDerivative();

        //- Update old correction. Needed for quasi-Newton Methods
        virtual void updateOldCorrection(const scalarField&);

        //- Write useful quantities to files
        virtual void write();

        //- Get source term
        const autoPtr<volScalarField>& sourcePtr();

        //- Get a reference to the line search object
        autoPtr<lineSearch>& getLineSearch();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
