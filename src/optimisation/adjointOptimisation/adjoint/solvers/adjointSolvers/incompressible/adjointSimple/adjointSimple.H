/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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
    Foam::adjointSimple

Description
    Solution of the adjoint PDEs for incompressible, steady-state flows

    Reference:
    \verbatim
        For the development of the adjoint PDEs

            Papoutsis-Kiachagias, E. M., & Giannakoglou, K. C. (2014).
            Continuous Adjoint Methods for Turbulent Flows, Applied to Shape
            and Topology Optimization: Industrial Applications.
            Archives of Computational Methods in Engineering, 23(2), 255-299.
            http://doi.org/10.1007/s11831-014-9141-9
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef adjointSimple_H
#define adjointSimple_H

#include "incompressibleAdjointSolver.H"
#include "SIMPLEControl.H"
#include "incompressibleVars.H"
#include "incompressibleAdjointVars.H"
#include "adjointSensitivityIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class adjointSimple Declaration
\*---------------------------------------------------------------------------*/

class adjointSimple
:
    public incompressibleAdjointSolver
{

private:

    // Private Member Functions

        //- No copy construct
        adjointSimple(const adjointSimple&) = delete;

        //- No copy assignment
        void operator=(const adjointSimple&) = delete;


protected:

    // Protected data

        //- Solver control
        autoPtr<SIMPLEControl> solverControl_;

        //- Cumulative continuity error
        scalar cumulativeContErr_;

        //- Sensitivity Derivatives engine
        autoPtr<incompressible::adjointSensitivity> adjointSensitivity_;


    // Protected Member Functions

        //- In case variable names are different than the base ones,
        //- add extra schemes and relaxation factors to the appropriate dicts
        //  Note: 160813: Changes in the baseline solution and fvSchemes classes
        //  have to be made in order to add schemes in the dict at run time.
        //  Not supported for now
        void addExtraSchemes();

        //- Compute continuity errors
        void continuityErrors();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("adjointSimple");


    // Constructors

        //- Construct from mesh and dictionary
        adjointSimple
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~adjointSimple() = default;


    // Member Functions

        //- Read dict if updated
        virtual bool readDict(const dictionary& dict);


        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter();

            //- Main control loop
            virtual void solve();

            //- Looper (advances iters, time step)
            virtual bool loop();

            //- Compute sensitivities of the underlaying objectives
            virtual void computeObjectiveSensitivities();

            //- Grab a reference to the computed sensitivities
            virtual const scalarField& getObjectiveSensitivities();

            //- Return the base sensitivity object
            virtual sensitivity& getSensitivityBase();

            //- Write average iteration
            virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
