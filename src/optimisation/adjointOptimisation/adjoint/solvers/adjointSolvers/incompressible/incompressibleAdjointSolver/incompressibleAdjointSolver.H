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
    Foam::incompressibleAdjointSolver

Description
    Base class for incompressibleAdjoint solvers

\*---------------------------------------------------------------------------*/

#ifndef incompressibleAdjointSolver_H
#define incompressibleAdjointSolver_H

#include "adjointSolver.H"
#include "incompressibleVars.H"
#include "incompressibleAdjointVars.H"
#include "ATCModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class incompressibleAdjointSolver Declaration
\*---------------------------------------------------------------------------*/

class incompressibleAdjointSolver
:
    public adjointSolver
{
private:

    // Private Member Functions

        //- No copy construct
        incompressibleAdjointSolver
        (
            const incompressibleAdjointSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const incompressibleAdjointSolver&) = delete;


protected:

    // Protected data

        //- Primal variable set
        incompressibleVars& primalVars_;

        //- Adjoint Transpose Convection options
        autoPtr<ATCModel> ATCModel_;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("incompressible");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            incompressibleAdjointSolver,
            dictionary,
            (
                fvMesh& mesh,
                const word& managerType,
                const dictionary& dict,
                const word& primalSolverName
            ),
            (mesh, managerType, dict, primalSolverName)
        );


    // Constructors

        //- Construct from mesh and dictionary
        incompressibleAdjointSolver
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName
        );


    // Selectors

        //- Return a reference to the selected incompressible adjoint solver
        static autoPtr<incompressibleAdjointSolver> New
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~incompressibleAdjointSolver() = default;


    // Member Functions

        //- Read dict if updated
        virtual bool readDict(const dictionary& dict);

        //- Should solver name be appended to fields
        virtual bool useSolverNameForFields() const;


        // Access

            //- Access to the incompressible primal variables set
            const incompressibleVars& getPrimalVars() const;

            //- Access to the incompressible adjoint variables set
            virtual const incompressibleAdjointVars& getAdjointVars() const;

            //- Access to the incompressible adjoint variables set
            virtual incompressibleAdjointVars& getAdjointVars();

            //- Access to the ATC model
            const autoPtr<ATCModel>& getATCModel() const;

            //- Access to the ATC model
            autoPtr<ATCModel>& getATCModel();


        // Evolution

            //- Update primal based quantities, e.g. the primal fields
            //- in adjoint turbulence models
            virtual void updatePrimalBasedQuantities();


        // IO

            //- In case of multi-point runs with turbulent flows,
            //- output dummy turbulence fields with the base names, to allow
            //- continuation
            virtual bool write(const bool valid = true) const
            {
                if (mesh_.time().writeTime())
                {
                    return primalVars_.write();
                }

                return false;
            }

            //- In case of multi-point runs with turbulent flows,
            //- output dummy turbulence fields with the base names, to allow
            //- continuation
            virtual bool writeNow() const
            {
                return primalVars_.write();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
