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



Description
    Shape optimisation support library

Class
    Foam::incompressible::shapeOptimisation

Description
    Calculates shape sensitivities using the adjoint approach,
    computes boundaryMesh movement and propagates it to the volume mesh

SourceFiles
    shapeOptimisation.C

\*---------------------------------------------------------------------------*/

#ifndef shapeOptimisationIncompressible_H
#define shapeOptimisationIncompressible_H

#include "optimisationTypeIncompressible.H"
#include "optMeshMovement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                      Class shapeOptimisation Declaration
\*---------------------------------------------------------------------------*/

class shapeOptimisation
:
    public optimisationType
{
protected:

    // Protected data

        autoPtr<optMeshMovement> optMeshMovement_;

        bool writeEachMesh_;
        bool updateGeometry_;


    // Protected Member Functions

        //- Update the design variables given their correction
        virtual void updateDesignVariables(scalarField& correction);

        //- Compute eta if not set in the first step
        virtual void computeEta(scalarField& correction);


private:

    // Private Member Functions

        //- No copy construct
        shapeOptimisation(const shapeOptimisation&) = delete;

        //- No copy assignment
        void operator=(const shapeOptimisation&) = delete;


public:

    //- Runtime type information
    TypeName("shapeOptimisation");


    // Constructors

        //- Construct from components
        shapeOptimisation
        (
            fvMesh& mesh,
            const dictionary& dict,
            PtrList<adjointSolverManager>& adjointSolverManagers
        );


    //- Destructor
    virtual ~shapeOptimisation() = default;


    // Member Functions

        //- Store design variables, as the starting point for line search
        virtual void storeDesignVariables();

        //- Store design variables, as the starting point for line search
        virtual void resetDesignVariables();

        //- Write useful quantities to files
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
