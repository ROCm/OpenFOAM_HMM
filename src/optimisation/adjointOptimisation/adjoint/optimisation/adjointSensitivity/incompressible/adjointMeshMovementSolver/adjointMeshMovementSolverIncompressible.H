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
    Foam::incompressible::adjointMeshMovementSolver

Description
    Solver of the adjoint to the Laplace grid displacement equation

    Reference:
    \verbatim
        Kavvadias, I., Papoutsis-Kiachagias, E., & Giannakoglou, K. (2015).
        On the proper treatment of grid sensitivities in continuous adjoint
        methods for shape optimization.
        Journal of Computational Physics, 301, 1–18.
        http://doi.org/10.1016/j.jcp.2015.08.012
    \endverbatim

SourceFiles
    adjointMeshMovementSolver.C

\*---------------------------------------------------------------------------*/

#ifndef adjointMeshMovementSolverIncompressible_H
#define adjointMeshMovementSolverIncompressible_H

#include "adjointSensitivityIncompressible.H"
#include "adjointEikonalSolverIncompressible.H"
#include "createZeroField.H"
#include "boundaryFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                  Class adjointMeshMovementSolver Declaration
\*---------------------------------------------------------------------------*/

class adjointMeshMovementSolver
{
protected:

    // Protected data

        const fvMesh& mesh_;
        dictionary dict_;
        Foam::incompressible::adjointSensitivity& adjointSensitivity_;
        const labelHashSet& sensitivityPatchIDs_;
        label nLaplaceIters_;
        scalar tolerance_;
        volVectorField ma_;
        volVectorField source_;

        //- Wall face sens w.r.t.(x, y.z) //wall face sens w.r.t. (x,y.z)
        autoPtr<boundaryVectorField> meshMovementSensPtr_;
        const autoPtr<adjointEikonalSolver>& adjointEikonalSolverPtr_;

        //- Read options each time a new solution is found
        void read();


private:

    // Private Member Functions

        //- No copy construct
        adjointMeshMovementSolver(const adjointMeshMovementSolver&) = delete;

        //- No copy assignment
        void operator=(const adjointMeshMovementSolver&) = delete;


public:

    //- Runtime type information
    TypeName("adjointMeshMovementSolver");


    // Constructors

        //- Construct from components
        adjointMeshMovementSolver
        (
            const fvMesh& mesh,
            const dictionary& dict,
            Foam::incompressible::adjointSensitivity& adjointSensitivity,
            const labelHashSet& sensitivityPatchIDs,
            const autoPtr<adjointEikonalSolver>& adjointEikonalSolverPtr
        );

    //- Destructor
    virtual ~adjointMeshMovementSolver() = default;


    // Member Functions

       //- Read dict if changed
       virtual bool readDict(const dictionary& dict);

       //- Accumulate source term
       void accumulateIntegrand(const scalar dt);

       //- Calculate the adjoint distance field
       void solve();

       //- Reset source term
       void reset();

       //- Return the sensitivity term depending on da
       boundaryVectorField& meshMovementSensitivities();

       //- Return the adjoint distance field
       const volVectorField& ma();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
