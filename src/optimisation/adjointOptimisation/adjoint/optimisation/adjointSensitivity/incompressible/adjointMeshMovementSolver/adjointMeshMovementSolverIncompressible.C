/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "adjointMeshMovementSolverIncompressible.H"
#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointMeshMovementSolver, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void adjointMeshMovementSolver::read()
{
    nLaplaceIters_ = dict_.getOrDefault<label>("iters", 1000);
    tolerance_ = dict_.getOrDefault<scalar>("tolerance", 1e-6);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointMeshMovementSolver::adjointMeshMovementSolver
(
    const fvMesh& mesh,
    const dictionary& dict,
    Foam::incompressible::adjointSensitivity& adjointSensitivity,
    const labelHashSet& sensitivityPatchIDs,
    const autoPtr<adjointEikonalSolver>& adjointEikonalSolverPtr
)
:
    mesh_(mesh),
    dict_(dict.subOrEmptyDict("adjointMeshMovementSolver")),
    adjointSensitivity_(adjointSensitivity),
    sensitivityPatchIDs_(sensitivityPatchIDs),
    nLaplaceIters_(-1),
    tolerance_(-1),
    ma_
    (
        variablesSet::autoCreateMeshMovementField
        (
            mesh,
            "ma",
            dimensionSet(pow3(dimLength/dimTime))
        )
    ),
    source_
    (
        IOobject
        (
            "sourceAdjointMeshMovement",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimLength/pow3(dimTime), Zero)
    ),
    meshMovementSensPtr_(createZeroBoundaryPtr<vector>(mesh_)),
    adjointEikonalSolverPtr_(adjointEikonalSolverPtr)
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool adjointMeshMovementSolver::readDict(const dictionary& dict)
{
    dict_ = dict.subOrEmptyDict("adjointMeshMovementSolver");

    return true;
}


void adjointMeshMovementSolver::accumulateIntegrand(const scalar dt)
{
    // Accumulate integrand from the current time step
    source_ += adjointSensitivity_.adjointMeshMovementSource()*dt;

    // Part of the source depending on the adjoint distance can be added only
    // after solving the adjoint eikonal equation. Added in solve()
}


void adjointMeshMovementSolver::solve()
{
    read();

    // Add source from the adjoint eikonal equation
    if (adjointEikonalSolverPtr_)
    {
        source_ -=
            fvc::div(adjointEikonalSolverPtr_().getFISensitivityTerm()().T());
    }

    // Iterate the adjoint to the mesh movement equation
    for (label iter = 0; iter < nLaplaceIters_; iter++)
    {
        Info<< "Adjoint Mesh Movement Iteration: " << iter << endl;

        fvVectorMatrix maEqn
        (
            fvm::laplacian(ma_)
          + source_
        );

        maEqn.boundaryManipulate(ma_.boundaryFieldRef());

        //scalar residual = max(maEqn.solve().initialResidual());
        scalar residual = mag(maEqn.solve().initialResidual());

        Info<< "Max ma " << gMax(mag(ma_)()) << endl;

        mesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance_)
        {
            Info<< "\n***Reached adjoint mesh movement convergence limit, "
                   "iteration " << iter << "***\n\n";
            break;
        }
    }
    ma_.write();
}


void adjointMeshMovementSolver::reset()
{
    source_ == dimensionedVector(source_.dimensions(), Zero);
    meshMovementSensPtr_() = vector::zero;
}


boundaryVectorField& adjointMeshMovementSolver::meshMovementSensitivities()
{
    Info<< "Calculating mesh movement sensitivities " << endl;

    boundaryVectorField& meshMovementSens = meshMovementSensPtr_();

    for (const label patchi : sensitivityPatchIDs_)
    {
        // No surface area included. Will be done by the actual sensitivity tool
        meshMovementSens[patchi] = -ma_.boundaryField()[patchi].snGrad();
    }

    return meshMovementSens;
}


const volVectorField& adjointMeshMovementSolver::ma()
{
    return ma_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
