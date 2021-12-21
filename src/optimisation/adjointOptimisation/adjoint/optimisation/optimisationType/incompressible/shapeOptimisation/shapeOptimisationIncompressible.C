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

#include "shapeOptimisationIncompressible.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(shapeOptimisation, 0);
addToRunTimeSelectionTable
(
    optimisationType,
    shapeOptimisation,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void shapeOptimisation::updateDesignVariables(scalarField& correction)
{
    // Communicate the movement to optMeshMovement
    optMeshMovement_->setCorrection(correction);

    if (updateGeometry_)
    {
        // Update the mesh
        optMeshMovement_->moveMesh();

        if (writeEachMesh_)
        {
            Info<< "  Writing new mesh points " << endl;
            pointIOField points
            (
                IOobject
                (
                   "points",
                    mesh_.pointsInstance(),
                    mesh_.meshSubDir,
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_.points()
            );
            points.write();
        }
    }
}


void shapeOptimisation::computeEta
(
    scalarField& correction
)
{
    if (!updateMethod_->initialEtaSet())
    {
        // In the unlikely event that eta is not set and the line search step
        // is not 1, multiply with it
        // if (lineSearch_) correction *= lineSearch_->step();

        // Compute eta based on desirable mesh movement size
        scalar eta = optMeshMovement_->computeEta(correction);
        correction *= eta;

        // Update eta known by the optimisation method and inform it that is
        // has been set
        updateMethod_->setStep(eta);
        updateMethod_->initialEtaSet() = true;

        // If a backtracking should be made at the first optimisation cycle,
        // the direction of the subsequent line searches of the same cycle
        // should also be scaled with the newly computed eta. We do this by
        // changing the line search step. This will happen only at the first
        // optimisation cycle since the updated value of eta will be included
        // in the line search direction in all subsequent optimisation cycles
        //correction *= eta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

shapeOptimisation::shapeOptimisation
(
    fvMesh& mesh,
    const dictionary& dict,
    PtrList<adjointSolverManager>& adjointSolverManagers
)
:
    optimisationType(mesh, dict, adjointSolverManagers),
    optMeshMovement_(nullptr),
    writeEachMesh_
    (
        dict.subDict("optimisationType").
            getOrDefault<bool>("writeEachMesh", false)
    ),
    updateGeometry_
    (
        dict.subDict("optimisationType").
            getOrDefault<bool>("updateGeometry", true)
    )
{
    // Note: to be updated
    labelHashSet patches
    (
        mesh_.boundaryMesh().patchSet
        (
            dict_.subDict("sensitivities").get<wordRes>("patches")
        )
    );
    if (patches.empty())
    {
        WarningInFunction
            << "There is no patch on which to compute sensitivities. "
            << "Check optimisationDict \n"
            << endl;
    }
    labelList sensitivityPatchIDs = patches.toc();
    optMeshMovement_.reset
    (
        optMeshMovement::New
        (
            mesh_,
            dict_.subDict("meshMovement"),
            sensitivityPatchIDs
        ).ptr()
    );

    // Sanity checks: at least one of eta or maxAllowedDisplacement must be set
    if
    (
        !updateMethod_->initialEtaSet()
     && !optMeshMovement_().maxAllowedDisplacementSet()
    )
    {
        FatalErrorInFunction
            << "Neither eta (updateMethod) "
            << "nor maxAllowedDisplacement (meshMovement) have been set"
            << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void shapeOptimisation::storeDesignVariables()
{
    optMeshMovement_->storeDesignVariables();
}


void shapeOptimisation::resetDesignVariables()
{
    optMeshMovement_->resetDesignVariables();
}


void shapeOptimisation::write()
{
    optimisationType::write();
    updateMethod_->writeCorrection();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
