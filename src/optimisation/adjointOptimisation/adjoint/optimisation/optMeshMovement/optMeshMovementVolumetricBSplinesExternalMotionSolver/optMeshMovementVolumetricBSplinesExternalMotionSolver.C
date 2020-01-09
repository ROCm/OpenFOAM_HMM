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

\*---------------------------------------------------------------------------*/

#include "optMeshMovementVolumetricBSplinesExternalMotionSolver.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        optMeshMovementVolumetricBSplinesExternalMotionSolver,
        0
    );
    addToRunTimeSelectionTable
    (
        optMeshMovement,
        optMeshMovementVolumetricBSplinesExternalMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::optMeshMovementVolumetricBSplinesExternalMotionSolver::
computeBoundaryMovement
(
    const scalarField& correction
)
{
    const label nCPs(volBSplinesBase_.getTotalControlPointsNumber());
    dx_.primitiveFieldRef() = vector::zero;
    cpMovement_ = vector::zero;

    for (label iCP = 0; iCP < nCPs; iCP++)
    {
        cpMovement_[iCP].x() = correction[3*iCP];
        cpMovement_[iCP].y() = correction[3*iCP + 1];
        cpMovement_[iCP].z() = correction[3*iCP + 2];
    }

    // Bound control point movement for non-active CPs
    volBSplinesBase_.boundControlPointMovement(cpMovement_);

    // Compute boundary movement
    label passedCPs(0);
    PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
    forAll(boxes, iNURB)
    {
        const label nb = boxes[iNURB].getControlPoints().size();
        for (label cpI = 0; cpI < nb; ++cpI)
        {
            const label globalCP = passedCPs + cpI;
            forAll(patchIDs_, pI)
            {
                const label patchI = patchIDs_[pI];
                vectorField boundaryMovement
                (
                    boxes[iNURB].patchDxDb(patchI, cpI)
                  & cpMovement_[globalCP]
                );
                dx_.boundaryField()[patchI].addToInternalField
                (
                    dx_.primitiveFieldRef(),
                    boundaryMovement
                );
            }
        }

        // Increment number of passed sensitivities
        passedCPs += nb;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optMeshMovementVolumetricBSplinesExternalMotionSolver::
optMeshMovementVolumetricBSplinesExternalMotionSolver
(
    fvMesh& mesh,
    const dictionary& dict,
    const labelList& patchIDs
)
:
    optMeshMovement(mesh, dict, patchIDs),
    volBSplinesBase_
    (
        const_cast<volBSplinesBase&>(volBSplinesBase::New(mesh))
    ),
    dx_
    (
        IOobject
        (
            "dx",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedVector(dimless, Zero)
    ),
    cpMovement_(volBSplinesBase_.getTotalControlPointsNumber(), Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optMeshMovementVolumetricBSplinesExternalMotionSolver::moveMesh()
{
    // Compute boundary movement
    computeBoundaryMovement(correction_);

    // Set boundary movement of motion solver
    displMethodPtr_->setMotionField(dx_);

    // Positions of control points have not changed since only the boundary dx
    // has been computed.
    // Use correction to update them
    volBSplinesBase_.moveControlPoints(cpMovement_);

    // Write control points to files
    volBSplinesBase_.writeControlPoints();

    // Move the mesh and check quality
    optMeshMovement::moveMesh();
}


Foam::scalar
Foam::optMeshMovementVolumetricBSplinesExternalMotionSolver::computeEta
(
    const scalarField& correction
)
{
    computeBoundaryMovement(correction);

    // Get maximum boundary movement
    scalar maxDisplacement = gMax(mag(dx_.primitiveField()));

    // Compute eta value
    Info<< "maxAllowedDisplacement/maxDisplacement \t"
        << getMaxAllowedDisplacement() << "/" << maxDisplacement << endl;
    scalar eta = getMaxAllowedDisplacement() / maxDisplacement;
    Info<< "Setting eta value to " << eta << endl;

    return eta;
}


Foam::labelList
Foam::optMeshMovementVolumetricBSplinesExternalMotionSolver::getActiveDesignVariables()
const
{
    return volBSplinesBase_.getActiveDesignVariables();
}


// ************************************************************************* //
