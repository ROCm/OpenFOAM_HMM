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

\*---------------------------------------------------------------------------*/

#include "optMeshMovementVolumetricBSplines.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(optMeshMovementVolumetricBSplines, 0);
    addToRunTimeSelectionTable
    (
        optMeshMovement,
        optMeshMovementVolumetricBSplines,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vectorField Foam::optMeshMovementVolumetricBSplines::controlPointMovement
(
    const scalarField& correction
)
{
    const label nControlPoints(correction.size()/3);
    vectorField cpMovement(nControlPoints, Zero);

    for (label iCP = 0; iCP < nControlPoints; ++iCP)
    {
        cpMovement[iCP].x() = correction[3*iCP];
        cpMovement[iCP].y() = correction[3*iCP + 1];
        cpMovement[iCP].z() = correction[3*iCP + 2];
    }
    displMethodPtr_->boundControlField(cpMovement);

    return cpMovement;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optMeshMovementVolumetricBSplines::optMeshMovementVolumetricBSplines
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
    cpsInit_(volBSplinesBase_.getNumberOfBoxes())
{
    PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();

    forAll(boxes, boxI)
    {
        cpsInit_[boxI].setSize
        (
            boxes[boxI].getControlPoints().size()
        );
        cpsInit_[boxI] = boxes[boxI].getControlPoints();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optMeshMovementVolumetricBSplines::moveMesh()
{
    // Get controlPoint movement from correction
    vectorField cpMovement = controlPointMovement(correction_);

    // Set movement of the B-Splines control points
    displMethodPtr_->setControlField(cpMovement);

    // Move the mesh and check quality
    optMeshMovement::moveMesh();
}


void Foam::optMeshMovementVolumetricBSplines::storeDesignVariables()
{
    optMeshMovement::storeDesignVariables();
    const PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxes();
    forAll(boxes, boxI)
    {
        cpsInit_[boxI] = boxes[boxI].getControlPoints();
    }
}


void Foam::optMeshMovementVolumetricBSplines::resetDesignVariables()
{
    // Reset mesh points
    optMeshMovement::resetDesignVariables();

    DebugInfo
        << "optMeshMovementVolumetricBSplines:: resetting control points"
        << endl;

    PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
    forAll(boxes, boxI)
    {
        boxes[boxI].setControlPoints(cpsInit_[boxI]);
    }
}


Foam::scalar Foam::optMeshMovementVolumetricBSplines::computeEta
(
    const scalarField& correction
)
{
    const vectorField cpMovement(controlPointMovement(correction));
    const scalar maxDisplacement
    (
        volBSplinesBase_.computeMaxBoundaryDisplacement
        (
            cpMovement,
            patchIDs_
        )
    );

    Info<< "maxAllowedDisplacement/maxDisplacement of boundary\t"
        << getMaxAllowedDisplacement() << "/" << maxDisplacement << endl;
    scalar eta = getMaxAllowedDisplacement() / maxDisplacement;
    Info<< "Setting eta value to " << eta << endl;

    return eta;
}


Foam::labelList
Foam::optMeshMovementVolumetricBSplines::getActiveDesignVariables() const
{
    return volBSplinesBase_.getActiveDesignVariables();
}


// ************************************************************************* //
