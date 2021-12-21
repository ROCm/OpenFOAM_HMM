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

#include "displacementMethodlaplacianMotionSolver.H"
#include "laplacianMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(displacementMethodlaplacianMotionSolver, 1);
addToRunTimeSelectionTable
(
    displacementMethod,
    displacementMethodlaplacianMotionSolver,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

displacementMethodlaplacianMotionSolver::displacementMethodlaplacianMotionSolver
(
    fvMesh& mesh,
    const labelList& patchIDs
)
:
    displacementMethod(mesh, patchIDs),
    pointMotionU_(refCast<laplacianMotionSolver>(motionPtr_()).pointMotionU()),
    cellMotionU_(refCast<laplacianMotionSolver>(motionPtr_()).cellMotionU()),
    /*
    // getting a reference from the database is dangerous since multiple
    // fields with the same name might be registered
    pointMotionU_
    (
        const_cast<pointVectorField&>
        (
            mesh_.lookupObject<pointVectorField>("pointMotionU")
        )
    ),
    cellMotionU_
    (
        const_cast<volVectorField&>
        (
            mesh_.lookupObject<volVectorField>("cellMotionU")
        )
    ),
    */
    resetFields_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::AUTO_WRITE,
                false
            )
        ).subDict("laplacianMotionSolverCoeffs").getOrDefault<bool>
        (
            "resetFields",
            true
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void displacementMethodlaplacianMotionSolver::setMotionField
(
    const pointVectorField& pointMovement
)
{
    if (resetFields_)
    {
        pointMotionU_.primitiveFieldRef() = vector::zero;
        cellMotionU_.primitiveFieldRef() = vector::zero;
        cellMotionU_.correctBoundaryConditions();
    }

    // Set boundary mesh movement and calculate
    // max current boundary displacement
    for (label patchI : patchIDs_)
    {
        // Set boundary field. Needed for the motionSolver class
        pointMotionU_.boundaryFieldRef()[patchI] ==
            pointMovement.boundaryField()[patchI].patchInternalField()();

        // Set boundary field values as seen from the internalField!
        // Needed for determining the max displacement
        pointMotionU_.boundaryFieldRef()[patchI].setInInternalField
        (
            pointMotionU_.primitiveFieldRef(),
            pointMovement.boundaryField()[patchI].patchInternalField()()
        );

        // Find max value
        maxDisplacement_ =
            max
            (
                maxDisplacement_,
                gMax
                (
                    mag
                    (
                        pointMotionU_.boundaryField()[patchI].
                            patchInternalField()
                    )
                )
            );
    }
}


void displacementMethodlaplacianMotionSolver::setMotionField
(
    const volVectorField& cellMovement
)
{
    NotImplemented;
    /*
    // Set boundary mesh movement and calculate max current boundary
    // displacement
    forAll(patchIDs_, pI)
    {
        label patchI = patchIDs_[pI];
        cellMotionU_.boundaryField()[patchI] ==
            cellMovement.boundaryField()[patchI];

        // Find max value
        maxDisplacement_ =
            max
            (
                maxDisplacement_,
                gMax
                (
                    mag(cellMovement.boundaryField()[patchI])
                )
            );
    }
    */
}


void displacementMethodlaplacianMotionSolver::setControlField
(
    const vectorField& controlField
)
{
    NotImplemented;
}


void displacementMethodlaplacianMotionSolver::setControlField
(
    const scalarField& controlField
)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
