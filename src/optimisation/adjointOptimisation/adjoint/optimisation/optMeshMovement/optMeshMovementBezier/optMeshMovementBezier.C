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

#include "optMeshMovementBezier.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(optMeshMovementBezier, 0);
    addToRunTimeSelectionTable
    (
        optMeshMovement,
        optMeshMovementBezier,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::optMeshMovementBezier::computeBoundaryMovement
(
    const scalarField& correction
)
{
    // Re-initialize movement to zero
    dx_.primitiveFieldRef() = vector::zero;

    // Compute boundary mesh movement using derivatives of the control points
    // and parameterization information
    const label nBezier = Bezier_.nBezier();
    const boolList& confineXmovement = Bezier_.confineXmovement();
    const boolList& confineYmovement = Bezier_.confineYmovement();
    const boolList& confineZmovement = Bezier_.confineZmovement();
    vectorField actualMovement(nBezier, Zero);
    for (label iCP = 0; iCP < nBezier; iCP++)
    {
        // Confine x movement
        if (!confineXmovement[iCP])
        {
            actualMovement[iCP].x() = correction[iCP];
        }
        // Confine y movement
        if (!confineYmovement[iCP])
        {
            actualMovement[iCP].y() = correction[iCP + nBezier];
        }
        // Confine z movement
        if (!confineZmovement[iCP])
        {
            actualMovement[iCP].z() = correction[iCP + 2*nBezier];
        }
        dx_ += Bezier_.dxidXj()[iCP] & actualMovement[iCP];
    }

    // Add to cumulative control point change (wrong in the first optimisation
    // cycle if initial eta not set)
    cumulativeChange_ += actualMovement;
    Info<< "Cumulative control point change " << cumulativeChange_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optMeshMovementBezier::optMeshMovementBezier
(
    fvMesh& mesh,
    const dictionary& dict,
    const labelList& patchIDs
)
:
    optMeshMovement(mesh, dict, patchIDs),
    Bezier_(mesh, mesh.lookupObject<IOdictionary>("optimisationDict")),
    dx_
    (
        IOobject
        (
            "dx",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedVector(dimless, Zero)
    ),
    cumulativeChange_(Bezier_.nBezier(), Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optMeshMovementBezier::moveMesh()
{
    // Update the boundary movement
    computeBoundaryMovement(correction_);

    // Set boundary movement of motion solver
    displMethodPtr_->setMotionField(dx_);

    // Move the mesh and check quality
    optMeshMovement::moveMesh();
}


Foam::scalar
Foam::optMeshMovementBezier::computeEta(const scalarField& correction)
{
    // Set unscaled correction
    computeBoundaryMovement(correction);

    // Get maximum boundary movement
    const scalar maxDisplacement = gMax(mag(dx_.primitiveField()));

    // Compute eta value
    Info<< "maxAllowedDisplacement/maxDisplacement \t"
        << getMaxAllowedDisplacement() << "/" << maxDisplacement << endl;
    const scalar eta = getMaxAllowedDisplacement()/maxDisplacement;
    Info<< "Setting eta value to " << eta << endl;

    return eta;
}


Foam::labelList Foam::optMeshMovementBezier::getActiveDesignVariables() const
{
    return Bezier_.getActiveDesignVariables();
}


// ************************************************************************* //
