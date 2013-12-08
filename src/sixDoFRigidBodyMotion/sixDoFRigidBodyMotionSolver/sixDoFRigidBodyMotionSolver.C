/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "sixDoFRigidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "patchDist.H"
#include "volPointInterpolation.H"
#include "uniformDimensionedFields.H"
#include "forces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFRigidBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        sixDoFRigidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionSolver::sixDoFRigidBodyMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    motion_(coeffDict()),
    patches_(wordReList(coeffDict().lookup("patches"))),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    distance_(readScalar(coeffDict().lookup("distance"))),
    rhoInf_(1.0),
    rhoName_(coeffDict().lookupOrDefault<word>("rhoName", "rho")),
    scale_
    (
        IOobject
        (
            "scale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    curTimeIndex_(-1)
{
    const fvMesh& fMesh = dynamic_cast<const fvMesh&>(mesh);

    if (rhoName_ == "rhoInf")
    {
        rhoInf_ = readScalar(dict.lookup("rhoInf"));
    }

    // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        patchDist dist(fMesh, patchSet_);

        // Update distance to object
        pointScalarField pDist
        (
            volPointInterpolation::New
            (
                fMesh
            ).interpolate(dist)
        );

        // Scaling. 1 at patches, 0 at distance_ away from patches
        scale_.internalField() = 1 - min(1.0, pDist.internalField()/distance_);
        scale_.correctBoundaryConditions();
        scale_.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionSolver::~sixDoFRigidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyMotionSolver::curPoints() const
{
    return points0() + pointDisplacement_.internalField();
}


void Foam::sixDoFRigidBodyMotionSolver::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorIn
        (
            "sixDoFRigidBodyMotionSolver::curPoints() const"
        )   << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    // Patch force data is valid for the current positions, so
    // calculate the forces on the motion object from this data, then
    // update the positions

    motion_.updatePosition(t.deltaTValue(), t.deltaT0Value());

    dictionary forcesDict;

    forcesDict.add("type", forces::typeName);
    forcesDict.add("patches", patches_);
    forcesDict.add("rhoInf", rhoInf_);
    forcesDict.add("rhoName", rhoName_);
    forcesDict.add("CofR", motion_.centreOfMass());

    forces f("forces", db(), forcesDict);

    f.calcForcesMoment();

    uniformDimensionedVectorField g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    // scalar ramp = min(max((this->db().time().value() - 5)/10, 0), 1);
    scalar ramp = 1.0;

    motion_.updateAcceleration
    (
        ramp*(f.forceEff() + g.value()*motion_.mass()),
        ramp*(f.momentEff()),
        t.deltaTValue()
    );

    // Update the displacements
    pointDisplacement_.internalField() =
        motion_.scaledPosition(points0(), scale_) - points0();

    pointDisplacement_.correctBoundaryConditions();
}


// ************************************************************************* //
