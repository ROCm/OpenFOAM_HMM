/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd
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

#include "addToRunTimeSelectionTable.H"
#include "dynamicMotionSolverTopoFvMesh.H"
#include "mapPolyMesh.H"
#include "OBJstream.H"
#include "Time.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverTopoFvMesh, 0);

    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverTopoFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverTopoFvMesh::dynamicMotionSolverTopoFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
    motionPtr_(motionSolver::New(*this))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverTopoFvMesh::~dynamicMotionSolverTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicMotionSolverTopoFvMesh::update()
{
    // Do mesh changes (not using inflation - points added directly into mesh)
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh(false);

    if (topoChangeMap.valid())
    {
        Info << "Executing mesh topology update" << endl;
        motionPtr_->updateMesh(topoChangeMap());

        setV0() = V();

        pointField newPoints(motionPtr_->newPoints());
        movePoints(newPoints);

        if (debug)
        {
            OBJstream osOld("oldPts_" + time().timeName() + ".obj");
            const pointField& oldPts = oldPoints();
            forAll(oldPts, i)
            {
                osOld.write(oldPts[i]);
            }

            OBJstream osNew("newPts_" + time().timeName() + ".obj");
            forAll(points(), i)
            {
                osNew.write(points()[i]);
            }
        }
    }
    else
    {
        // Calculate the new point positions using the motion solver
        pointField newPoints(motionPtr_->newPoints());

        // The mesh now contains the cells with zero volume
        Info << "Executing mesh motion" << endl;
        movePoints(newPoints);
    }

    return true;
}


// ************************************************************************* //
