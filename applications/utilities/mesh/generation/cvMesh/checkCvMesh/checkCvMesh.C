/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    checkCvMesh

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "autoSnapDriver.H"
#include "faceSet.H"
#include "motionSmoother.H"
#include "timeSelector.H"


using namespace Foam;


int main(int argc, char *argv[])
{
#   include "addOverwriteOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    runTime.functionObjects().off();

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Create mesh for time = " << runTime.timeName()
            << nl << endl;

        mesh.readUpdate();

        Info<< "Read mesh in = "
            << runTime.cpuTimeIncrement() << " s" << endl;

        // Check patches and faceZones are synchronised
        mesh.boundaryMesh().checkParallelSync(true);
        meshRefinement::checkCoupledFaceZones(mesh);

        // Read meshing dictionary
        IOdictionary cvMeshDict
        (
           IOobject
           (
                "cvMeshDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
           )
        );

        // mesh motion and mesh quality parameters
        const dictionary& meshQualityDict
            = cvMeshDict.subDict("meshQualityControls");


        Info<< "Checking initial mesh ..." << endl;
        faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh, meshQualityDict, wrongFaces);

        const label nInitErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        Info<< "Detected " << nInitErrors << " illegal faces"
            << " (concave, zero area or negative cell pyramid volume)"
            << endl;

        if (nInitErrors > 0)
        {
            Info<< "Writing " << nInitErrors
                << " faces in error to set "
                << wrongFaces.name() << endl;

            wrongFaces.instance() = mesh.pointsInstance();
            wrongFaces.write();
        }

        Info<< nl << "End of time " << runTime.timeName() << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;

}

