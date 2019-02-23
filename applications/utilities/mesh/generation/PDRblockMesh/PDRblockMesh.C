/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    PDRblockMesh

Group
    grpMeshGenerationUtilities

Description
    A specialized single-block mesh generator for a rectilinear mesh
    in x-y-z.

    Uses the mesh description found in
      - \c system/PDRblockMeshDict

Usage
    \b PDRblockMesh [OPTION]

    Options:
      - \par -dict \<filename\>
        Alternative dictionary for the mesh description.

      - \par -noClean
        Do not remove any existing polyMesh/ directory or files

      - \par -time
        Write resulting mesh to a time directory (instead of constant)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "PDRblock.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OSspecific.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A block mesh generator for a rectilinear mesh in x-y-z.\n"
        "  The ordering of vertex and face labels within a block as shown "
        "below.\n"
        "  For the local vertex numbering in the sequence 0 to 7:\n"
        "    Faces 0, 1  ==  x-min, x-max.\n"
        "    Faces 2, 3  ==  y-min, y-max.\n"
        "    Faces 4, 5  ==  z-min, z-max.\n"
        "\n"
        "                        7 ---- 6\n"
        "                 f5     |\\     |\\     f3\n"
        "                 |      | 4 ---- 5     \\\n"
        "                 |      3 |--- 2 |      \\\n"
        "                 |       \\|     \\|      f2\n"
        "                 f4       0 ---- 1\n"
        "    Y  Z\n"
        "     \\ |                f0 ------ f1\n"
        "      \\|\n"
        "       O--- X\n"
    );

    argList::noParallel();
    argList::noFunctionObjects();

    argList::addBoolOption
    (
        "noClean",
        "Do not remove any existing polyMesh/ directory or files"
    );
    argList::addOption
    (
        "dict",
        "file",
        "Alternative dictionary for the PDRblockMesh description"
    );
    argList::addOption
    (
        "time",
        "time",
        "Specify a time to write mesh to (default: constant)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    // Instance for resulting mesh
    bool useTime = false;
    word meshInstance(runTime.constant());

    if
    (
        args.readIfPresent("time", meshInstance)
     && runTime.constant() != meshInstance
    )
    {
        // Verify that the value is actually good
        scalar timeValue;

        useTime = readScalar(meshInstance, timeValue);
        if (!useTime)
        {
            FatalErrorInFunction
                << "Bad input value: " << meshInstance
                << "Should be a scalar or 'constant'"
                << nl << endl
                << exit(FatalError);
        }
    }


    const word dictName("PDRblockMeshDict");

    #include "setSystemRunTimeDictionaryIO.H"

    IOdictionary meshDict(dictIO);

    Info<< "Creating PDRblockMesh from "
        << runTime.relativePath(dictIO.objectPath()) << endl;

    PDRblock blkMesh(meshDict, true);


    // Instance for resulting mesh
    if (useTime)
    {
        Info<< "Writing polyMesh to " << meshInstance << nl << endl;

        // Make sure that the time is seen to be the current time.
        // This is the logic inside regIOobject that resets the instance
        // to the current time before writing
        runTime.setTime(instant(meshInstance), 0);
    }

    if (!args.found("noClean"))
    {
        const fileName polyMeshPath
        (
            runTime.path()/meshInstance/polyMesh::meshSubDir
        );

        if (exists(polyMeshPath))
        {
            Info<< "Deleting polyMesh directory "
                << runTime.relativePath(polyMeshPath) << endl;
            rmDir(polyMeshPath);
        }
    }


    Info<< nl << "Creating polyMesh from PDRblockMesh" << endl;

    auto meshPtr = blkMesh.mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            meshInstance,
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );


    const polyMesh& mesh = *meshPtr;

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing polyMesh with "
        << mesh.cellZones().size() << " cellZones" << endl;

    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    // Mesh summary
    {
        Info<< "----------------" << nl
            << "Mesh Information" << nl
            << "----------------" << nl
            << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
            << "  " << "nPoints: " << mesh.nPoints() << nl
            << "  " << "nCells: " << mesh.nCells() << nl
            << "  " << "nFaces: " << mesh.nFaces() << nl
            << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

        Info<< "----------------" << nl
            << "Patches" << nl
            << "----------------" << nl;

        for (const polyPatch& p : mesh.boundaryMesh())
        {
            Info<< "  " << "patch " << p.index()
                << " (start: " << p.start()
                << " size: " << p.size()
                << ") name: " << p.name()
                << nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
