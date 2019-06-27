/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011, 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    blockMesh

Group
    grpMeshGenerationUtilities

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
      - \c system/blockMeshDict
      - \c system/\<region\>/blockMeshDict
      - \c constant/polyMesh/blockMeshDict
      - \c constant/\<region\>/polyMesh/blockMeshDict

Usage
    \b blockMesh [OPTION]

    Options:
      - \par -blockTopology
        Write the topology as a set of edges in OBJ format and exit.

      - \par -region \<name\>
        Specify alternative mesh region.

      - \par -dict \<filename\>
        Alternative dictionary for the block mesh description.

      - \par -sets
        Write cellZones as cellSets too (for processing purposes)

      - \par -noClean
        Do not remove any existing polyMesh/ directory or files

      - \par -time
        Write resulting mesh to a time directory (instead of constant)

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "attachPolyTopoChanger.H"
#include "polyTopoChange.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "slidingInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Block mesh generator.\n"
        "  The ordering of vertex and face labels within a block as shown "
        "below.\n"
        "  For the local vertex numbering in the sequence 0 to 7:\n"
        "    Faces 0, 1 (x-direction) are left, right.\n"
        "    Faces 2, 3 (y-direction) are front, back.\n"
        "    Faces 4, 5 (z-direction) are bottom, top.\n"
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
        "blockTopology",
        "Write block edges and centres as obj files and exit"
    );
    argList::addBoolOption
    (
        "noClean",
        "Do not remove any existing polyMesh/ directory or files"
    );
    argList::addOption
    (
        "dict",
        "file",
        "Alternative dictionary for the blockMesh description"
    );
    argList::addBoolOption
    (
        "sets",
        "Write cellZones as cellSets too (for processing purposes)"
    );
    argList::addOption
    (
        "time",
        "time",
        "Specify a time to write mesh to (default: constant)"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    word regionName = polyMesh::defaultRegion;
    word regionPath;

    // Check if the region is specified otherwise mesh the default region
    if (args.readIfPresent("region", regionName))
    {
        Info<< nl << "Generating mesh for region " << regionName << endl;
        regionPath = regionName;
    }


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


    // Locate appropriate blockMeshDict
    #include "findBlockMeshDict.H"

    blockMesh blocks(meshDict, regionName);

    if (!blocks.valid())
    {
        // Could/should be Fatal?

        WarningIn(args.executable())
            << "Did not generate any blocks. Stopping." << nl << endl;

        return 1;
    }


    if (args.found("blockTopology"))
    {
        Info<< nl;

        // Write mesh as edges.
        {
            OFstream os(runTime.path()/"blockTopology.obj");

            Info<< "Writing block structure as obj format: "
                << os.name().name() << endl;

            blocks.writeTopology(os);
        }

        // Write centres of blocks
        {
            OFstream os(runTime.path()/"blockCentres.obj");

            Info<< "Writing block centres as obj format: "
                << os.name().name() << endl;

            const polyMesh& topo = blocks.topology();

            for (const point& cc :  topo.cellCentres())
            {
                os << "v " << cc.x() << ' ' << cc.y() << ' ' << cc.z() << nl;
            }
        }

        Info<< "\nEnd\n" << endl;

        return 0;
    }



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
            runTime.path()/meshInstance/regionPath/polyMesh::meshSubDir
        );

        if (exists(polyMeshPath))
        {
            if (exists(polyMeshPath/dictName))
            {
                Info<< "Not deleting polyMesh directory "
                    << runTime.relativePath(polyMeshPath) << nl
                    << "    because it contains " << dictName << endl;
            }
            else
            {
                Info<< "Deleting polyMesh directory "
                    << runTime.relativePath(polyMeshPath) << endl;
                rmDir(polyMeshPath);
            }
        }
    }


    Info<< nl << "Creating polyMesh from blockMesh" << endl;

    polyMesh mesh
    (
        IOobject
        (
            regionName,
            meshInstance,
            runTime
        ),
        pointField(blocks.points()),  // Copy, could we re-use space?
        blocks.cells(),
        blocks.patches(),
        blocks.patchNames(),
        blocks.patchDicts(),
        "defaultFaces",               // Default patch name
        emptyPolyPatch::typeName      // Default patch type
    );


    // Handle merging of patch pairs. Dictionary entry "mergePatchPairs"
    #include "mergePatchPairs.H"

    // Set any cellZones
    #include "addCellZones.H"


    // Detect any cyclic patches and force re-ordering of the faces
    {
        bool hasCyclic = false;
        for (const polyPatch& pp : mesh.boundaryMesh())
        {
            if (isA<cyclicPolyPatch>(pp))
            {
                hasCyclic = true;
                break;
            }
        }

        if (hasCyclic)
        {
            Info<< nl << "Detected cyclic patches; ordering boundary faces"
                << endl;
            const word oldInstance = mesh.instance();
            polyTopoChange meshMod(mesh);
            meshMod.changeMesh(mesh, false);
            mesh.setInstance(oldInstance);
        }
    }


    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing polyMesh with "
        << mesh.cellZones().size() << " cellZones";

    if (args.found("sets") && !mesh.cellZones().empty())
    {
        Info<< " (written as cellSets too)";
    }
    Info<< endl;

    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    if (args.found("sets"))
    {
        for (const cellZone& cz : mesh.cellZones())
        {
            cellSet(mesh, cz.name(), cz).write();
        }
    }

    // Write summary
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
