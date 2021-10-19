/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "loadOrCreateMesh.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "Time.H"
#include "polyBoundaryMeshEntries.H"
#include "IOobjectList.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "basicFvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Read mesh if available. Otherwise create empty mesh with same non-proc
// patches as proc0 mesh. Requires:
//  - all processors to have all patches (and in same order).
//  - io.instance() set to facesInstance
Foam::autoPtr<Foam::fvMesh> Foam::loadOrCreateMesh
(
    const bool decompose,
    const IOobject& io
)
{
    // Region name
    // ~~~~~~~~~~~

    fileName meshSubDir;

    if (io.name() == polyMesh::defaultRegion)
    {
        meshSubDir = polyMesh::meshSubDir;
    }
    else
    {
        meshSubDir = io.name()/polyMesh::meshSubDir;
    }


    // Patch types
    // ~~~~~~~~~~~
    // Read and scatter master patches (without reading master mesh!)

    PtrList<entry> patchEntries;
    if (Pstream::master())
    {
        const bool oldParRun = Pstream::parRun(false);

        patchEntries = polyBoundaryMeshEntries
        (
            IOobject
            (
                "boundary",
                io.instance(),  //facesInstance,
                meshSubDir,
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        Pstream::parRun(oldParRun);

        // Send patches
        for (const int slave : Pstream::subProcs())
        {
            OPstream toSlave(Pstream::commsTypes::scheduled, slave);
            toSlave << patchEntries;
        }
    }
    else
    {
        // Receive patches
        IPstream fromMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );
        fromMaster >> patchEntries;
    }



    // Dummy meshes
    // ~~~~~~~~~~~~

    // Check who has a mesh
    bool haveMesh;
    if (decompose)
    {
        // Mesh needs to be present on the master only
        haveMesh = Pstream::master();
    }
    else
    {
        const fileName facesFile
        (
            io.time().path()
           /io.instance()   //facesInstance
           /meshSubDir
           /"faces"
        );

        // Check presence of the searched-for faces file
        haveMesh = fileHandler().isFile(fileHandler().filePath(facesFile));
    }

    if (!haveMesh)
    {
        const bool oldParRun = Pstream::parRun(false);


        // Create dummy mesh. Only used on procs that don't have mesh.
        IOobject noReadIO(io);
        noReadIO.readOpt(IOobject::NO_READ);
        noReadIO.writeOpt(IOobject::AUTO_WRITE);
        fvMesh dummyMesh(noReadIO, Zero, false);

        // Add patches
        List<polyPatch*> patches(patchEntries.size());
        label nPatches = 0;

        forAll(patchEntries, patchi)
        {
            const entry& e = patchEntries[patchi];
            const word type(e.dict().get<word>("type"));
            const word& name = e.keyword();

            if
            (
                type != processorPolyPatch::typeName
             && type != processorCyclicPolyPatch::typeName
            )
            {
                dictionary patchDict(e.dict());
                patchDict.set("nFaces", 0);
                patchDict.set("startFace", 0);

                patches[patchi] = polyPatch::New
                (
                    name,
                    patchDict,
                    nPatches++,
                    dummyMesh.boundaryMesh()
                ).ptr();
            }
        }
        patches.setSize(nPatches);
        dummyMesh.addFvPatches(patches, false);  // no parallel comms


        // Add some dummy zones so upon reading it does not read them
        // from the undecomposed case. Should be done as extra argument to
        // regIOobject::readStream?
        List<pointZone*> pz
        (
            1,
            new pointZone
            (
                "dummyPointZone",
                0,
                dummyMesh.pointZones()
            )
        );
        List<faceZone*> fz
        (
            1,
            new faceZone
            (
                "dummyFaceZone",
                0,
                dummyMesh.faceZones()
            )
        );
        List<cellZone*> cz
        (
            1,
            new cellZone
            (
                "dummyCellZone",
                0,
                dummyMesh.cellZones()
            )
        );
        dummyMesh.addZones(pz, fz, cz);
        dummyMesh.pointZones().clear();
        dummyMesh.faceZones().clear();
        dummyMesh.cellZones().clear();
        //Pout<< "Writing dummy mesh to " << dummyMesh.polyMesh::objectPath()
        //    << endl;
        dummyMesh.write();

        Pstream::parRun(oldParRun);  // Restore parallel state
    }



    // Read mesh
    // ~~~~~~~~~
    // Now all processors have a (possibly zero size) mesh so read in
    // parallel

    //Pout<< "Reading mesh from " << io.objectPath() << endl;
    auto meshPtr = autoPtr<fvMesh>::New(io);
    fvMesh& mesh = *meshPtr;

    // Make sure to use a non-parallel geometry calculation method
    {
        tmp<fvGeometryScheme> basicGeometry
        (
            fvGeometryScheme::New
            (
                mesh,
                dictionary(),
                basicFvGeometryScheme::typeName
            )
        );
        mesh.geometry(basicGeometry);
    }


    // Sync patches
    // ~~~~~~~~~~~~

    if (!Pstream::master() && haveMesh)
    {
        // Check master names against mine

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(patchEntries, patchi)
        {
            const entry& e = patchEntries[patchi];
            const word type(e.dict().get<word>("type"));
            const word& name = e.keyword();

            if
            (
                type == processorPolyPatch::typeName
             || type == processorCyclicPolyPatch::typeName
            )
            {
                break;
            }

            if (patchi >= patches.size())
            {
                FatalErrorInFunction
                    << "Non-processor patches not synchronised."
                    << endl
                    << "Processor " << Pstream::myProcNo()
                    << " has only " << patches.size()
                    << " patches, master has "
                    << patchi
                    << exit(FatalError);
            }

            if
            (
                type != patches[patchi].type()
             || name != patches[patchi].name()
            )
            {
                FatalErrorInFunction
                    << "Non-processor patches not synchronised."
                    << endl
                    << "Master patch " << patchi
                    << " name:" << type
                    << " type:" << type << endl
                    << "Processor " << Pstream::myProcNo()
                    << " patch " << patchi
                    << " has name:" << patches[patchi].name()
                    << " type:" << patches[patchi].type()
                    << exit(FatalError);
            }
        }
    }


    // Determine zones
    // ~~~~~~~~~~~~~~~

    wordList pointZoneNames(mesh.pointZones().names());
    Pstream::scatter(pointZoneNames);
    wordList faceZoneNames(mesh.faceZones().names());
    Pstream::scatter(faceZoneNames);
    wordList cellZoneNames(mesh.cellZones().names());
    Pstream::scatter(cellZoneNames);

    if (!haveMesh)
    {
        // Add the zones. Make sure to remove the old dummy ones first
        mesh.pointZones().clear();
        mesh.faceZones().clear();
        mesh.cellZones().clear();

        List<pointZone*> pz(pointZoneNames.size());
        forAll(pointZoneNames, i)
        {
            pz[i] = new pointZone
            (
                pointZoneNames[i],
                i,
                mesh.pointZones()
            );
        }
        List<faceZone*> fz(faceZoneNames.size());
        forAll(faceZoneNames, i)
        {
            fz[i] = new faceZone
            (
                faceZoneNames[i],
                i,
                mesh.faceZones()
            );
        }
        List<cellZone*> cz(cellZoneNames.size());
        forAll(cellZoneNames, i)
        {
            cz[i] = new cellZone
            (
                cellZoneNames[i],
                i,
                mesh.cellZones()
            );
        }
        mesh.addZones(pz, fz, cz);
    }


    // Determine sets
    // ~~~~~~~~~~~~~~

    wordList pointSetNames;
    wordList faceSetNames;
    wordList cellSetNames;
    if (Pstream::master())
    {
        // Read sets
        const bool oldParRun = Pstream::parRun(false);
        IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
        Pstream::parRun(oldParRun);

        pointSetNames = objects.sortedNames(pointSet::typeName);
        faceSetNames = objects.sortedNames(faceSet::typeName);
        cellSetNames = objects.sortedNames(cellSet::typeName);
    }
    Pstream::scatter(pointSetNames);
    Pstream::scatter(faceSetNames);
    Pstream::scatter(cellSetNames);

    if (!haveMesh)
    {
        for (const word& setName : pointSetNames)
        {
            pointSet(mesh, setName, 0).write();
        }
        for (const word& setName : faceSetNames)
        {
            faceSet(mesh, setName, 0).write();
        }
        for (const word& setName : cellSetNames)
        {
            cellSet(mesh, setName, 0).write();
        }
    }


    // Force recreation of globalMeshData.
    mesh.globalData();


    // Do some checks.

    // Check if the boundary definition is unique
    mesh.boundaryMesh().checkDefinition(true);
    // Check if the boundary processor patches are correct
    mesh.boundaryMesh().checkParallelSync(true);
    // Check names of zones are equal
    mesh.cellZones().checkDefinition(true);
    mesh.cellZones().checkParallelSync(true);
    mesh.faceZones().checkDefinition(true);
    mesh.faceZones().checkParallelSync(true);
    mesh.pointZones().checkDefinition(true);
    mesh.pointZones().checkParallelSync(true);

    return meshPtr;
}


bool Foam::removeEmptyDir(const fileName& path)
{
    // Return true if empty directory. Note bypass of fileHandler to be
    // consistent with polyMesh.removeFiles for now.

    {
        fileNameList files
        (
            Foam::readDir
            (
                path,
                fileName::FILE,
                false,                  // filterGz
                false                   // followLink
            )
        );
        if (files.size())
        {
            return false;
        }
    }
    {
        fileNameList dirs
        (
            Foam::readDir
            (
                path,
                fileName::DIRECTORY,
                false,                  // filterGz
                false                   // followLink
            )
        );
        if (dirs.size())
        {
            return false;
        }
    }
    {
        fileNameList links
        (
            Foam::readDir
            (
                path,
                fileName::LINK,
                false,                  // filterGz
                false                   // followLink
            )
        );
        if (links.size())
        {
            return false;
        }
    }
    {
        fileNameList other
        (
            Foam::readDir
            (
                path,
                fileName::UNDEFINED,
                false,                  // filterGz
                false                   // followLink
            )
        );
        if (other.size())
        {
            return false;
        }
    }

    // Avoid checking success of deletion since initial path might not
    // exist (e.g. contain 'region0'). Will stop when trying to delete
    // parent directory anyway since now not empty.
    Foam::rm(path);
    return true;
}


void Foam::removeEmptyDirs(const fileName& meshPath)
{
    // Delete resulting directory if empty
    fileName path(meshPath);
    path.clean();

    // Do subdirectories
    {
        const fileNameList dirs
        (
            Foam::readDir
            (
                path,
                fileName::DIRECTORY,
                false,                  // filterGz
                false                   // followLink
            )
        );
        for (const auto& dir : dirs)
        {
            removeEmptyDirs(path/dir);
        }
    }

    removeEmptyDir(path);
}


// ************************************************************************* //
