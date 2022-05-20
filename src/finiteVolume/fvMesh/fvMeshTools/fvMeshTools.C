/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "fvMeshTools.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "IOobjectList.H"
#include "basicFvGeometryScheme.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "polyBoundaryMeshEntries.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Adds patch if not yet there. Returns patchID.
Foam::label Foam::fvMeshTools::addPatch
(
    fvMesh& mesh,
    const polyPatch& patch,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const bool validBoundary
)
{
    auto& polyPatches = const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    label patchi = polyPatches.findPatchID(patch.name());
    if (patchi != -1)
    {
        // Already there
        return patchi;
    }


    // Append at end unless there are processor patches
    label insertPatchi = polyPatches.size();
    label startFacei = mesh.nFaces();

    if (!isA<processorPolyPatch>(patch))
    {
        forAll(polyPatches, patchi)
        {
            const polyPatch& pp = polyPatches[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                insertPatchi = patchi;
                startFacei = pp.start();
                break;
            }
        }
    }


    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    label sz = polyPatches.size();

    auto& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Add polyPatch at the end
    polyPatches.resize(sz+1);
    polyPatches.set
    (
        sz,
        patch.clone
        (
            polyPatches,
            insertPatchi,   //index
            0,              //size
            startFacei      //start
        )
    );
    fvPatches.resize(sz+1);
    fvPatches.set
    (
        sz,
        fvPatch::New
        (
            polyPatches[sz],  // point to newly added polyPatch
            mesh.boundary()
        )
    );


    do
    {
        #undef  doLocalCode
        #define doLocalCode(FieldType)                                \
        {                                                             \
            addPatchFields<FieldType>                                 \
            (                                                         \
                mesh, patchFieldDict, defaultPatchFieldType, Zero     \
            );                                                        \
        }

        // Volume fields
        doLocalCode(volScalarField);
        doLocalCode(volVectorField);
        doLocalCode(volSphericalTensorField);
        doLocalCode(volSymmTensorField);
        doLocalCode(volTensorField);

        // Surface fields
        doLocalCode(surfaceScalarField);
        doLocalCode(surfaceVectorField);
        doLocalCode(surfaceSphericalTensorField);
        doLocalCode(surfaceSymmTensorField);
        doLocalCode(surfaceTensorField);

        #undef doLocalCode
    }
    while (false);


    // Create reordering list
    // - patches before insert position stay as is
    // - patches after insert position move one up
    // - appended patch gets moved to insert position

    labelList oldToNew(identity(sz+1));
    for (label i = insertPatchi; i < sz; ++i)
    {
        oldToNew[i] = i+1;
    }
    oldToNew[sz] = insertPatchi;


    // Shuffle into place
    polyPatches.reorder(oldToNew, validBoundary);
    fvPatches.reorder(oldToNew);

    do
    {
        #undef  doLocalCode
        #define doLocalCode(FieldType)                                \
        {                                                             \
            reorderPatchFields<FieldType>(mesh, oldToNew);            \
        }

        // Volume fields
        doLocalCode(volScalarField);
        doLocalCode(volVectorField);
        doLocalCode(volSphericalTensorField);
        doLocalCode(volSymmTensorField);
        doLocalCode(volTensorField);

        // Surface fields
        doLocalCode(surfaceScalarField);
        doLocalCode(surfaceVectorField);
        doLocalCode(surfaceSphericalTensorField);
        doLocalCode(surfaceSymmTensorField);
        doLocalCode(surfaceTensorField);

        #undef doLocalCode
    }
    while (false);

    return insertPatchi;
}


void Foam::fvMeshTools::setPatchFields
(
    fvMesh& mesh,
    const label patchi,
    const dictionary& patchFieldDict
)
{
    do
    {
        #undef  doLocalCode
        #define doLocalCode(FieldType)                                \
        {                                                             \
            setPatchFields<FieldType>(mesh, patchi, patchFieldDict);  \
        }

        // Volume fields
        doLocalCode(volScalarField);
        doLocalCode(volVectorField);
        doLocalCode(volSphericalTensorField);
        doLocalCode(volSymmTensorField);
        doLocalCode(volTensorField);

        // Surface fields
        doLocalCode(surfaceScalarField);
        doLocalCode(surfaceVectorField);
        doLocalCode(surfaceSphericalTensorField);
        doLocalCode(surfaceSymmTensorField);
        doLocalCode(surfaceTensorField);

        #undef doLocalCode
    }
    while (false);
}


void Foam::fvMeshTools::zeroPatchFields(fvMesh& mesh, const label patchi)
{
    do
    {
        #undef  doLocalCode
        #define doLocalCode(FieldType)                                \
        {                                                             \
            setPatchFields<FieldType>(mesh, patchi, Zero);            \
        }

        // Volume fields
        doLocalCode(volScalarField);
        doLocalCode(volVectorField);
        doLocalCode(volSphericalTensorField);
        doLocalCode(volSymmTensorField);
        doLocalCode(volTensorField);

        // Surface fields
        doLocalCode(surfaceScalarField);
        doLocalCode(surfaceVectorField);
        doLocalCode(surfaceSphericalTensorField);
        doLocalCode(surfaceSymmTensorField);
        doLocalCode(surfaceTensorField);

        #undef doLocalCode
    }
    while (false);
}


// Deletes last patch
void Foam::fvMeshTools::trimPatches(fvMesh& mesh, const label nPatches)
{
    // Clear local fields and e.g. polyMesh globalMeshData.
    mesh.clearOut();

    auto& polyPatches = const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    auto& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    if (polyPatches.empty())
    {
        FatalErrorInFunction
            << "No patches in mesh"
            << abort(FatalError);
    }

    label nFaces = 0;
    for (label patchi = nPatches; patchi < polyPatches.size(); patchi++)
    {
        nFaces += polyPatches[patchi].size();
    }
    reduce(nFaces, sumOp<label>());

    if (nFaces)
    {
        FatalErrorInFunction
            << "There are still " << nFaces
            << " faces in " << polyPatches.size()-nPatches
            << " patches to be deleted" << abort(FatalError);
    }

    // Remove actual patches
    polyPatches.resize(nPatches);
    fvPatches.resize(nPatches);

    do
    {
        #undef  doLocalCode
        #define doLocalCode(FieldType)                                \
        {                                                             \
            trimPatchFields<FieldType>(mesh, nPatches);               \
        }

        // Volume fields
        doLocalCode(volScalarField);
        doLocalCode(volVectorField);
        doLocalCode(volSphericalTensorField);
        doLocalCode(volSymmTensorField);
        doLocalCode(volTensorField);

        // Surface fields
        doLocalCode(surfaceScalarField);
        doLocalCode(surfaceVectorField);
        doLocalCode(surfaceSphericalTensorField);
        doLocalCode(surfaceSymmTensorField);
        doLocalCode(surfaceTensorField);

        #undef doLocalCode
    }
    while (false);
}


void Foam::fvMeshTools::reorderPatches
(
    fvMesh& mesh,
    const labelList& oldToNew,
    const label nNewPatches,
    const bool validBoundary
)
{
    auto& polyPatches = const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    auto& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Shuffle into place
    polyPatches.reorder(oldToNew, validBoundary);
    fvPatches.reorder(oldToNew);

    do
    {
        #undef  doLocalCode
        #define doLocalCode(FieldType)                                \
        {                                                             \
            reorderPatchFields<FieldType>(mesh, oldToNew);            \
        }

        // Volume fields
        doLocalCode(volScalarField);
        doLocalCode(volVectorField);
        doLocalCode(volSphericalTensorField);
        doLocalCode(volSymmTensorField);
        doLocalCode(volTensorField);

        // Surface fields
        doLocalCode(surfaceScalarField);
        doLocalCode(surfaceVectorField);
        doLocalCode(surfaceSphericalTensorField);
        doLocalCode(surfaceSymmTensorField);
        doLocalCode(surfaceTensorField);

        #undef doLocalCode
    }
    while (false);

    // Remove last.
    trimPatches(mesh, nNewPatches);
}


Foam::labelList Foam::fvMeshTools::removeEmptyPatches
(
    fvMesh& mesh,
    const bool validBoundary
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList newToOld(pbm.size());
    labelList oldToNew(pbm.size(), -1);
    label newI = 0;


    // Assumes all non-coupled boundaries are on all processors!
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (!isA<processorPolyPatch>(pp))
        {
            label nFaces = pp.size();
            if (validBoundary)
            {
                reduce(nFaces, sumOp<label>());
            }

            if (nFaces > 0)
            {
                newToOld[newI] = patchI;
                oldToNew[patchI] = newI++;
            }
        }
    }

    // Same for processor patches (but need no reduction)
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (isA<processorPolyPatch>(pp) && pp.size())
        {
            newToOld[newI] = patchI;
            oldToNew[patchI] = newI++;
        }
    }

    newToOld.resize(newI);

    // Move all deletable patches to the end
    forAll(oldToNew, patchI)
    {
        if (oldToNew[patchI] == -1)
        {
            oldToNew[patchI] = newI++;
        }
    }

    reorderPatches(mesh, oldToNew, newToOld.size(), validBoundary);

    return newToOld;
}


void Foam::fvMeshTools::setBasicGeometry(fvMesh& mesh)
{
    // Set the fvGeometryScheme to basic since it does not require
    // any parallel communication to construct the geometry. During
    // redistributePar one might temporarily end up with processors
    // with zero procBoundaries. Normally procBoundaries trigger geometry
    // calculation (e.g. send over cellCentres) so on the processors without
    // procBoundaries this will not happen. The call to the geometry calculation
    // is not synchronised and this might lead to a hang for geometry schemes
    // that do require synchronisation

    tmp<fvGeometryScheme> basicGeometry
    (
        new basicFvGeometryScheme(mesh, dictionary())
    );
    mesh.geometry(basicGeometry);
}


Foam::autoPtr<Foam::fvMesh>
Foam::fvMeshTools::newMesh
(
    const IOobject& io,
    const bool masterOnlyReading,
    const bool verbose
)
{
    // Region name
    // ~~~~~~~~~~~

    const fileName meshSubDir
    (
        polyMesh::regionName(io.name()) / polyMesh::meshSubDir
    );


    fileName facesInstance;
    fileName pointsInstance;

    // Patch types
    // ~~~~~~~~~~~
    // Read and scatter master patches (without reading master mesh!)

    PtrList<entry> patchEntries;
    if (Pstream::master())
    {
        const bool oldParRun = Pstream::parRun(false);

        facesInstance = io.time().findInstance
        (
            meshSubDir,
            "faces",
            IOobject::MUST_READ
        );
        pointsInstance = io.time().findInstance
        (
            meshSubDir,
            "points",
            IOobject::MUST_READ
        );

        patchEntries = polyBoundaryMeshEntries
        (
            IOobject
            (
                "boundary",
                facesInstance,
                meshSubDir,
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        Pstream::parRun(oldParRun);
    }

    // Broadcast information to all
    Pstream::broadcasts
    (
        UPstream::worldComm,
        patchEntries,
        facesInstance,
        pointsInstance
    );


    // Dummy meshes
    // ~~~~~~~~~~~~

    // Set up to read-if-present.
    // Note: does not search for mesh so set instance explicitly

    IOobject meshIO(io);
    meshIO.instance() = facesInstance;
    meshIO.readOpt(IOobject::READ_IF_PRESENT);

    // For mesh components (points, faces, ...)
    IOobject cmptIO(meshIO, "points", meshSubDir);
    cmptIO.readOpt(IOobject::MUST_READ);
    cmptIO.writeOpt(IOobject::NO_WRITE);
    cmptIO.registerObject(false);


    // Check who has a mesh
    const fileName meshDir = io.time().path()/facesInstance/meshSubDir;
    bool haveMesh = isDir(meshDir);
    if (masterOnlyReading && !Pstream::master())
    {
        haveMesh = false;
        meshIO.readOpt(IOobject::NO_READ);
    }

    if (!haveMesh)
    {
        cmptIO.readOpt(IOobject::NO_READ);
    }


    // Read mesh
    // ~~~~~~~~~
    // Now all processors use supplied points,faces etc
    // Note: solution, schemes are also using the supplied IOobject so
    //       on slave will be NO_READ, on master READ_IF_PRESENT. This will
    //       conflict with e.g. timeStampMaster reading so switch off.
    // Note: v2006 used the READ_IF_PRESENT flag in the meshIO.readOpt(). v2012
    //       (correctly) does no longer so below code explicitly addFvPatches
    //       using the separately read boundary file.

    const auto oldCheckType = IOobject::fileModificationChecking;
    IOobject::fileModificationChecking = IOobject::timeStamp;


    // Points
    cmptIO.instance() = pointsInstance;
    cmptIO.rename("points");
    pointIOField points(cmptIO);

    // All other mesh components use facesInstance ...
    cmptIO.instance() = facesInstance;

    // Faces
    cmptIO.rename("faces");
    faceCompactIOList faces(cmptIO);

    // Face owner
    cmptIO.rename("owner");
    labelIOList owner(cmptIO);

    // Face neighbour
    cmptIO.rename("neighbour");
    labelIOList neighbour(cmptIO);


    auto meshPtr = autoPtr<fvMesh>::New
    (
        meshIO,
        std::move(points),
        std::move(faces),
        std::move(owner),
        std::move(neighbour)
    );
    fvMesh& mesh = *meshPtr;

    IOobject::fileModificationChecking = oldCheckType;


    // Some processors without patches? - add patches

    if (returnReduce(mesh.boundary().empty(), orOp<bool>()))
    {
        DynamicList<polyPatch*> newPatches(patchEntries.size());

        if (mesh.boundary().size() == patchEntries.size())
        {
            // Assumably we have correctly read the boundary and can clone.
            // Note: for
            // v2012 onwards this is probably never the case and this whole
            // section can be removed.
            forAll(mesh.boundary(), patchi)
            {
                newPatches.append
                (
                    mesh.boundaryMesh()[patchi].clone(mesh.boundaryMesh()).ptr()
                );
            }
        }
        else
        {
            // Use patchEntries, which were read on master and broadcast
            label nPatches = 0;

            const bool isEmptyMesh = (mesh.nInternalFaces() == 0);

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
                    // Unlikely to work with inter-mixed proc-patches anyhow
                    // - could break out of loop here
                }
                else
                {
                    dictionary patchDict(e.dict());

                    if (isEmptyMesh)
                    {
                        patchDict.set("nFaces", 0);
                        patchDict.set("startFace", 0);
                    }

                    newPatches.append
                    (
                        polyPatch::New
                        (
                            name,
                            patchDict,
                            nPatches++,
                            mesh.boundaryMesh()
                        ).ptr()
                    );
                }
            }
        }

        mesh.removeFvBoundary();
        mesh.addFvPatches(newPatches);
    }


    // Determine zones
    // ~~~~~~~~~~~~~~~

    wordList pointZoneNames(mesh.pointZones().names());
    wordList faceZoneNames(mesh.faceZones().names());
    wordList cellZoneNames(mesh.cellZones().names());

    Pstream::broadcasts
    (
        UPstream::worldComm,
        pointZoneNames,
        faceZoneNames,
        cellZoneNames
    );

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

        if (pz.size() || fz.size() || cz.size())
        {
            mesh.addZones(pz, fz, cz);
        }
    }

    return meshPtr;
}


Foam::autoPtr<Foam::fvMesh>
Foam::fvMeshTools::loadOrCreateMesh
(
    const IOobject& io,
    const bool decompose,
    const bool verbose
)
{

    // Region name
    // ~~~~~~~~~~~

    const fileName meshSubDir
    (
        polyMesh::regionName(io.name()) / polyMesh::meshSubDir
    );


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
                io.instance(),
                meshSubDir,
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        Pstream::parRun(oldParRun);
    }

    // Broadcast: send patches to all
    Pstream::broadcast(patchEntries);  // == worldComm;


    // Dummy meshes
    // ~~~~~~~~~~~~

    // Check who has or needs a mesh.
    // For 'decompose', only need mesh on master.
    // Otherwise check for presence of the "faces" file

    bool haveMesh =
    (
        decompose
      ? Pstream::master()
      : fileHandler().isFile
        (
            fileHandler().filePath
            (
                io.time().path()/io.instance()/meshSubDir/"faces"
            )
        )
    );


    if (!haveMesh)
    {
        const bool oldParRun = Pstream::parRun(false);

        // Create dummy mesh - on procs that don't already have a mesh
        fvMesh dummyMesh
        (
            IOobject(io, IOobject::NO_READ, IOobject::AUTO_WRITE),
            Foam::zero{},
            false
        );

        // Add patches
        polyPatchList patches(patchEntries.size());
        label nPatches = 0;

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
                // Stop at the first processor patch.
                // - logic fails with inter-mixed proc-patches anyhow
                break;
            }
            else
            {
                dictionary patchDict(e.dict());
                patchDict.set("nFaces", 0);
                patchDict.set("startFace", 0);

                patches.set
                (
                    patchi,
                    polyPatch::New
                    (
                        name,
                        patchDict,
                        nPatches++,
                        dummyMesh.boundaryMesh()
                    )
                );
            }
        }
        patches.resize(nPatches);
        dummyMesh.addFvPatches(patches, false);  // No parallel comms


        // Add some dummy zones so upon reading it does not read them
        // from the undecomposed case. Should be done as extra argument to
        // regIOobject::readStream?

        List<pointZone*> pz
        (
            1,
            new pointZone("dummyPointZone", 0, dummyMesh.pointZones())
        );
        List<faceZone*> fz
        (
            1,
            new faceZone("dummyFaceZone", 0, dummyMesh.faceZones())
        );
        List<cellZone*> cz
        (
            1,
            new cellZone("dummyCellZone", 0, dummyMesh.cellZones())
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

    //Pout<< "Reading mesh from " << io.objectRelPath() << endl;
    auto meshPtr = autoPtr<fvMesh>::New(io);
    fvMesh& mesh = *meshPtr;

    // Make sure to use a non-parallel geometry calculation method
    fvMeshTools::setBasicGeometry(mesh);


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
    wordList faceZoneNames(mesh.faceZones().names());
    wordList cellZoneNames(mesh.cellZones().names());
    Pstream::broadcasts
    (
        UPstream::worldComm,
        pointZoneNames,
        faceZoneNames,
        cellZoneNames
    );

    if (!haveMesh)
    {
        // Add the zones. Make sure to remove the old dummy ones first
        mesh.pointZones().clear();
        mesh.faceZones().clear();
        mesh.cellZones().clear();

        List<pointZone*> pz(pointZoneNames.size());
        forAll(pointZoneNames, i)
        {
            pz[i] = new pointZone(pointZoneNames[i], i, mesh.pointZones());
        }
        List<faceZone*> fz(faceZoneNames.size());
        forAll(faceZoneNames, i)
        {
            fz[i] = new faceZone(faceZoneNames[i], i, mesh.faceZones());
        }
        List<cellZone*> cz(cellZoneNames.size());
        forAll(cellZoneNames, i)
        {
            cz[i] = new cellZone(cellZoneNames[i], i, mesh.cellZones());
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

        pointSetNames = objects.sortedNames<pointSet>();
        faceSetNames = objects.sortedNames<faceSet>();
        cellSetNames = objects.sortedNames<cellSet>();
    }
    Pstream::broadcasts
    (
        UPstream::worldComm,
        pointSetNames,
        faceSetNames,
        cellSetNames
    );

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
    // and processor patches are correct
    mesh.boundaryMesh().checkDefinition(verbose);
    mesh.boundaryMesh().checkParallelSync(verbose);

    // Check names of zones are equal
    mesh.cellZones().checkDefinition(verbose);
    mesh.cellZones().checkParallelSync(verbose);
    mesh.faceZones().checkDefinition(verbose);
    mesh.faceZones().checkParallelSync(verbose);
    mesh.pointZones().checkDefinition(verbose);
    mesh.pointZones().checkParallelSync(verbose);

    return meshPtr;
}


void Foam::fvMeshTools::createDummyFvMeshFiles
(
    const objectRegistry& mesh,
    const word& regionName,
    const bool verbose
)
{
    // Create dummy system/fv*
    {
        IOobject io
        (
            "fvSchemes",
            mesh.time().system(),
            regionName,
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>(false))
        {
            if (verbose)
            {
                Info<< "Writing dummy " << regionName/io.name() << endl;
            }
            dictionary dict;
            dict.add("divSchemes", dictionary());
            dict.add("gradSchemes", dictionary());
            dict.add("laplacianSchemes", dictionary());

            IOdictionary(io, dict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>(false))
        {
            if (verbose)
            {
                Info<< "Writing dummy " << regionName/io.name() << endl;
            }
            dictionary dict;
            IOdictionary(io, dict).regIOobject::write();
        }
    }
}


// ************************************************************************* //
