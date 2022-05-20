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

#include "faMeshTools.H"
#include "faBoundaryMeshEntries.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "polyMesh.H"
#include "processorFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::faMeshTools::unregisterMesh(const faMesh& mesh)
{
    auto& obr = const_cast<objectRegistry&>(mesh.thisDb());

    // Checkout by name (casting ambiguity)
    obr.checkOut(faMesh::typeName);
    obr.checkOut("faBoundaryMesh");
    obr.checkOut("faSchemes");
    obr.checkOut("faSolution");
}


void Foam::faMeshTools::forceDemandDriven(faMesh& mesh)
{
    (void)mesh.globalData();

    (void)mesh.Le();
    (void)mesh.magLe();
    (void)mesh.areaCentres();
    (void)mesh.edgeCentres();

    (void)mesh.faceAreaNormals();
    (void)mesh.edgeAreaNormals();
    (void)mesh.pointAreaNormals();
    (void)mesh.faceCurvatures();
    (void)mesh.edgeTransformTensors();
}


Foam::autoPtr<Foam::faMesh> Foam::faMeshTools::newMesh
(
    const IOobject& io,
    const polyMesh& pMesh,
    const bool masterOnlyReading,
    const bool verbose
)
{
    // Region name
    // ~~~~~~~~~~~

    const fileName meshSubDir
    (
        pMesh.regionName() / faMesh::meshSubDir
    );


    fileName facesInstance;

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
            "faceLabels",
            IOobject::MUST_READ
        );

        patchEntries = faBoundaryMeshEntries
        (
            IOobject
            (
                "faBoundary",
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
        facesInstance
    );


    // Dummy meshes
    // ~~~~~~~~~~~~

    // Set up to read-if-present. Note: does not search for mesh so set
    // instance explicitly

    IOobject meshIO(io);
    meshIO.instance() = facesInstance;
    meshIO.readOpt(IOobject::READ_IF_PRESENT);

    // For mesh components (faceLabels, ...)
    IOobject cmptIO(meshIO, "faceLabels", meshSubDir);
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

    const auto oldCheckType = IOobject::fileModificationChecking;
    IOobject::fileModificationChecking = IOobject::timeStamp;


    // faceLabels
    cmptIO.rename("faceLabels");
    labelIOList faceLabels(cmptIO);


    auto meshPtr = autoPtr<faMesh>::New
    (
        pMesh,
        std::move(faceLabels),
        meshIO
    );
    auto& mesh = *meshPtr;

    IOobject::fileModificationChecking = oldCheckType;


    // Some processors without patches? - add patches

    if (returnReduce(mesh.boundary().empty(), orOp<bool>()))
    {
        // Use patchEntries, which were read on master and broadcast

        faPatchList patches(patchEntries.size());
        label nPatches = 0;

        const bool isEmptyMesh = (mesh.faceLabels().empty());

        forAll(patchEntries, patchi)
        {
            const entry& e = patchEntries[patchi];
            const word type(e.dict().get<word>("type"));
            const word& name = e.keyword();

            if
            (
                type == processorFaPatch::typeName
            )
            {
                // Stop at the first processor patch.
                // - logic will not work with inter-mixed proc-patches anyhow
                break;
            }
            else
            {
                dictionary patchDict(e.dict());

                if (isEmptyMesh)
                {
                    patchDict.set("edgeLabels", labelList());
                }

                patches.set
                (
                    patchi,
                    faPatch::New
                    (
                        name,
                        patchDict,
                        nPatches++,
                        mesh.boundary()
                    )
                );
            }
        }

        patches.resize(nPatches);
        mesh.addFaPatches(patches, false);  // No parallel comms
    }

    // Recreate basic geometry, globalMeshData etc.
    mesh.init(false);
    (void)mesh.globalData();

    return meshPtr;
}


Foam::autoPtr<Foam::faMesh> Foam::faMeshTools::loadOrCreateMesh
(
    const IOobject& io,
    const polyMesh& pMesh,
    const bool decompose,
    const bool verbose
)
{
    // Region name
    // ~~~~~~~~~~~

    const fileName meshSubDir
    (
        pMesh.regionName() / faMesh::meshSubDir
    );


    // Patch types
    // ~~~~~~~~~~~
    // Read and scatter master patches (without reading master mesh!)

    PtrList<entry> patchEntries;
    if (Pstream::master())
    {
        const bool oldParRun = Pstream::parRun(false);

        patchEntries = faBoundaryMeshEntries
        (
            IOobject
            (
                "faBoundary",
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

    /// Info<< patchEntries << nl;

    // Dummy meshes
    // ~~~~~~~~~~~~

    // Check who has or needs a mesh.
    // For 'decompose', only need mesh on master.
    // Otherwise check for presence of the "faceLabels" file

    bool haveMesh =
    (
        decompose
      ? Pstream::master()
      : fileHandler().isFile
        (
            fileHandler().filePath
            (
                io.time().path()/io.instance()/meshSubDir/"faceLabels"
            )
        )
    );


    if (!haveMesh)
    {
        const bool oldParRun = Pstream::parRun(false);

        // Create dummy mesh - on procs that don't already have a mesh
        faMesh dummyMesh
        (
            pMesh,
            labelList(),
            IOobject(io, IOobject::NO_READ, IOobject::AUTO_WRITE)
        );

        // Add patches
        faPatchList patches(patchEntries.size());
        label nPatches = 0;

        forAll(patchEntries, patchi)
        {
            const entry& e = patchEntries[patchi];
            const word type(e.dict().get<word>("type"));
            const word& name = e.keyword();

            if
            (
                type == processorFaPatch::typeName
            )
            {
                // Stop at the first processor patch.
                // - logic will not work with inter-mixed proc-patches anyhow
                break;
            }
            else
            {
                dictionary patchDict(e.dict());
                patchDict.set("edgeLabels", labelList());

                patches.set
                (
                    patchi,
                    faPatch::New
                    (
                        name,
                        patchDict,
                        nPatches++,
                        dummyMesh.boundary()
                    )
                );
            }
        }

        patches.resize(nPatches);
        dummyMesh.addFaPatches(patches, false);  // No parallel comms

        // Bad hack, but the underlying polyMesh is NO_WRITE
        // so it does not create the faMesh subDir for us...
        Foam::mkDir(dummyMesh.boundary().path());

        //Pout<< "Writing dummy mesh to " << dummyMesh.boundary().path()
        //    << endl;
        dummyMesh.write();

        Pstream::parRun(oldParRun);  // Restore parallel state
    }

    // Read mesh
    // ~~~~~~~~~
    // Now all processors have a (possibly zero size) mesh so read in
    // parallel

    /// Pout<< "Reading area mesh from " << io.objectRelPath() << endl;
    auto meshPtr = autoPtr<faMesh>::New(pMesh, false);
    faMesh& mesh = *meshPtr;

    // Sync patches
    // ~~~~~~~~~~~~

    if (!Pstream::master() && haveMesh)
    {
        // Check master names against mine

        const faBoundaryMesh& patches = mesh.boundary();

        forAll(patchEntries, patchi)
        {
            const entry& e = patchEntries[patchi];
            const word type(e.dict().get<word>("type"));
            const word& name = e.keyword();

            if
            (
                type == processorFaPatch::typeName
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


    // Recreate basic geometry, globalMeshData etc.
    mesh.init(false);
    (void)mesh.globalData();

    /// #if 0
    /// faMeshTools::forceDemandDriven(mesh);
    /// faMeshTools::unregisterMesh(mesh);
    /// #endif

    // Do some checks.

    // Check if the boundary definition is unique
    // and processor patches are correct
    mesh.boundary().checkDefinition(verbose);
    mesh.boundary().checkParallelSync(verbose);

    return meshPtr;
}


// ************************************************************************* //
