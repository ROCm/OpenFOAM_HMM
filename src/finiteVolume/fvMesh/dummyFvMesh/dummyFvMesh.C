/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "dummyFvMesh.H"
#include "Time.H"
#include "polyBoundaryMeshEntries.H"
#include "IOobjectList.H"
#include "fieldDictionary.H"
#include "topoSet.H"

#include "calculatedFvPatchField.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicACMIPolyPatch.H"

#include "fvPatchField.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dummyFvMesh, 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dummyFvMesh::fvPatchFieldExists(const word& patchType)
{
    if
    (
        fvPatchField<scalar>::dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<vector>::dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<sphericalTensor>::
            dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<symmTensor>::
            dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<tensor>::dictionaryConstructorTablePtr_->found(patchType)
    )
    {
        return true;
    }

    return false;
}

Foam::autoPtr<Foam::fvMesh> Foam::dummyFvMesh::singleCellMesh
(
    const Time& runTime,
    const scalar d
)
{
    pointField points(8);
    points[0] = vector(0, 0, 0);
    points[1] = vector(d, 0, 0);
    points[2] = vector(d, d, 0);
    points[3] = vector(0, d, 0);
    points[4] = vector(0, 0, d);
    points[5] = vector(d, 0, d);
    points[6] = vector(d, d, d);
    points[7] = vector(0, d, d);

    faceList faces = cellModel::ref(cellModel::HEX).modelFaces();

    autoPtr<fvMesh> meshPtr
    (
        new Foam::fvMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            std::move(points),
            std::move(faces),
            labelList(6, Zero),
            labelList()
        )
    );

    fvMesh& mesh = meshPtr();

    List<polyPatch*> patches(1);

    patches[0] = new emptyPolyPatch
    (
        "boundary",
        6,
        0,
        0,
        mesh.boundaryMesh(),
        emptyPolyPatch::typeName
    );

    mesh.addFvPatches(patches);

    return meshPtr;
}


bool Foam::dummyFvMesh::setPatchEntries
(
    const Time& runTime,
    const word& instance,
    dictionary& patchEntries,
    label& nPatchWithFace
)
{
    IOobject boundaryIO
    (
        "boundary",
        instance,
        polyMesh::meshSubDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    const wordHashSet constraintPatches(polyPatch::constraintTypes());

    if (boundaryIO.typeHeaderOk<polyBoundaryMesh>(true))
    {
        polyBoundaryMeshEntries allPatchEntries(boundaryIO);

        Info<< "Creating dummy mesh using " << allPatchEntries.path() << endl;

        for (const entry& e : allPatchEntries)
        {
            const word type(e.dict().lookup("type"));

            if (!constraintPatches.found(type))
            {
                if (readLabel(e.dict().lookup("nFaces")))
                {
                    ++nPatchWithFace;
                }
                patchEntries.add(e.keyword(), e.dict());
            }
        }

        return true;
    }
    else
    {
        // No boundary file - try reading from a field
        IOobjectList objects(runTime, runTime.timeName());

        if (objects.empty())
        {
            // No fields - cannot create a mesh
            return false;
        }

        const IOobject& io = *objects.begin()();

        if (io.instance() == runTime.constant())
        {
            FatalErrorInFunction
                << "No time directories found for field reading"
                << exit(FatalError);
        }

        const fieldDictionary fieldDict(io, io.headerClassName());

        Info<< "Creating dummy mesh from field "
            << fieldDict.objectPath()
            << endl;

        const dictionary& boundaryFieldDict =
            fieldDict.subDict("boundaryField");

        for (const entry& e : boundaryFieldDict)
        {
            const word type(e.dict().lookup("type"));

            if (fvPatchFieldExists(type))
            {
                if (!constraintPatches.found(type))
                {
                    ++nPatchWithFace;
                    dictionary dummyEntries;
                    dummyEntries.add("startFace", 0);
                    dummyEntries.add("nFaces", 1);
                    dummyEntries.add("type", "wall"); // default to wall type

                    patchEntries.add(e.keyword(), dummyEntries);
                }
            }
            else
            {
                Info<< "Ignoring unknown patch type " << type << endl;
            }
        }

        return false;
    }
}


Foam::autoPtr<Foam::fvMesh> Foam::dummyFvMesh::equivalent1DMesh
(
    const Time& runTime
)
{
    DebugInfo << "Constructing 1-D mesh" << nl << endl;

    const word instance =
        runTime.findInstance
        (
            polyMesh::meshSubDir,
            "boundary",
            IOobject::READ_IF_PRESENT
        );

    // Read patches - filter out proc patches for parallel runs
    dictionary patchEntries;
    label nPatchWithFace = 0;
    bool createFromMesh =
        setPatchEntries(runTime, instance, patchEntries, nPatchWithFace);

    const label nPatch = patchEntries.size();

    DebugPout << "Read " << nPatch << " patches" << endl;

    /*
               3 ---- 7
          f5   |\     |\   f2
          |    | 2 ---- 6   \
          |    0 |--- 4 |    \
          |     \|     \|    f3
          f4     1 ---- 5

                f0 ----- f1
    */

    // Create mesh with nPatchWithFace hex cells
    // - hex face 1: set internal face to link to next hex cell
    // - hex face 2,3,4,5: set to boundary condition
    //
    // TODO: Coupled patches:
    // - separated/collocated cyclic - faces 2 and 3
    // - rotational cyclic - convert to separated and warn not handled?
    // - proc patches - faces need to be collocated?

    pointField points(nPatchWithFace*4 + 4);
    faceList faces(nPatchWithFace*5 + 1);

    labelList owner(faces.size(), label(-1));
    labelList neighbour(owner);

    vector dx(Zero);
    {
        // Determine the mesh bounds
        IOobject pointsIO
        (
            "points",
            instance,
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        scalar dxi = 0;
        scalar dyi = 0;
        scalar dzi = 0;
        point origin(point::zero);

        if (pointsIO.typeHeaderOk<vectorIOField>(true))
        {
            const pointIOField meshPoints(pointsIO);

            boundBox meshBb(meshPoints, true);

            Info<< "Mesh bounds: " << meshBb << endl;

            origin = meshBb.min();
            vector span = meshBb.span();
            dxi = span.x()/scalar(nPatchWithFace);
            dyi = span.y();
            dzi = span.z();
        }
        else
        {
            scalar Lref = GREAT;
            origin = point(-Lref, -Lref, -Lref);
            dxi = 2.0*Lref/scalar(nPatchWithFace);
            dyi = Lref;
            dzi = Lref;
        }

        dx = vector(dxi, 0, 0);
        const vector dy(0, dyi, 0);
        const vector dz(0, 0, dzi);

        // End face
        points[0] = origin;
        points[1] = origin + dy;
        points[2] = origin + dy + dz;
        points[3] = origin + dz;
    }

    label n = 4;
    for (label i = 1; i <= nPatchWithFace; ++i)
    {
        const vector idx(i*dx);
        points[i*n] = points[0] + idx;
        points[i*n + 1] = points[1] + idx;
        points[i*n + 2] = points[2] + idx;
        points[i*n + 3] = points[3] + idx;
    }

    if (debug) Pout<< "points:" << points << endl;

    label facei = 0;

    // Internal faces first
    for (label i = 0; i < nPatchWithFace - 1; ++i)
    {
        label o = i*4;
        faces[facei] = face({4 + o, 5 + o, 6 + o, 7 + o});
        owner[facei] = i;
        neighbour[facei] = i + 1;
        ++facei;
    }

    // Boundary conditions
    for (label i = 0; i < nPatchWithFace; ++i)
    {
        label o = i*4;
        faces[facei] = face({0 + o, 4 + o, 7 + o, 3 + o});
        owner[facei] = i;
        ++facei;

        faces[facei] = face({0 + o, 1 + o, 5 + o, 4 + o});
        owner[facei] = i;
        ++facei;

        faces[facei] = face({1 + o, 2 + o, 6 + o, 5 + o});
        owner[facei] = i;
        ++facei;

        faces[facei] = face({3 + o, 7 + o, 6 + o, 2 + o});
        owner[facei] = i;
        ++facei;
    }
    {
        // End caps
        faces[facei] = face({0, 3, 2, 1});
        owner[facei] = 0;
        ++facei;

        label o = 4*nPatchWithFace;
        faces[facei] = face({0 + o, 1 + o, 2 + o, 3 + o});
        owner[facei] = nPatchWithFace - 1;
        ++facei;
    }

    DebugPout
        << "faces:" << faces << nl
        << "owner:" << owner << nl
        << "neighbour:" << neighbour
        << endl;

    autoPtr<fvMesh> meshPtr
    (
        new Foam::fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.constant(),
                runTime,
                IOobject::NO_READ, // Do not read any existing mesh
                IOobject::NO_WRITE
            ),
            std::move(points),
            std::move(faces),
            std::move(owner),
            std::move(neighbour)
        )
    );

    Foam::fvMesh& mesh = meshPtr();

    // Workaround to read fvSchemes and fvSolution
    {
        mesh.fvSchemes::readOpt() = IOobject::MUST_READ;
        mesh.fvSchemes::read();
        mesh.fvSolution::readOpt() = IOobject::MUST_READ;
        mesh.fvSolution::read();
    }

    List<polyPatch*> patches(nPatch + 1);

    label nInternalFace = nPatchWithFace - 1;
    label startFace = nInternalFace;
    label entryi = 0;
    for (const entry& e : patchEntries)
    {
        // Re-create boundary types, but reset nFaces and startFace settings
        dictionary patchDict = e.dict();
        const word& patchName = e.keyword();

        DebugPout << "Setting " << patchName << endl;

        label nFaces0 = readLabel(patchDict.lookup("nFaces"));

        if (nFaces0)
        {
            // Only set to 4 faces if there were faces in the original patch
            nFaces0 = 4;
            patchDict.set("nFaces", nFaces0);
        }

        patchDict.set("startFace", startFace);
        patches[entryi] =
            polyPatch::New
            (
                patchName,
                patchDict,
                entryi,
                mesh.boundaryMesh()
            ).ptr();

        ++entryi;
        startFace += nFaces0;
    }

    patches.last() = new emptyPolyPatch
    (
        typeName + ":default",              // name
        2,                                  // number of faces
        nInternalFace + 4*nPatchWithFace,   // start face
        nPatch - 1,                         // index in boundary list
        mesh.boundaryMesh(),                // polyBoundaryMesh
        emptyPolyPatch::typeName            // patchType
    );

    mesh.addFvPatches(patches);

    if (debug)
    {
        Pout<< "patches:" << nl << endl;
        forAll(patches, patchi)
        {
            Pout<< "patch: " << patches[patchi]->name() << nl
                << *patches[patchi] << endl;
        }
    }

    if (createFromMesh)
    {
        // Initialise the zones
        initialiseZone<pointZoneMesh>("point", instance, mesh.pointZones());
        initialiseZone<faceZoneMesh>("face", instance, mesh.faceZones());
        initialiseZone<cellZoneMesh>("cell", instance, mesh.cellZones());

        // Dummy sets created on demand
        topoSet::disallowGenericSets = 1;
    }
    else
    {
        // Dummy zones and sets created on demand
        cellZoneMesh::disallowGenericZones = 1;
        topoSet::disallowGenericSets = 1;
    }

    if (debug)
    {
        mesh.setInstance(runTime.timeName());
        mesh.objectRegistry::write();
    }

    return meshPtr;
}


// ************************************************************************* //
