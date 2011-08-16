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

Description
    Extrude faceZones into separate mesh (as a different region).

    - used to e.g. extrude baffles (extrude internal faces) or create
    liquid film regions.
    - if extruding internal faces:
        - create baffles in original mesh with directMappedWall patches
    - if extruding boundary faces:
        - convert boundary faces to directMappedWall patches
    - extrude edges of faceZone as a \<zone\>_sidePatch
    - extrude edges inbetween different faceZones as a
      (nonuniformTransform)cyclic \<zoneA\>_\<zoneB\>
    - extrudes into master direction (i.e. away from the owner cell
      if flipMap is false)
    - not parallel


Internal face extrusion
-----------------------

    +-------------+
    |             |
    |             |
    +---AAAAAAA---+
    |             |
    |             |
    +-------------+

    AAA=faceZone to extrude.


For the case of no flipMap the extrusion starts at owner and extrudes
into the space of the neighbour:

      +CCCCCCC+
      |       |         <= extruded mesh
      +BBBBBBB+

    +-------------+
    |             |
    | (neighbour) |
    |___CCCCCCC___|       <= original mesh (with 'baffles' added)
    |   BBBBBBB   |
    |(owner side) |
    |             |
    +-------------+

    BBB=directMapped between owner on original mesh and new extrusion.
        (zero offset)
    CCC=directMapped between neighbour on original mesh and new extrusion
        (offset due to the thickness of the extruded mesh)

For the case of flipMap the extrusion is the other way around: from the
neighbour side into the owner side.


Boundary face extrusion
-----------------------

    +--AAAAAAA--+
    |           |
    |           |
    +-----------+

    AAA=faceZone to extrude. E.g. slave side is owner side (no flipmap)

becomes

      +CCCCCCC+
      |       |         <= extruded mesh
      +BBBBBBB+

    +--BBBBBBB--+
    |           |       <= original mesh
    |           |
    +-----------+

    BBB=directMapped between original mesh and new extrusion
    CCC=polypatch




Usage

    - extrudeToRegionMesh \<regionName\> \<faceZones\> \<thickness\>

    \param \<regionName\> \n
    Name of mesh to create.

    \param \<faceZones\> \n
    List of faceZones to extrude

    \param \<thickness\> \n
    Thickness of extruded mesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "polyTopoChange.H"
#include "patchPointEdgeCirculator.H"
#include "OFstream.H"
#include "meshTools.H"
#include "directMappedWallPolyPatch.H"
#include "createShellMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "cyclicPolyPatch.H"
#include "wedgePolyPatch.H"
#include "nonuniformTransformCyclicPolyPatch.H"
#include "extrudeModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void addPatchFields(fvMesh& mesh, const word& patchFieldType)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            const_cast<typename GeoField::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        label sz = bfld.size();
        bfld.setSize(sz+1);
        bfld.set
        (
            sz,
            GeoField::PatchFieldType::New
            (
                patchFieldType,
                mesh.boundary()[sz],
                fld.dimensionedInternalField()
            )
        );
    }
}


// Remove last patch field
template<class GeoField>
void trimPatchFields(fvMesh& mesh, const label nPatches)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        const_cast<typename GeoField::GeometricBoundaryField&>
        (
            fld.boundaryField()
        ).setSize(nPatches);
    }
}


// Reorder patch field
template<class GeoField>
void reorderPatchFields(fvMesh& mesh, const labelList& oldToNew)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            const_cast<typename GeoField::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        bfld.reorder(oldToNew);
    }
}


void addAllPatchFields(fvMesh& mesh, const label insertPatchI)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    label sz = polyPatches.size();

    addPatchFields<volScalarField>
    (
        mesh,
        calculatedFvPatchField<scalar>::typeName
    );
    addPatchFields<volVectorField>
    (
        mesh,
        calculatedFvPatchField<vector>::typeName
    );
    addPatchFields<volSphericalTensorField>
    (
        mesh,
        calculatedFvPatchField<sphericalTensor>::typeName
    );
    addPatchFields<volSymmTensorField>
    (
        mesh,
        calculatedFvPatchField<symmTensor>::typeName
    );
    addPatchFields<volTensorField>
    (
        mesh,
        calculatedFvPatchField<tensor>::typeName
    );

    // Surface fields

    addPatchFields<surfaceScalarField>
    (
        mesh,
        calculatedFvPatchField<scalar>::typeName
    );
    addPatchFields<surfaceVectorField>
    (
        mesh,
        calculatedFvPatchField<vector>::typeName
    );
    addPatchFields<surfaceSphericalTensorField>
    (
        mesh,
        calculatedFvPatchField<sphericalTensor>::typeName
    );
    addPatchFields<surfaceSymmTensorField>
    (
        mesh,
        calculatedFvPatchField<symmTensor>::typeName
    );
    addPatchFields<surfaceTensorField>
    (
        mesh,
        calculatedFvPatchField<tensor>::typeName
    );

    // Create reordering list
    // patches before insert position stay as is
    labelList oldToNew(sz);
    for (label i = 0; i < insertPatchI; i++)
    {
        oldToNew[i] = i;
    }
    // patches after insert position move one up
    for (label i = insertPatchI; i < sz-1; i++)
    {
        oldToNew[i] = i+1;
    }
    // appended patch gets moved to insert position
    oldToNew[sz-1] = insertPatchI;

    // Shuffle into place
    polyPatches.reorder(oldToNew);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);
}


// Adds patch if not yet there. Returns patchID.
template<class PatchType>
label addPatch(fvMesh& mesh, const word& patchName, const dictionary& dict)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    label patchI = polyPatches.findPatchID(patchName);
    if (patchI != -1)
    {
        if (isA<PatchType>(polyPatches[patchI]))
        {
            // Already there
            return patchI;
        }
        else
        {
            FatalErrorIn("addPatch<PatchType>(fvMesh&, const word&)")
                << "Already have patch " << patchName
                << " but of type " << PatchType::typeName
                << exit(FatalError);
        }
    }


    label insertPatchI = polyPatches.size();
    label startFaceI = mesh.nFaces();

    forAll(polyPatches, patchI)
    {
        const polyPatch& pp = polyPatches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            insertPatchI = patchI;
            startFaceI = pp.start();
            break;
        }
    }

    dictionary patchDict(dict);
    patchDict.set("type", PatchType::typeName);
    patchDict.set("nFaces", 0);
    patchDict.set("startFace", startFaceI);


    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    label sz = polyPatches.size();

    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Add polyPatch at the end
    polyPatches.setSize(sz+1);
    polyPatches.set
    (
        sz,
        polyPatch::New
        (
            patchName,
            patchDict,
            insertPatchI,
            polyPatches
        )
    );
    fvPatches.setSize(sz+1);
    fvPatches.set
    (
        sz,
        fvPatch::New
        (
            polyPatches[sz],  // point to newly added polyPatch
            mesh.boundary()
        )
    );

    addAllPatchFields(mesh, insertPatchI);

    return insertPatchI;
}


template<class PatchType>
label addPatch(fvMesh& mesh, const word& patchName)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    label patchI = polyPatches.findPatchID(patchName);
    if (patchI != -1)
    {
        if (isA<PatchType>(polyPatches[patchI]))
        {
            // Already there
            return patchI;
        }
        else
        {
            FatalErrorIn("addPatch<PatchType>(fvMesh&, const word&)")
                << "Already have patch " << patchName
                << " but of type " << PatchType::typeName
                << exit(FatalError);
        }
    }


    label insertPatchI = polyPatches.size();
    label startFaceI = mesh.nFaces();

    forAll(polyPatches, patchI)
    {
        const polyPatch& pp = polyPatches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            insertPatchI = patchI;
            startFaceI = pp.start();
            break;
        }
    }

    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    label sz = polyPatches.size();

    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Add polyPatch at the end
    polyPatches.setSize(sz+1);
    polyPatches.set
    (
        sz,
        polyPatch::New
        (
            PatchType::typeName,
            patchName,
            0,              // size
            startFaceI,
            insertPatchI,
            polyPatches
        )
    );
    fvPatches.setSize(sz+1);
    fvPatches.set
    (
        sz,
        fvPatch::New
        (
            polyPatches[sz],  // point to newly added polyPatch
            mesh.boundary()
        )
    );

    addAllPatchFields(mesh, insertPatchI);

    return insertPatchI;
}


// Reorder and delete patches.
void reorderPatches
(
    fvMesh& mesh,
    const labelList& oldToNew,
    const label nNewPatches
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Shuffle into place
    polyPatches.reorder(oldToNew);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);

    // Remove last.
    polyPatches.setSize(nNewPatches);
    fvPatches.setSize(nNewPatches);
    trimPatchFields<volScalarField>(mesh, nNewPatches);
    trimPatchFields<volVectorField>(mesh, nNewPatches);
    trimPatchFields<volSphericalTensorField>(mesh, nNewPatches);
    trimPatchFields<volSymmTensorField>(mesh, nNewPatches);
    trimPatchFields<volTensorField>(mesh, nNewPatches);
    trimPatchFields<surfaceScalarField>(mesh, nNewPatches);
    trimPatchFields<surfaceVectorField>(mesh, nNewPatches);
    trimPatchFields<surfaceSphericalTensorField>(mesh, nNewPatches);
    trimPatchFields<surfaceSymmTensorField>(mesh, nNewPatches);
    trimPatchFields<surfaceTensorField>(mesh, nNewPatches);
}


// Remove zero-sized patches
void deleteEmptyPatches(fvMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList oldToNew(patches.size());
    label usedI = 0;
    label notUsedI = patches.size();
    forAll(patches, patchI)
    {
        if (returnReduce(patches[patchI].size(), sumOp<label>()) == 0)
        {
            oldToNew[patchI] = --notUsedI;
        }
        else
        {
            oldToNew[patchI] = usedI++;
        }
    }

    reorderPatches(mesh, oldToNew, usedI);
}


void createDummyFvMeshFiles(const polyMesh& mesh, const word& regionName)
{
    // Create dummy system/fv*
    {
        IOobject io
        (
            "fvSchemes",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        Info<< "Testing:" << io.objectPath() << endl;

        if (!io.headerOk())
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.headerOk())
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
}


// Find a patch face that is not extruded. Return -1 if not found.
label findUncoveredPatchFace
(
    const fvMesh& mesh,
    const UIndirectList<label>& extrudeMeshFaces,// mesh faces that are extruded
    const label meshEdgeI                       // mesh edge
)
{
    // Make set of extruded faces.
    labelHashSet extrudeFaceSet(extrudeMeshFaces.size());
    forAll(extrudeMeshFaces, i)
    {
        extrudeFaceSet.insert(extrudeMeshFaces[i]);
    }

    const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];
    forAll(eFaces, i)
    {
        label faceI = eFaces[i];
        if (!mesh.isInternalFace(faceI) && !extrudeFaceSet.found(faceI))
        {
            return faceI;
        }
    }
    return -1;
}


// Count the number of faces in patches that need to be created
void countExtrudePatches
(
    const fvMesh& mesh,
    const primitiveFacePatch& extrudePatch,
    const label nZones,
    const labelList& zoneID,
    const labelList& extrudeMeshFaces,
    const labelList& extrudeMeshEdges,

    labelList& zoneSidePatch,
    labelList& zoneZonePatch
)
{
    const labelListList& edgeFaces = extrudePatch.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];
        if (eFaces.size() == 2)
        {
            label zone0 = zoneID[eFaces[0]];
            label zone1 = zoneID[eFaces[1]];

            if (zone0 != zone1)
            {
                label minZone = min(zone0,zone1);
                label maxZone = max(zone0,zone1);
                zoneZonePatch[minZone*nZones+maxZone]++;
            }
        }
        else
        {
            // Check whether we are on a mesh edge with external patches. If
            // so choose any uncovered one. If none found put face in
            // undetermined zone 'side' patch

            label faceI = findUncoveredPatchFace
            (
                mesh,
                UIndirectList<label>(extrudeMeshFaces, eFaces),
                extrudeMeshEdges[edgeI]
            );

            if (faceI == -1)
            {
                // Determine the min zone of all connected zones.
                label minZone = zoneID[eFaces[0]];
                for (label i = 1; i < eFaces.size(); i++)
                {
                    minZone = min(minZone, zoneID[eFaces[i]]);
                }
                zoneSidePatch[minZone]++;
            }
        }
    }
    Pstream::listCombineGather(zoneSidePatch, plusEqOp<label>());
    Pstream::listCombineScatter(zoneSidePatch);
    Pstream::listCombineGather(zoneZonePatch, plusEqOp<label>());
    Pstream::listCombineScatter(zoneZonePatch);
}


// Constrain&sync normals on points that are on coupled patches to make sure
// the face extruded from the edge has a valid normal with its coupled
// equivalent.
// Note that only points on cyclic edges need to be constrained and not
// all points touching cyclics since only edges become faces.
void constrainCoupledNormals
(
    const fvMesh& mesh,
    const primitiveFacePatch& extrudePatch,
    const labelList& meshEdges,
    const faceList& pointRegions,   // per face, per index the region

    vectorField& regionNormals
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Mark edges that are on boundary of extrusion.
    Map<label> meshToExtrudEdge
    (
        2*(extrudePatch.nEdges()-extrudePatch.nInternalEdges())
    );
    for
    (
        label extrudeEdgeI = extrudePatch.nInternalEdges();
        extrudeEdgeI < extrudePatch.nEdges();
        extrudeEdgeI++
    )
    {
        meshToExtrudEdge.insert(meshEdges[extrudeEdgeI], extrudeEdgeI);
    }


    // For owner: normal at first point of edge when walking through faces
    // in order.
    vectorField edgeNormals0(mesh.nEdges(), vector::zero);
    vectorField edgeNormals1(mesh.nEdges(), vector::zero);

    // Loop through all edges of patch. If they are to be extruded mark the
    // point normals in order.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<cyclicPolyPatch>(pp))
        {
            bool isOwner = refCast<const cyclicPolyPatch>(pp).owner();

            forAll(pp.faceEdges(), faceI)
            {
                const labelList& fEdges = pp.faceEdges()[faceI];
                forAll(fEdges, fp)
                {
                    label meshEdgeI = pp.meshEdges()[fEdges[fp]];
                    if (meshToExtrudEdge.found(meshEdgeI))
                    {
                        // Edge corresponds to a extrusion edge. Store extrusion
                        // normals on edge so we can syncTools it.

                        //const edge& ppE = pp.edges()[fEdges[fp]];
                        //Pout<< "ppedge:" << pp.localPoints()[ppE[0]]
                        //    << pp.localPoints()[ppE[1]]
                        //    << endl;

                        const face& f = pp.localFaces()[faceI];
                        label fp0 = fp;
                        label fp1 = f.fcIndex(fp0);
                        label mp0 = pp[faceI][fp0];
                        label mp1 = pp[faceI][fp1];

                        // Find corresponding face and indices.
                        vector regionN0;
                        vector regionN1;
                        {
                            label exEdgeI = meshToExtrudEdge[meshEdgeI];
                            const labelList& eFaces =
                                extrudePatch.edgeFaces()[exEdgeI];
                            // Use 0th face.
                            label exFaceI = eFaces[0];
                            const face& exF = extrudePatch[exFaceI];
                            const face& exRegions = pointRegions[exFaceI];
                            // Find points
                            label r0 = exRegions[findIndex(exF, mp0)];
                            regionN0 = regionNormals[r0];
                            label r1 = exRegions[findIndex(exF, mp1)];
                            regionN1 = regionNormals[r1];
                        }

                        vector& nA =
                        (
                            isOwner
                          ? edgeNormals0[meshEdgeI]
                          : edgeNormals1[meshEdgeI]
                        );

                        nA = regionN0;
                        const vector& cyc0 = pp.pointNormals()[f[fp0]];
                        nA -= (nA&cyc0)*cyc0;

                        vector& nB =
                        (
                            isOwner
                          ? edgeNormals1[meshEdgeI]
                          : edgeNormals0[meshEdgeI]
                        );

                        nB = regionN1;
                        const vector& cyc1 = pp.pointNormals()[f[fp1]];
                        nB -= (nB&cyc1)*cyc1;
                    }
                }
            }
        }
    }


    // Synchronise regionNormals
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Synchronise
    syncTools::syncEdgeList
    (
        mesh,
        edgeNormals0,
        maxMagSqrEqOp<vector>(),
        vector::zero            // nullValue
    );
    syncTools::syncEdgeList
    (
        mesh,
        edgeNormals1,
        maxMagSqrEqOp<vector>(),
        vector::zero            // nullValue
    );


    // Re-work back into regionNormals
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<cyclicPolyPatch>(pp))
        {
            bool isOwner = refCast<const cyclicPolyPatch>(pp).owner();

            forAll(pp.faceEdges(), faceI)
            {
                const labelList& fEdges = pp.faceEdges()[faceI];
                forAll(fEdges, fp)
                {
                    label meshEdgeI = pp.meshEdges()[fEdges[fp]];
                    if (meshToExtrudEdge.found(meshEdgeI))
                    {
                        const face& f = pp.localFaces()[faceI];
                        label fp0 = fp;
                        label fp1 = f.fcIndex(fp0);
                        label mp0 = pp[faceI][fp0];
                        label mp1 = pp[faceI][fp1];


                        const vector& nA =
                        (
                            isOwner
                          ? edgeNormals0[meshEdgeI]
                          : edgeNormals1[meshEdgeI]
                        );

                        const vector& nB =
                        (
                            isOwner
                          ? edgeNormals1[meshEdgeI]
                          : edgeNormals0[meshEdgeI]
                        );

                        // Find corresponding face and indices.
                        {
                            label exEdgeI = meshToExtrudEdge[meshEdgeI];
                            const labelList& eFaces =
                                extrudePatch.edgeFaces()[exEdgeI];
                            // Use 0th face.
                            label exFaceI = eFaces[0];
                            const face& exF = extrudePatch[exFaceI];
                            const face& exRegions = pointRegions[exFaceI];
                            // Find points
                            label r0 = exRegions[findIndex(exF, mp0)];
                            regionNormals[r0] = nA;
                            label r1 = exRegions[findIndex(exF, mp1)];
                            regionNormals[r1] = nB;
                        }
                    }
                }
            }
        }
    }
}


tmp<pointField> calcOffset
(
    const primitiveFacePatch& extrudePatch,
    const createShellMesh& extruder,
    const polyPatch& pp
)
{
    vectorField::subField fc = pp.faceCentres();

    tmp<pointField> toffsets(new pointField(fc.size()));
    pointField& offsets = toffsets();

    forAll(fc, i)
    {
        label meshFaceI = pp.start()+i;
        label patchFaceI = mag(extruder.faceToFaceMap()[meshFaceI])-1;
        point patchFc = extrudePatch[patchFaceI].centre
        (
            extrudePatch.points()
        );
        offsets[i] = patchFc - fc[i];
    }
    return toffsets;
}



// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addNote
    (
        "Create region mesh by extruding a faceZone"
    );

    #include "addRegionOption.H"
    #include "addOverwriteOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();
    bool overwrite = args.optionFound("overwrite");

    IOdictionary dict
    (
        IOobject
        (
            "extrudeToRegionMeshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );

    // Point generator
    autoPtr<extrudeModel> model(extrudeModel::New(dict));


    // Region
    const word shellRegionName(dict.lookup("region"));
    const wordList zoneNames(dict.lookup("faceZones"));
    wordList zoneShadowNames(0);
    if (dict.found("faceZonesShadow"))
    {
        dict.lookup("faceZonesShadow") >> zoneShadowNames;
    }

    const Switch oneD(dict.lookup("oneD"));
    const Switch adaptMesh(dict.lookup("adaptMesh"));

    Info<< "Extruding zones " << zoneNames
        << " on mesh " << regionName
        << " into shell mesh " << shellRegionName
        << endl;

    if (shellRegionName == regionName)
    {
        FatalErrorIn(args.executable())
            << "Cannot extrude into same region as mesh." << endl
            << "Mesh region : " << regionName << endl
            << "Shell region : " << shellRegionName
            << exit(FatalError);
    }



    // Create dummy fv* files
    createDummyFvMeshFiles(mesh, shellRegionName);


    word meshInstance;
    if (!overwrite)
    {
        runTime++;
        meshInstance = runTime.timeName();
    }
    else
    {
        meshInstance = oldInstance;
    }
    Info<< "Writing meshes to " << meshInstance << nl << endl;


    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const faceZoneMesh& faceZones = mesh.faceZones();


    // Check zones
    // ~~~~~~~~~~~

    labelList zoneIDs(zoneNames.size());
    forAll(zoneNames, i)
    {
        zoneIDs[i] = faceZones.findZoneID(zoneNames[i]);
        if (zoneIDs[i] == -1)
        {
            FatalErrorIn(args.executable())
                << "Cannot find zone " << zoneNames[i] << endl
                << "Valid zones are " << faceZones.names()
                << exit(FatalError);
        }
    }

    labelList zoneShadowIDs;
    if (zoneShadowNames.size())
    {
        zoneShadowIDs.setSize(zoneShadowNames.size());
        forAll(zoneShadowNames, i)
        {
            zoneShadowIDs[i] = faceZones.findZoneID(zoneShadowNames[i]);
            if (zoneShadowIDs[i] == -1)
            {
                FatalErrorIn(args.executable())
                    << "Cannot find zone " << zoneShadowNames[i] << endl
                    << "Valid zones are " << faceZones.names()
                    << exit(FatalError);
            }
        }
    }

    label nShadowFaces = 0;
    forAll(zoneShadowIDs, i)
    {
        nShadowFaces += faceZones[zoneShadowIDs[i]].size();
    }

    labelList extrudeMeshShadowFaces(nShadowFaces);
    boolList zoneShadowFlipMap(nShadowFaces);
    labelList zoneShadowID(nShadowFaces);

    nShadowFaces = 0;
    forAll(zoneShadowIDs, i)
    {
        const faceZone& fz = faceZones[zoneShadowIDs[i]];
        forAll(fz, j)
        {
            extrudeMeshShadowFaces[nShadowFaces] = fz[j];
            zoneShadowFlipMap[nShadowFaces] = fz.flipMap()[j];
            zoneShadowID[nShadowFaces] = zoneShadowIDs[i];
            nShadowFaces++;
        }
    }



    // Collect faces to extrude and per-face information
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nExtrudeFaces = 0;
    forAll(zoneIDs, i)
    {
        nExtrudeFaces += faceZones[zoneIDs[i]].size();
    }
    labelList extrudeMeshFaces(nExtrudeFaces);
    faceList zoneFaces(nExtrudeFaces);
    labelList zoneID(nExtrudeFaces);
    boolList zoneFlipMap(nExtrudeFaces);
    nExtrudeFaces = 0;
    forAll(zoneIDs, i)
    {
        const faceZone& fz = faceZones[zoneIDs[i]];
        const primitiveFacePatch& fzp = fz();
        forAll(fz, j)
        {
            extrudeMeshFaces[nExtrudeFaces] = fz[j];
            zoneFaces[nExtrudeFaces] = fzp[j];
            zoneID[nExtrudeFaces] = zoneIDs[i];
            zoneFlipMap[nExtrudeFaces] = fz.flipMap()[j];
            nExtrudeFaces++;
        }
    }
    primitiveFacePatch extrudePatch(zoneFaces.xfer(), mesh.points());
    const pointField& extrudePoints = extrudePatch.localPoints();
    const faceList& extrudeFaces = extrudePatch.localFaces();
    const labelListList& edgeFaces = extrudePatch.edgeFaces();


    Info<< "extrudePatch :"
        << " faces:" << extrudePatch.size()
        << " points:" << extrudePatch.nPoints()
        << " edges:" << extrudePatch.nEdges()
        << nl
        << endl;

     // Check nExtrudeFaces = nShadowFaces
    if (zoneShadowNames.size())
    {
        if (nExtrudeFaces != nShadowFaces)
        {
            FatalErrorIn(args.executable())
                << "Extruded faces " << nExtrudeFaces << endl
                << "is different from shadow faces. " << nShadowFaces
                << "This is not permitted " << endl
                << exit(FatalError);
        }
    }


    // Determine corresponding mesh edges
    const labelList extrudeMeshEdges
    (
        extrudePatch.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );


    // Check whether the zone is internal or external faces to determine
    // what patch type to insert. Cannot be mixed
    // since then how to couple? - directMapped only valid for a whole patch.
    boolList isInternal(zoneIDs.size(), false);
    forAll(zoneIDs, i)
    {
        const faceZone& fz = faceZones[zoneIDs[i]];
        forAll(fz, j)
        {
            if (mesh.isInternalFace(fz[j]))
            {
                isInternal[i] = true;
                break;
            }
        }
    }
    Pstream::listCombineGather(isInternal, orEqOp<bool>());
    Pstream::listCombineScatter(isInternal);

    forAll(zoneIDs, i)
    {
        const faceZone& fz = faceZones[zoneIDs[i]];
        if (isInternal[i])
        {
            Info<< "FaceZone " << fz.name() << " has internal faces" << endl;
        }
        else
        {
            Info<< "FaceZone " << fz.name() << " has boundary faces" << endl;
        }
    }
    Info<< endl;



    // Add interface patches
    // ~~~~~~~~~~~~~~~~~~~~~
    // Note that these actually get added to the original mesh
    // so the shell mesh creation copies them. They then get removed
    // from the original mesh.

    Info<< "Adding coupling patches:" << nl << nl
        << "patchID\tpatch\ttype" << nl
        << "-------\t-----\t----"
        << endl;
    labelList interRegionTopPatch(zoneNames.size());
    labelList interRegionBottomPatch(zoneNames.size());
    label nCoupled = 0;
    forAll(zoneIDs, i)
    {
        word interName(regionName+"_to_"+shellRegionName+'_'+zoneNames[i]);

        if (isInternal[i])
        {
            interRegionTopPatch[i] = addPatch<directMappedWallPolyPatch>
            (
                mesh,
                interName + "_top"
            );
            nCoupled++;
            Info<< interRegionTopPatch[i]
                << '\t' << patches[interRegionTopPatch[i]].name()
                << '\t' << patches[interRegionTopPatch[i]].type()
                << nl;

            interRegionBottomPatch[i] = addPatch<directMappedWallPolyPatch>
            (
                mesh,
                interName + "_bottom"
            );
            nCoupled++;
            Info<< interRegionBottomPatch[i]
                << '\t' << patches[interRegionBottomPatch[i]].name()
                << '\t' << patches[interRegionBottomPatch[i]].type()
                << nl;
        }
        else if (zoneShadowNames.size() == 0)
        {
            interRegionTopPatch[i] = addPatch<polyPatch>
            (
                mesh,
                zoneNames[i] + "_top"
            );
            nCoupled++;
            Info<< interRegionTopPatch[i]
                << '\t' << patches[interRegionTopPatch[i]].name()
                << '\t' << patches[interRegionTopPatch[i]].type()
                << nl;

            interRegionBottomPatch[i] = addPatch<directMappedWallPolyPatch>
            (
                mesh,
                interName
            );
            nCoupled++;
            Info<< interRegionBottomPatch[i]
                << '\t' << patches[interRegionBottomPatch[i]].name()
                << '\t' << patches[interRegionBottomPatch[i]].type()
                << nl;
        }
        else if (zoneShadowNames.size() > 0) //patch using shadow face zones.
        {
            interRegionTopPatch[i] = addPatch<directMappedWallPolyPatch>
            (
                mesh,
                zoneShadowNames[i] + "_top"
            );
            nCoupled++;
            Info<< interRegionTopPatch[i]
                << '\t' << patches[interRegionTopPatch[i]].name()
                << '\t' << patches[interRegionTopPatch[i]].type()
                << nl;

            interRegionBottomPatch[i] = addPatch<directMappedWallPolyPatch>
            (
                mesh,
                interName
            );
            nCoupled++;
            Info<< interRegionBottomPatch[i]
                << '\t' << patches[interRegionBottomPatch[i]].name()
                << '\t' << patches[interRegionBottomPatch[i]].type()
                << nl;
        }

    }
    Info<< "Added " << nCoupled << " inter-region patches." << nl
        << endl;


    labelList extrudeTopPatchID(extrudePatch.size());
    labelList extrudeBottomPatchID(extrudePatch.size());

    nExtrudeFaces = 0;
    forAll(zoneNames, i)
    {
        const faceZone& fz = faceZones[zoneNames[i]];
        forAll(fz, j)
        {
            extrudeTopPatchID[nExtrudeFaces] = interRegionTopPatch[i];
            extrudeBottomPatchID[nExtrudeFaces] = interRegionBottomPatch[i];
            nExtrudeFaces++;
        }
    }



    // Count how many patches on special edges of extrudePatch are necessary
    // - zoneXXX_sides
    // - zoneXXX_zoneYYY
    labelList zoneSidePatch(faceZones.size(), 0);
    // Patch to use for minZone
    labelList zoneZonePatch_min(faceZones.size()*faceZones.size(), 0);
    // Patch to use for maxZone
    labelList zoneZonePatch_max(faceZones.size()*faceZones.size(), 0);

    countExtrudePatches
    (
        mesh,
        extrudePatch,
        faceZones.size(),
        zoneID,
        extrudeMeshFaces,
        extrudeMeshEdges,

        zoneSidePatch,      // reuse for counting
        zoneZonePatch_min   // reuse for counting
    );

    // Now check which patches to add.
    Info<< "Adding patches for edges on zones:" << nl << nl
        << "patchID\tpatch" << nl
        << "-------\t-----"
        << endl;

    label nSide = 0;

    forAll(zoneSidePatch, zoneI)
    {
        if (oneD)
        {
            // Reuse single empty patch.
            word patchType = dict.lookup("oneDPolyPatchType");
            word patchName;
            if (patchType == "emptyPolyPatch")
            {
                patchName = "oneDEmptyPatch";
                zoneSidePatch[zoneI] = addPatch<emptyPolyPatch>
                (
                    mesh,
                    patchName
                );
            }
            else if (patchType == "wedgePolyPatch")
            {
                patchName = "oneDWedgePatch";
                zoneSidePatch[zoneI] = addPatch<wedgePolyPatch>
                (
                    mesh,
                    patchName
                );
            }
            else
            {
                FatalErrorIn(args.executable())
                << "Type " << patchType << " does not exist "
                << exit(FatalError);
            }

            Info<< zoneSidePatch[zoneI] << '\t' << patchName << nl;

            nSide++;
        }
        else if (zoneSidePatch[zoneI] > 0)
        {
            word patchName = faceZones[zoneI].name() + "_" + "side";

            zoneSidePatch[zoneI] = addPatch<polyPatch>
            (
                mesh,
                patchName
            );

            Info<< zoneSidePatch[zoneI] << '\t' << patchName << nl;

            nSide++;
        }
    }
    Info<< "Added " << nSide << " zone-edge patches." << nl
        << endl;



    Info<< "Adding inter-zone patches:" << nl << nl
        << "patchID\tpatch" << nl
        << "-------\t-----"
        << endl;

    dictionary transformDict;
    transformDict.add
    (
        "transform",
        cyclicPolyPatch::transformTypeNames[cyclicPolyPatch::NOORDERING]
    );

    label nInter = 0;
    if (!oneD)
    {
        forAll(zoneZonePatch_min, minZone)
        {
            for (label maxZone = minZone; maxZone < faceZones.size(); maxZone++)
            {
                label index = minZone*faceZones.size()+maxZone;

                if (zoneZonePatch_min[index] > 0)
                {
                    word minToMax =
                        faceZones[minZone].name()
                      + "_to_"
                      + faceZones[maxZone].name();
                    word maxToMin =
                        faceZones[maxZone].name()
                      + "_to_"
                      + faceZones[minZone].name();

                    {
                        transformDict.set("neighbourPatch", maxToMin);
                        zoneZonePatch_min[index] =
                        addPatch<nonuniformTransformCyclicPolyPatch>
                        (
                            mesh,
                            minToMax,
                            transformDict
                        );
                        Info<< zoneZonePatch_min[index] << '\t' << minToMax
                            << nl;
                        nInter++;
                    }
                    {
                        transformDict.set("neighbourPatch", minToMax);
                        zoneZonePatch_max[index] =
                        addPatch<nonuniformTransformCyclicPolyPatch>
                        (
                            mesh,
                            maxToMin,
                            transformDict
                        );
                        Info<< zoneZonePatch_max[index] << '\t' << maxToMin
                            << nl;
                        nInter++;
                    }

                }
            }
        }
    }
    Info<< "Added " << nInter << " inter-zone patches." << nl
        << endl;



    // Set patches to use for edges to be extruded into boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // In order of edgeFaces: per edge, per originating face the
    // patch to use for the side face (from the extruded edge).
    // If empty size create an internal face.
    labelListList extrudeEdgePatches(extrudePatch.nEdges());

    // Is edge an non-manifold edge
    PackedBoolList nonManifoldEdge(extrudePatch.nEdges());

    // Note: logic has to be same as in countExtrudePatches.
    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        labelList& ePatches = extrudeEdgePatches[edgeI];

        if (oneD)
        {
            //nonManifoldEdge[edgeI] = 1; //To fill the space
            ePatches.setSize(eFaces.size());
            forAll(eFaces, i)
            {
                ePatches[i] = zoneSidePatch[zoneID[eFaces[i]]];
            }
            if (eFaces.size() != 2)
            {
                nonManifoldEdge[edgeI] = 1;
            }
        }
        else if (eFaces.size() == 2)
        {
            label zone0 = zoneID[eFaces[0]];
            label zone1 = zoneID[eFaces[1]];

            if (zone0 != zone1) // || (cos(angle) > blabla))
            {
                label minZone = min(zone0,zone1);
                label maxZone = max(zone0,zone1);
                label index = minZone*faceZones.size()+maxZone;

                ePatches.setSize(eFaces.size());

                if (zone0 == minZone)
                {
                    ePatches[0] = zoneZonePatch_min[index];
                    ePatches[1] = zoneZonePatch_max[index];
                }
                else
                {
                    ePatches[0] = zoneZonePatch_max[index];
                    ePatches[1] = zoneZonePatch_min[index];
                }

                nonManifoldEdge[edgeI] = 1;
            }
        }
        else
        {
            label faceI = findUncoveredPatchFace
            (
                mesh,
                UIndirectList<label>(extrudeMeshFaces, eFaces),
                extrudeMeshEdges[edgeI]
            );

            if (faceI != -1)
            {
                label patchI = mesh.boundaryMesh().whichPatch(faceI);
                ePatches.setSize(eFaces.size(), patchI);
            }
            else
            {
                ePatches.setSize(eFaces.size());
                forAll(eFaces, i)
                {
                    ePatches[i] = zoneSidePatch[zoneID[eFaces[i]]];
                }
            }
            nonManifoldEdge[edgeI] = 1;
        }
    }



    // Assign point regions
    // ~~~~~~~~~~~~~~~~~~~~

    // Per face, per point the region number.
    faceList pointRegions(extrudePatch.size());
    // Per region the originating point
    labelList regionPoints;

    createShellMesh::calcPointRegions
    (
        extrudePatch,
        nonManifoldEdge,

        pointRegions,
        regionPoints
    );


    // Calculate a normal per region
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vectorField regionNormals(regionPoints.size(), vector::zero);
    vectorField regionCentres(regionPoints.size(), vector::zero);
    labelList nRegionFaces(regionPoints.size(), 0);

    forAll(pointRegions, faceI)
    {
        const face& f = extrudeFaces[faceI];

        forAll(f, fp)
        {
            label region = pointRegions[faceI][fp];
            regionNormals[region] += extrudePatch.faceNormals()[faceI];
            regionCentres[region] += extrudePatch.faceCentres()[faceI];
            nRegionFaces[region]++;
        }
    }
    regionNormals /= mag(regionNormals);
    forAll(regionCentres, regionI)
    {
        regionCentres[regionI] /= nRegionFaces[regionI];
    }


    // Constrain&sync normals on points that are on coupled patches.
    constrainCoupledNormals
    (
        mesh,
        extrudePatch,
        extrudeMeshEdges,
        pointRegions,

        regionNormals
    );

    // For debugging: dump hedgehog plot of normals
    if (false)
    {
        OFstream str(runTime.path()/"regionNormals.obj");
        label vertI = 0;

        scalar thickness = model().sumThickness(1);

        forAll(pointRegions, faceI)
        {
            const face& f = extrudeFaces[faceI];

            forAll(f, fp)
            {
                label region = pointRegions[faceI][fp];
                const point& pt = extrudePoints[f[fp]];

                meshTools::writeOBJ(str, pt);
                vertI++;
                meshTools::writeOBJ(str, pt+thickness*regionNormals[region]);
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Use model to create displacements of first layer
    vectorField firstDisp(regionNormals.size());
    forAll(firstDisp, regionI)
    {
        const point& regionPt = regionCentres[regionI];
        const vector& n = regionNormals[regionI];
        firstDisp[regionI] = model()(regionPt, n, 1) - regionPt;
    }



    // Create a new mesh
    // ~~~~~~~~~~~~~~~~~

    createShellMesh extruder(extrudePatch, pointRegions, regionPoints);


    autoPtr<fvMesh> regionMeshPtr;
    autoPtr<mapPolyMesh> shellMap;
    {
        polyTopoChange meshMod(mesh.boundaryMesh().size());

        extruder.setRefinement
        (
            firstDisp,                              // first displacement
            model().expansionRatio(),
            model().nLayers(),                      // nLayers
            extrudeTopPatchID,
            extrudeBottomPatchID,
            extrudeEdgePatches,
            meshMod
        );

        shellMap = meshMod.makeMesh
        (
            regionMeshPtr,     // mesh to create
            IOobject
            (
                shellRegionName,
                meshInstance,
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh            // mesh to get original patch info from
        );
    }
    fvMesh& regionMesh = regionMeshPtr();

    // Necessary?
    regionMesh.setInstance(meshInstance);


    // Update numbering on extruder.
    extruder.updateMesh(shellMap);


    // Change top and bottom boundary conditions on regionMesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Save offsets from shell mesh back to original mesh
    List<pointField> topOffsets(zoneIDs.size());
    List<pointField> bottomOffsets(zoneIDs.size());

    {
        const polyBoundaryMesh& regionPatches = regionMesh.boundaryMesh();
        List<polyPatch*> newPatches(regionPatches.size());
        forAll(regionPatches, patchI)
        {
            const polyPatch& pp = regionPatches[patchI];

            if
            (
                isA<directMappedWallPolyPatch>(pp)
             && (findIndex(interRegionTopPatch, patchI) != -1)
            )
            {
                label index = findIndex(interRegionTopPatch, patchI);

                topOffsets[index] = calcOffset(extrudePatch, extruder, pp);

                newPatches[patchI] = new directMappedWallPolyPatch
                (
                    pp.name(),
                    pp.size(),
                    pp.start(),
                    patchI,
                    regionName,                             // sampleRegion
                    directMappedPatchBase::NEARESTPATCHFACE,// sampleMode
                    pp.name(),                              // samplePatch
                    topOffsets[index],                      // offset
                    patches
                );
            }
            else if
            (
                isA<directMappedWallPolyPatch>(pp)
             && (findIndex(interRegionBottomPatch, patchI) != -1)
            )
            {
                label index = findIndex(interRegionBottomPatch, patchI);

                bottomOffsets[index] = calcOffset(extrudePatch, extruder, pp);

                newPatches[patchI] = new directMappedWallPolyPatch
                (
                    pp.name(),
                    pp.size(),
                    pp.start(),
                    patchI,
                    regionName,                             // sampleRegion
                    directMappedPatchBase::NEARESTPATCHFACE,// sampleMode
                    pp.name(),                              // samplePatch
                    bottomOffsets[index],                   // offset
                    patches
                );
            }
            else
            {
                newPatches[patchI] = pp.clone
                (
                    regionPatches,
                    patchI,
                    pp.size(),
                    pp.start()
                ).ptr();
            }
        }
        regionMesh.removeFvBoundary();
        regionMesh.addFvPatches(newPatches, true);
        deleteEmptyPatches(regionMesh);
    }



    // Write addressing from region mesh back to originating patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelIOList cellToPatchFaceAddressing
    (
        IOobject
        (
            "cellToPatchFaceAddressing",
            regionMesh.facesInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        extruder.cellToFaceMap()
    );
    cellToPatchFaceAddressing.note() = "cell to patch face addressing";

    labelIOList faceToPatchFaceAddressing
    (
        IOobject
        (
            "faceToPatchFaceAddressing",
            regionMesh.facesInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        extruder.faceToFaceMap()
    );
    faceToPatchFaceAddressing.note() =
        "front/back face + turning index to patch face addressing";

    labelIOList faceToPatchEdgeAddressing
    (
        IOobject
        (
            "faceToPatchEdgeAddressing",
            regionMesh.facesInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        extruder.faceToEdgeMap()
    );
    faceToPatchEdgeAddressing.note() =
        "side face to patch edge addressing";

    labelIOList pointToPatchPointAddressing
    (
        IOobject
        (
            "pointToPatchPointAddressing",
            regionMesh.facesInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        extruder.pointToPointMap()
    );
    pointToPatchPointAddressing.note() =
        "point to patch point addressing";


    Info<< "Writing mesh " << regionMesh.name()
        << " to " << regionMesh.facesInstance() << nl
        << endl;

    bool ok =
        regionMesh.write()
     && cellToPatchFaceAddressing.write()
     && faceToPatchFaceAddressing.write()
     && faceToPatchEdgeAddressing.write()
     && pointToPatchPointAddressing.write();

    if (!ok)
    {
        FatalErrorIn(args.executable())
            << "Failed writing mesh " << regionMesh.name()
            << " at location " << regionMesh.facesInstance()
            << exit(FatalError);
    }




    // Insert baffles into original mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<mapPolyMesh> addBafflesMap;

    if (adaptMesh)
    {
        polyTopoChange meshMod(mesh);

        // Modify faces to be in bottom (= always coupled) patch
        forAll(extrudeMeshFaces, zoneFaceI)
        {
            label meshFaceI = extrudeMeshFaces[zoneFaceI];
            label zoneI = zoneID[zoneFaceI];
            bool flip = zoneFlipMap[zoneFaceI];
            const face& f = mesh.faces()[meshFaceI];

            if (!flip)
            {
                meshMod.modifyFace
                (
                    f,                          // modified face
                    meshFaceI,                  // label of face being modified
                    mesh.faceOwner()[meshFaceI],// owner
                    -1,                         // neighbour
                    false,                      // face flip
                    extrudeBottomPatchID[zoneFaceI],// patch for face
                    zoneI,                      // zone for face
                    flip                        // face flip in zone
                );
            }
            else if (mesh.isInternalFace(meshFaceI))
            {
                meshMod.modifyFace
                (
                    f.reverseFace(),                // modified face
                    meshFaceI,                      // label of modified face
                    mesh.faceNeighbour()[meshFaceI],// owner
                    -1,                             // neighbour
                    true,                           // face flip
                    extrudeBottomPatchID[zoneFaceI],// patch for face
                    zoneI,                          // zone for face
                    !flip                           // face flip in zone
                );
            }
        }

        if (zoneShadowNames.size() > 0) //if there is a top faceZone specified
        {
            forAll(extrudeMeshFaces, zoneFaceI)
            {
                label meshFaceI = extrudeMeshShadowFaces[zoneFaceI];
                label zoneI = zoneShadowID[zoneFaceI];
                bool flip = zoneShadowFlipMap[zoneFaceI];
                const face& f = mesh.faces()[meshFaceI];

                if (!flip)
                {
                    meshMod.modifyFace
                    (
                        f,                          // modified face
                        meshFaceI,                  // face being modified
                        mesh.faceOwner()[meshFaceI],// owner
                        -1,                         // neighbour
                        false,                      // face flip
                        extrudeTopPatchID[zoneFaceI],// patch for face
                        zoneI,                      // zone for face
                        flip                        // face flip in zone
                    );
                }
                else if (mesh.isInternalFace(meshFaceI))
                {
                    meshMod.modifyFace
                    (
                        f.reverseFace(),                // modified face
                        meshFaceI,                      // label modified face
                        mesh.faceNeighbour()[meshFaceI],// owner
                        -1,                             // neighbour
                        true,                           // face flip
                        extrudeTopPatchID[zoneFaceI],   // patch for face
                        zoneI,                          // zone for face
                        !flip                           // face flip in zone
                    );
                }
            }

        }
        else
        {
            // Add faces (using same points) to be in top patch
            forAll(extrudeMeshFaces, zoneFaceI)
            {
                label meshFaceI = extrudeMeshFaces[zoneFaceI];
                bool flip = zoneFlipMap[zoneFaceI];
                const face& f = mesh.faces()[meshFaceI];

                if (!flip)
                {
                    if (mesh.isInternalFace(meshFaceI))
                    {
                        meshMod.addFace
                        (
                            f.reverseFace(),                // modified face
                            mesh.faceNeighbour()[meshFaceI],// owner
                            -1,                             // neighbour
                            -1,                             // master point
                            -1,                             // master edge
                            meshFaceI,                      // master face
                            true,                           // flip flux
                            extrudeTopPatchID[zoneFaceI],   // patch for face
                            -1,                             // zone for face
                            false                           //face flip in zone
                        );
                    }
                }
                else
                {
                    meshMod.addFace
                    (
                        f,                              // face
                        mesh.faceOwner()[meshFaceI],    // owner
                        -1,                             // neighbour
                        -1,                             // master point
                        -1,                             // master edge
                        meshFaceI,                      // master face
                        false,                          // flip flux
                        extrudeTopPatchID[zoneFaceI],   // patch for face
                        -1,                             // zone for face
                        false                           // zone flip
                    );
                }
            }
        }

        // Change the mesh. Change points directly (no inflation).
        addBafflesMap = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(addBafflesMap);


//XXXXXX
// Update maps! e.g. faceToPatchFaceAddressing
//XXXXXX

        // Move mesh (since morphing might not do this)
        if (addBafflesMap().hasMotionPoints())
        {
            mesh.movePoints(addBafflesMap().preMotionPoints());
        }

        mesh.setInstance(meshInstance);
    }



    // Change master and slave boundary conditions on originating mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (adaptMesh)
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        List<polyPatch*> newPatches(patches.size());
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if
            (
                isA<directMappedWallPolyPatch>(pp)
             && (findIndex(interRegionTopPatch, patchI) != -1)
            )
            {
                label index = findIndex(interRegionTopPatch, patchI);
                newPatches[patchI] = new directMappedWallPolyPatch
                (
                    pp.name(),
                    pp.size(),
                    pp.start(),
                    patchI,
                    shellRegionName,                        // sampleRegion
                    directMappedPatchBase::NEARESTPATCHFACE,// sampleMode
                    pp.name(),                              // samplePatch
                    -topOffsets[index],                     // offset
                    patches
                );
            }
            else if
            (
                isA<directMappedWallPolyPatch>(pp)
             && (findIndex(interRegionBottomPatch, patchI) != -1)
            )
            {
                label index = findIndex(interRegionBottomPatch, patchI);

                newPatches[patchI] = new directMappedWallPolyPatch
                (
                    pp.name(),
                    pp.size(),
                    pp.start(),
                    patchI,
                    shellRegionName,                        // sampleRegion
                    directMappedPatchBase::NEARESTPATCHFACE,// sampleMode
                    pp.name(),                              // samplePatch
                    -bottomOffsets[index],                  // offset
                    patches
                );
            }
            else
            {
                newPatches[patchI] = pp.clone
                (
                    patches,
                    patchI,
                    pp.size(),
                    pp.start()
                ).ptr();
            }
        }
        mesh.removeFvBoundary();
        mesh.addFvPatches(newPatches, true);

        // Remove any unused patches
        deleteEmptyPatches(mesh);
    }


    if (adaptMesh)
    {
        Info<< "Writing mesh " << mesh.name()
            << " to " << mesh.facesInstance() << nl
            << endl;

        if (!mesh.write())
        {
            FatalErrorIn(args.executable())
                << "Failed writing mesh " << mesh.name()
                << " at location " << mesh.facesInstance()
                << exit(FatalError);
        }
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
