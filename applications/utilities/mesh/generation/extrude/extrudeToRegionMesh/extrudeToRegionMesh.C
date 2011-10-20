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
        - create baffles in original mesh with mappedWall patches
    - if extruding boundary faces:
        - convert boundary faces to mappedWall patches
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

    BBB=mapped between owner on original mesh and new extrusion.
        (zero offset)
    CCC=mapped between neighbour on original mesh and new extrusion
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

    BBB=mapped between original mesh and new extrusion
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
#include "fvMesh.H"
#include "polyTopoChange.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mappedWallPolyPatch.H"
#include "createShellMesh.H"
#include "syncTools.H"
#include "cyclicPolyPatch.H"
#include "wedgePolyPatch.H"
#include "nonuniformTransformCyclicPolyPatch.H"
#include "extrudeModel.H"
#include "globalIndex.H"
#include "addPatchCellLayer.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
//#include "ReadFields.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//template<class GeoField>
//void addPatchFields(const fvMesh& mesh, const word& patchFieldType)
//{
//    HashTable<const GeoField*> flds
//    (
//        mesh.objectRegistry::lookupClass<GeoField>()
//    );
//
//    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
//    {
//        const GeoField& fld = *iter();
//
//        typename GeoField::GeometricBoundaryField& bfld =
//            const_cast<typename GeoField::GeometricBoundaryField&>
//            (
//                fld.boundaryField()
//            );
//
//        label sz = bfld.size();
//
//        for (label i = 0; i < sz; i++)
//        {
//            bfld.set
//            (
//                i,
//                bfld.clone(GeoField::PatchFieldType::New
//                (
//                    patchFieldType,
//                    fld.mesh().boundary()[sz],
//                    fld.dimensionedInternalField()
//                )
//            );
//
//
//
//        Pout<< "fld:" << fld.name() << " had " << sz << " patches." << endl;
//        Pout<< "fld before:" << fld << endl;
//        Pout<< "adding on patch:" << fld.mesh().boundary()[sz].name() << endl;
//
//        bfld.setSize(sz+1);
//        bfld.set
//        (
//            sz,
//            GeoField::PatchFieldType::New
//            (
//                patchFieldType,
//                fld.mesh().boundary()[sz],
//                fld.dimensionedInternalField()
//            )
//        );
//
//        bfld[sz].operator=(pTraits<typename GeoField::value_type>::zero);
//
//        Pout<< "fld:" << fld.name() << " now " << bfld.size() << " patches."
//            << endl;
//
//        const typename GeoField::PatchFieldType& pfld = bfld[sz];
//        Pout<< "pfld:" << pfld << endl;
//
//
//        Pout<< "fld value:" << fld << endl;
//    }
//}


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


//void addCalculatedPatchFields(const fvMesh& mesh)
//{
//    addPatchFields<volScalarField>
//    (
//        mesh,
//        calculatedFvPatchField<scalar>::typeName
//    );
//    addPatchFields<volVectorField>
//    (
//        mesh,
//        calculatedFvPatchField<vector>::typeName
//    );
//    addPatchFields<volSphericalTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<sphericalTensor>::typeName
//    );
//    addPatchFields<volSymmTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<symmTensor>::typeName
//    );
//    addPatchFields<volTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<tensor>::typeName
//    );
//
//    // Surface fields
//
//    addPatchFields<surfaceScalarField>
//    (
//        mesh,
//        calculatedFvPatchField<scalar>::typeName
//    );
//    addPatchFields<surfaceVectorField>
//    (
//        mesh,
//        calculatedFvPatchField<vector>::typeName
//    );
//    addPatchFields<surfaceSphericalTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<sphericalTensor>::typeName
//    );
//    addPatchFields<surfaceSymmTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<symmTensor>::typeName
//    );
//    addPatchFields<surfaceTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<tensor>::typeName
//    );
//
//    // Point fields
//
//    addPatchFields<pointScalarField>
//    (
//        mesh,
//        calculatedFvPatchField<scalar>::typeName
//    );
//    addPatchFields<pointVectorField>
//    (
//        mesh,
//        calculatedFvPatchField<vector>::typeName
//    );
//}
//
//
//void addAllPatchFields(fvMesh& mesh, const label insertPatchI)
//{
//    polyBoundaryMesh& polyPatches =
//        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
//    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());
//
//    label sz = polyPatches.size();
//
//    addPatchFields<volScalarField>
//    (
//        mesh,
//        calculatedFvPatchField<scalar>::typeName
//    );
//    addPatchFields<volVectorField>
//    (
//        mesh,
//        calculatedFvPatchField<vector>::typeName
//    );
//    addPatchFields<volSphericalTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<sphericalTensor>::typeName
//    );
//    addPatchFields<volSymmTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<symmTensor>::typeName
//    );
//    addPatchFields<volTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<tensor>::typeName
//    );
//
//    // Surface fields
//
//    addPatchFields<surfaceScalarField>
//    (
//        mesh,
//        calculatedFvPatchField<scalar>::typeName
//    );
//    addPatchFields<surfaceVectorField>
//    (
//        mesh,
//        calculatedFvPatchField<vector>::typeName
//    );
//    addPatchFields<surfaceSphericalTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<sphericalTensor>::typeName
//    );
//    addPatchFields<surfaceSymmTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<symmTensor>::typeName
//    );
//    addPatchFields<surfaceTensorField>
//    (
//        mesh,
//        calculatedFvPatchField<tensor>::typeName
//    );
//
//    // Create reordering list
//    // patches before insert position stay as is
//    labelList oldToNew(sz);
//    for (label i = 0; i < insertPatchI; i++)
//    {
//        oldToNew[i] = i;
//    }
//    // patches after insert position move one up
//    for (label i = insertPatchI; i < sz-1; i++)
//    {
//        oldToNew[i] = i+1;
//    }
//    // appended patch gets moved to insert position
//    oldToNew[sz-1] = insertPatchI;
//
//    // Shuffle into place
//    polyPatches.reorder(oldToNew);
//    fvPatches.reorder(oldToNew);
//
//    reorderPatchFields<volScalarField>(mesh, oldToNew);
//    reorderPatchFields<volVectorField>(mesh, oldToNew);
//    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
//    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
//    reorderPatchFields<volTensorField>(mesh, oldToNew);
//    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
//    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
//    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
//    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
//    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);
//}


//// Adds patch if not yet there. Returns patchID.
//template<class PatchType>
//label addPatch(fvMesh& mesh, const word& patchName, const dictionary& dict)
//{
//    polyBoundaryMesh& polyPatches =
//        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
//
//    label patchI = polyPatches.findPatchID(patchName);
//    if (patchI != -1)
//    {
//        if (isA<PatchType>(polyPatches[patchI]))
//        {
//            // Already there
//            return patchI;
//        }
//        else
//        {
//            FatalErrorIn("addPatch<PatchType>(fvMesh&, const word&)")
//                << "Already have patch " << patchName
//                << " but of type " << PatchType::typeName
//                << exit(FatalError);
//        }
//    }
//
//
//    label insertPatchI = polyPatches.size();
//    label startFaceI = mesh.nFaces();
//
//    forAll(polyPatches, patchI)
//    {
//        const polyPatch& pp = polyPatches[patchI];
//
//        if (isA<processorPolyPatch>(pp))
//        {
//            insertPatchI = patchI;
//            startFaceI = pp.start();
//            break;
//        }
//    }
//
//    dictionary patchDict(dict);
//    patchDict.set("type", PatchType::typeName);
//    patchDict.set("nFaces", 0);
//    patchDict.set("startFace", startFaceI);
//
//
//    // Below is all quite a hack. Feel free to change once there is a better
//    // mechanism to insert and reorder patches.
//
//    // Clear local fields and e.g. polyMesh parallelInfo.
//    mesh.clearOut();
//
//    label sz = polyPatches.size();
//
//    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());
//
//    // Add polyPatch at the end
//    polyPatches.setSize(sz+1);
//    polyPatches.set
//    (
//        sz,
//        polyPatch::New
//        (
//            patchName,
//            patchDict,
//            insertPatchI,
//            polyPatches
//        )
//    );
//    fvPatches.setSize(sz+1);
//    fvPatches.set
//    (
//        sz,
//        fvPatch::New
//        (
//            polyPatches[sz],  // point to newly added polyPatch
//            mesh.boundary()
//        )
//    );
//
//    addAllPatchFields(mesh, insertPatchI);
//
//    return insertPatchI;
//}
//
//
//template<class PatchType>
//label addPatch(fvMesh& mesh, const word& patchName)
//{
//Pout<< "addPatch:" << patchName << endl;
//
//    polyBoundaryMesh& polyPatches =
//        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
//
//    label patchI = polyPatches.findPatchID(patchName);
//    if (patchI != -1)
//    {
//        if (isA<PatchType>(polyPatches[patchI]))
//        {
//            // Already there
//            return patchI;
//        }
//        else
//        {
//            FatalErrorIn("addPatch<PatchType>(fvMesh&, const word&)")
//                << "Already have patch " << patchName
//                << " but of type " << PatchType::typeName
//                << exit(FatalError);
//        }
//    }
//
//
//    label insertPatchI = polyPatches.size();
//    label startFaceI = mesh.nFaces();
//
//    forAll(polyPatches, patchI)
//    {
//        const polyPatch& pp = polyPatches[patchI];
//
//        if (isA<processorPolyPatch>(pp))
//        {
//            insertPatchI = patchI;
//            startFaceI = pp.start();
//            break;
//        }
//    }
//
//    // Below is all quite a hack. Feel free to change once there is a better
//    // mechanism to insert and reorder patches.
//
//    // Clear local fields and e.g. polyMesh parallelInfo.
//    mesh.clearOut();
//
//    label sz = polyPatches.size();
//
//    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());
//
//    // Add polyPatch at the end
//    polyPatches.setSize(sz+1);
//    polyPatches.set
//    (
//        sz,
//        polyPatch::New
//        (
//            PatchType::typeName,
//            patchName,
//            0,              // size
//            startFaceI,
//            insertPatchI,
//            polyPatches
//        )
//    );
//    fvPatches.setSize(sz+1);
//    fvPatches.set
//    (
//        sz,
//        fvPatch::New
//        (
//            polyPatches[sz],  // point to newly added polyPatch
//            mesh.boundary()
//        )
//    );
//
//    addAllPatchFields(mesh, insertPatchI);
//
//    return insertPatchI;
//}


label findPatchID(const List<polyPatch*>& newPatches, const word& name)
{
    forAll(newPatches, i)
    {
        if (newPatches[i]->name() == name)
        {
            return i;
        }
    }
    return -1;
}


template<class PatchType>
label addPatch
(
    const polyBoundaryMesh& patches,
    const word& patchName,
    DynamicList<polyPatch*>& newPatches
)
{
    label patchI = findPatchID(newPatches, patchName);

    if (patchI != -1)
    {
        if (isA<PatchType>(*newPatches[patchI]))
        {
            // Already there
            return patchI;
        }
        else
        {
            FatalErrorIn
            (
                "addPatch<PatchType>(const polyBoundaryMesh&,"
                " const word&, DynamicList<polyPatch*>)"
            )   << "Already have patch " << patchName
                << " but of type " << newPatches[patchI]->type()
                << exit(FatalError);
        }
    }


    patchI = newPatches.size();

    label startFaceI = 0;
    if (patchI > 0)
    {
        const polyPatch& pp = *newPatches.last();
        startFaceI = pp.start()+pp.size();
    }


    newPatches.append
    (
        polyPatch::New
        (
            PatchType::typeName,
            patchName,
            0,                          // size
            startFaceI,                 // nFaces
            patchI,
            patches
        ).ptr()
    );

    return patchI;
}


template<class PatchType>
label addPatch
(
    const polyBoundaryMesh& patches,
    const word& patchName,
    const dictionary& dict,
    DynamicList<polyPatch*>& newPatches
)
{
    label patchI = findPatchID(newPatches, patchName);

    if (patchI != -1)
    {
        if (isA<PatchType>(*newPatches[patchI]))
        {
            // Already there
            return patchI;
        }
        else
        {
            FatalErrorIn
            (
                "addPatch<PatchType>(const polyBoundaryMesh&,"
                " const word&, DynamicList<polyPatch*>)"
            )   << "Already have patch " << patchName
                << " but of type " << newPatches[patchI]->type()
                << exit(FatalError);
        }
    }


    patchI = newPatches.size();

    label startFaceI = 0;
    if (patchI > 0)
    {
        const polyPatch& pp = *newPatches.last();
        startFaceI = pp.start()+pp.size();
    }

    dictionary patchDict(dict);
    patchDict.set("type", PatchType::typeName);
    patchDict.set("nFaces", 0);
    patchDict.set("startFace", startFaceI);

    newPatches.append
    (
        polyPatch::New
        (
            patchName,
            patchDict,
            patchI,
            patches
        ).ptr()
    );

    return patchI;
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
    reorderPatchFields<pointScalarField>(mesh, oldToNew);
    reorderPatchFields<pointVectorField>(mesh, oldToNew);

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
    trimPatchFields<pointScalarField>(mesh, nNewPatches);
    trimPatchFields<pointVectorField>(mesh, nNewPatches);
}


// Remove zero-sized patches
void deleteEmptyPatches(fvMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList oldToNew(patches.size());
    label usedI = 0;
    label notUsedI = patches.size();

    //Pout<< "deleteEmptyPatches:" << endl;
    //forAll(patches, patchI)
    //{
    //    Pout<< "    patch:" << patchI << " name:" << patches[patchI].name()
    //        << " start:" << patches[patchI].start()
    //        << " nFaces:" << patches[patchI].size()
    //        << " index:" << patches[patchI].index()
    //        << endl;
    //}
    //Pout<< endl;


    // Add all the non-empty, non-processor patches
    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
            // Processor patches are unique per processor so look at local
            // size only
            if (patches[patchI].size() == 0)
            {
                Pout<< "Deleting processor patch " << patchI
                    << " name:" << patches[patchI].name()
                    << endl;
                oldToNew[patchI] = --notUsedI;
            }
            else
            {
                oldToNew[patchI] = usedI++;
            }
        }
        else
        {
            // All non-processor patches are present everywhere to reduce
            // size
            if (returnReduce(patches[patchI].size(), sumOp<label>()) == 0)
            {
                Pout<< "Deleting patch " << patchI
                    << " name:" << patches[patchI].name()
                    << endl;
                oldToNew[patchI] = --notUsedI;
            }
            else
            {
                oldToNew[patchI] = usedI++;
            }
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


// Count the number of faces in patches that need to be created. Calculates:
//  zoneSidePatch[zoneI]         : the number of side faces to be created
//  zoneZonePatch[zoneA,zoneB]   : the number of faces inbetween zoneA and B
// Since this only counts we're not taking the processor patches into
// account.
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


void addCouplingPatches
(
    const fvMesh& mesh,
    const word& regionName,
    const word& shellRegionName,
    const wordList& zoneNames,
    const wordList& zoneShadowNames,
    const boolList& isInternal,
    DynamicList<polyPatch*>& newPatches,
    labelList& interRegionTopPatch,
    labelList& interRegionBottomPatch
)
{
    Pout<< "Adding coupling patches:" << nl << nl
        << "patchID\tpatch\ttype" << nl
        << "-------\t-----\t----"
        << endl;

    interRegionTopPatch.setSize(zoneNames.size());
    interRegionBottomPatch.setSize(zoneNames.size());

    label nCoupled = 0;
    forAll(zoneNames, i)
    {
        word interName(regionName+"_to_"+shellRegionName+'_'+zoneNames[i]);

        if (isInternal[i])
        {
            interRegionTopPatch[i] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                interName + "_top",
                newPatches
            );
            nCoupled++;
            Pout<< interRegionTopPatch[i]
                << '\t' << newPatches[interRegionTopPatch[i]]->name()
                << '\t' << newPatches[interRegionTopPatch[i]]->type()
                << nl;

            interRegionBottomPatch[i] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                interName + "_bottom",
                newPatches
            );
            nCoupled++;
            Pout<< interRegionBottomPatch[i]
                << '\t' << newPatches[interRegionBottomPatch[i]]->name()
                << '\t' << newPatches[interRegionBottomPatch[i]]->type()
                << nl;
        }
        else if (zoneShadowNames.size() == 0)
        {
            interRegionTopPatch[i] = addPatch<polyPatch>
            (
                mesh.boundaryMesh(),
                zoneNames[i] + "_top",
                newPatches
            );
            nCoupled++;
            Pout<< interRegionTopPatch[i]
                << '\t' << newPatches[interRegionTopPatch[i]]->name()
                << '\t' << newPatches[interRegionTopPatch[i]]->type()
                << nl;

            interRegionBottomPatch[i] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                interName,
                newPatches
            );
            nCoupled++;
            Pout<< interRegionBottomPatch[i]
                << '\t' << newPatches[interRegionBottomPatch[i]]->name()
                << '\t' << newPatches[interRegionBottomPatch[i]]->type()
                << nl;
        }
        else    //patch using shadow face zones.
        {
            interRegionTopPatch[i] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                zoneShadowNames[i] + "_top",
                newPatches
            );
            nCoupled++;
            Pout<< interRegionTopPatch[i]
                << '\t' << newPatches[interRegionTopPatch[i]]->name()
                << '\t' << newPatches[interRegionTopPatch[i]]->type()
                << nl;

            interRegionBottomPatch[i] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                interName,
                newPatches
            );
            nCoupled++;
            Pout<< interRegionBottomPatch[i]
                << '\t' << newPatches[interRegionBottomPatch[i]]->name()
                << '\t' << newPatches[interRegionBottomPatch[i]]->type()
                << nl;
        }
    }
    Pout<< "Added " << nCoupled << " inter-region patches." << nl
        << endl;
}


void addProcPatches
(
    const fvMesh& mesh,
    const primitiveFacePatch& extrudePatch,
    const labelList& extrudeMeshFaces,

    labelList& sidePatchID,
    DynamicList<polyPatch*>& newPatches
)
{
    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    indirectPrimitivePatch pp
    (
        IndirectList<face>(mesh.faces(), extrudeMeshFaces),
        mesh.points()
    );

    forAll(extrudePatch.edges(), edgeI)
    {
        const edge& extrudeEdge = extrudePatch.edges()[edgeI];
        const edge& ppEdge = pp.edges()[edgeI];

        if (extrudeEdge != ppEdge)
        {
            FatalErrorIn("addProcPatches()")
                << "Problem: extrudeEdge:" << extrudeEdge
                << " ppEdge:" << ppEdge
                << exit(FatalError);
        }
    }


    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches).
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            pp
        )
    );

    // Calculate for every edge the patch to use. This will be an existing
    // patch for uncoupled edge and a possibly new patch (in patchToNbrProc)
    // for processor patches.
    label nNewPatches;
    Map<label> nbrProcToPatch;
    Map<label> patchToNbrProc;

    addPatchCellLayer::calcSidePatch
    (
        mesh,
        globalFaces,
        edgeGlobalFaces,
        pp,

        sidePatchID,
        nNewPatches,
        nbrProcToPatch,
        patchToNbrProc
    );


    // All patchIDs calcSidePatch are in mesh.boundaryMesh() numbering.
    // Redo in newPatches numbering.

    Pout<< "Adding inter-processor patches:" << nl << nl
        << "patchID\tpatch" << nl
        << "-------\t-----"
        << endl;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(sidePatchID, edgeI)
    {
        label meshPatchI = sidePatchID[edgeI];
        if (meshPatchI != -1)
        {
            label newPatchI = -1;
            if (meshPatchI < patches.size())
            {
                // Existing mesh patch. See if I already have it.
                newPatchI = findPatchID
                (
                    newPatches,
                    patches[meshPatchI].name()
                );
            }

            if (newPatchI == -1)
            {
                // Newly added processor patch
                label nbrProcI = patchToNbrProc[meshPatchI];
                word name =
                        "procBoundary"
                      + Foam::name(Pstream::myProcNo())
                      + "to"
                      + Foam::name(nbrProcI);

                dictionary patchDict;
                patchDict.add("myProcNo", Pstream::myProcNo());
                patchDict.add("neighbProcNo", nbrProcI);

                newPatchI = addPatch<processorPolyPatch>
                (
                    mesh.boundaryMesh(),
                    name,
                    patchDict,
                    newPatches
                );

                Pout<< newPatchI << '\t' << name
                    << nl;
            }

            sidePatchID[edgeI] = newPatchI;
        }
    }


    // Clear out sidePatchID for uncoupled edges. Just so we don't have
    // to expose all the globalEdgeFaces info.
    forAll(sidePatchID, edgeI)
    {
        if
        (
            edgeGlobalFaces[edgeI].size() == 2
         && pp.edgeFaces()[edgeI].size() == 1
        )
        {}
        else
        {
            sidePatchID[edgeI] = -1;
        }
    }
}


void addZoneSidePatches
(
    const fvMesh& mesh,
    const word& oneDPolyPatchType,
    const wordList& zoneNames,

    DynamicList<polyPatch*>& newPatches,
    labelList& zoneSidePatch
)
{
    Pout<< "Adding patches for sides on zones:" << nl << nl
        << "patchID\tpatch" << nl
        << "-------\t-----"
        << endl;

    label nSide = 0;

    forAll(zoneNames, zoneI)
    {
        if (oneDPolyPatchType != word::null)
        {
            // Reuse single empty patch.
            word patchName;
            if (oneDPolyPatchType == "emptyPolyPatch")
            {
                patchName = "oneDEmptyPatch";
                zoneSidePatch[zoneI] = addPatch<emptyPolyPatch>
                (
                    mesh.boundaryMesh(),
                    patchName,
                    newPatches
                );
            }
            else if (oneDPolyPatchType == "wedgePolyPatch")
            {
                patchName = "oneDWedgePatch";
                zoneSidePatch[zoneI] = addPatch<wedgePolyPatch>
                (
                    mesh.boundaryMesh(),
                    patchName,
                    newPatches
                );
            }
            else
            {
                FatalErrorIn("addZoneSidePatches()")
                    << "Type " << oneDPolyPatchType << " does not exist "
                    << exit(FatalError);
            }

            Pout<< zoneSidePatch[zoneI] << '\t' << patchName << nl;

            nSide++;
        }
        else if (zoneSidePatch[zoneI] > 0)
        {
            word patchName = zoneNames[zoneI] + "_" + "side";

            zoneSidePatch[zoneI] = addPatch<polyPatch>
            (
                mesh.boundaryMesh(),
                patchName,
                newPatches
            );

            Pout<< zoneSidePatch[zoneI] << '\t' << patchName << nl;

            nSide++;
        }
    }
    Pout<< "Added " << nSide << " zone-side patches." << nl
        << endl;
}


void addInterZonePatches
(
    const fvMesh& mesh,
    const wordList& zoneNames,
    const bool oneD,

    labelList& zoneZonePatch_min,
    labelList& zoneZonePatch_max,
    DynamicList<polyPatch*>& newPatches
)
{
    Pout<< "Adding inter-zone patches:" << nl << nl
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
            for (label maxZone = minZone; maxZone < zoneNames.size(); maxZone++)
            {
                label index = minZone*zoneNames.size()+maxZone;

                if (zoneZonePatch_min[index] > 0)
                {
                    word minToMax =
                        zoneNames[minZone]
                      + "_to_"
                      + zoneNames[maxZone];
                    word maxToMin =
                        zoneNames[maxZone]
                      + "_to_"
                      + zoneNames[minZone];

                    {
                        transformDict.set("neighbourPatch", maxToMin);
                        zoneZonePatch_min[index] =
                        addPatch<nonuniformTransformCyclicPolyPatch>
                        (
                            mesh.boundaryMesh(),
                            minToMax,
                            transformDict,
                            newPatches
                        );
                        Pout<< zoneZonePatch_min[index] << '\t' << minToMax
                            << nl;
                        nInter++;
                    }
                    {
                        transformDict.set("neighbourPatch", minToMax);
                        zoneZonePatch_max[index] =
                        addPatch<nonuniformTransformCyclicPolyPatch>
                        (
                            mesh.boundaryMesh(),
                            maxToMin,
                            transformDict,
                            newPatches
                        );
                        Pout<< zoneZonePatch_max[index] << '\t' << maxToMin
                            << nl;
                        nInter++;
                    }

                }
            }
        }
    }
    Pout<< "Added " << nInter << " inter-zone patches." << nl
        << endl;
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


void setCouplingInfo
(
    fvMesh& mesh,
    const labelList& zoneToPatch,
    const word& sampleRegion,
    const mappedWallPolyPatch::sampleMode mode,
    const List<pointField>& offsets
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    List<polyPatch*> newPatches(patches.size(), NULL);

    forAll(zoneToPatch, zoneI)
    {
        label patchI = zoneToPatch[zoneI];

        const polyPatch& pp = patches[patchI];

        if (isA<mappedWallPolyPatch>(pp))
        {
            newPatches[patchI] = new mappedWallPolyPatch
            (
                pp.name(),
                pp.size(),
                pp.start(),
                patchI,
                sampleRegion,                           // sampleRegion
                mode,                                   // sampleMode
                pp.name(),                              // samplePatch
                offsets[zoneI],                         // offset
                patches
            );
        }
    }

    forAll(newPatches, patchI)
    {
        if (!newPatches[patchI])
        {
            newPatches[patchI] = patches[patchI].clone(patches).ptr();
        }
    }

    mesh.removeFvBoundary();
    mesh.addFvPatches(newPatches, true);
}


// Main program:

int main(int argc, char *argv[])
{
    argList::addNote("Create region mesh by extruding a faceZone");

    #include "addRegionOption.H"
    #include "addOverwriteOption.H"
    argList::addNote("Create region mesh by extruding a faceZone");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    if (mesh.boundaryMesh().checkParallelSync(true))
    {
        List<wordList> allNames(Pstream::nProcs());
        allNames[Pstream::myProcNo()] = mesh.boundaryMesh().names();
        Pstream::gatherList(allNames);
        Pstream::scatterList(allNames);

        FatalErrorIn(args.executable())
            << "Patches are not synchronised on all processors."
            << " Per processor patches " << allNames
            << exit(FatalError);
    }


    const word oldInstance = mesh.pointsInstance();
    bool overwrite = args.optionFound("overwrite");
    const word dictName
        (args.optionLookupOrDefault<word>("dict", "extrudeToRegionMeshDict"));

    IOdictionary dict
    (
        IOobject
        (
            dictName,
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

    mappedPatchBase::sampleMode sampleMode =
        mappedPatchBase::sampleModeNames_[dict.lookup("sampleMode")];

    const Switch oneD(dict.lookup("oneD"));
    const Switch adaptMesh(dict.lookup("adaptMesh"));

    Pout<< "Extruding zones " << zoneNames
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



    //// Read objects in time directory
    //IOobjectList objects(mesh, runTime.timeName());
    //
    //// Read vol fields.
    //
    //PtrList<volScalarField> vsFlds;
    //ReadFields(mesh, objects, vsFlds);
    //
    //PtrList<volVectorField> vvFlds;
    //ReadFields(mesh, objects, vvFlds);
    //
    //PtrList<volSphericalTensorField> vstFlds;
    //ReadFields(mesh, objects, vstFlds);
    //
    //PtrList<volSymmTensorField> vsymtFlds;
    //ReadFields(mesh, objects, vsymtFlds);
    //
    //PtrList<volTensorField> vtFlds;
    //ReadFields(mesh, objects, vtFlds);
    //
    //// Read surface fields.
    //
    //PtrList<surfaceScalarField> ssFlds;
    //ReadFields(mesh, objects, ssFlds);
    //
    //PtrList<surfaceVectorField> svFlds;
    //ReadFields(mesh, objects, svFlds);
    //
    //PtrList<surfaceSphericalTensorField> sstFlds;
    //ReadFields(mesh, objects, sstFlds);
    //
    //PtrList<surfaceSymmTensorField> ssymtFlds;
    //ReadFields(mesh, objects, ssymtFlds);
    //
    //PtrList<surfaceTensorField> stFlds;
    //ReadFields(mesh, objects, stFlds);
    //
    //// Read point fields.
    //
    //PtrList<pointScalarField> psFlds;
    //ReadFields(pointMesh::New(mesh), objects, psFlds);
    //
    //PtrList<pointVectorField> pvFlds;
    //ReadFields(pointMesh::New(mesh), objects, pvFlds);



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
    Pout<< "Writing meshes to " << meshInstance << nl << endl;


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


    Pout<< "extrudePatch :"
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
    // since then how to couple? - mapped only valid for a whole patch.
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
            Pout<< "FaceZone " << fz.name() << " has internal faces" << endl;
        }
        else
        {
            Pout<< "FaceZone " << fz.name() << " has boundary faces" << endl;
        }
    }
    Pout<< endl;



    DynamicList<polyPatch*> regionPatches(patches.size());
    // Copy all non-local patches since these are used on boundary edges of
    // the extrusion
    forAll(patches, patchI)
    {
        if (!isA<processorPolyPatch>(patches[patchI]))
        {
            label newPatchI = regionPatches.size();
            regionPatches.append
            (
                patches[patchI].clone
                (
                    patches,
                    newPatchI,
                    0,              // size
                    0               // start
                ).ptr()
            );
        }
    }


    // Add interface patches
    // ~~~~~~~~~~~~~~~~~~~~~

    // From zone to interface patch (region side)
    labelList interRegionTopPatch(zoneNames.size());
    labelList interRegionBottomPatch(zoneNames.size());

    addCouplingPatches
    (
        mesh,
        regionName,
        shellRegionName,
        zoneNames,
        zoneShadowNames,
        isInternal,

        regionPatches,
        interRegionTopPatch,
        interRegionBottomPatch
    );

    // From zone to interface patch (mesh side)
    labelList interMeshTopPatch(zoneNames.size());
    labelList interMeshBottomPatch(zoneNames.size());

    if (adaptMesh)
    {
        DynamicList<polyPatch*> newPatches(patches.size());
        forAll(patches, patchI)
        {
            newPatches.append(patches[patchI].clone(patches).ptr());
        }

        addCouplingPatches
        (
            mesh,
            regionName,
            shellRegionName,
            zoneNames,
            zoneShadowNames,
            isInternal,

            newPatches,
            interMeshTopPatch,
            interMeshBottomPatch
        );
        mesh.clearOut();
        mesh.removeFvBoundary();
        mesh.addFvPatches(newPatches, true);
    }


    // Patch per extruded face
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

    // Now we'll have:
    //  zoneSidePatch[zoneA] : number of faces needed on the side of zoneA
    //  zoneZonePatch_min[zoneA,zoneB] : number of faces needed inbetween A,B


    // Add the zone-side patches.
    addZoneSidePatches
    (
        mesh,
        (oneD ? dict.lookup("oneDPolyPatchType") : word::null),
        zoneNames,

        regionPatches,
        zoneSidePatch
    );


    // Add the patches inbetween zones
    addInterZonePatches
    (
        mesh,
        zoneNames,
        oneD,

        zoneZonePatch_min,
        zoneZonePatch_max,
        regionPatches
    );


    // Additionally check if any interprocessor patches need to be added.
    // Reuses addPatchCellLayer functionality.
    // Note: does not handle edges with > 2 faces?
    labelList sidePatchID;
    addProcPatches
    (
        mesh,
        extrudePatch,
        extrudeMeshFaces,

        sidePatchID,
        regionPatches
    );


//    // Add all the newPatches to the mesh and fields
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    {
//        forAll(newPatches, patchI)
//        {
//            Pout<< "Adding patch " << patchI
//                << " name:" << newPatches[patchI]->name()
//                << endl;
//        }
//        //label nOldPatches = mesh.boundary().size();
//        mesh.clearOut();
//        mesh.removeFvBoundary();
//        mesh.addFvPatches(newPatches, true);
//        //// Add calculated fvPatchFields for the added patches
//        //for
//        //(
//        //    label patchI = nOldPatches;
//        //    patchI < mesh.boundary().size();
//        //    patchI++
//        //)
//        //{
//        //    Pout<< "ADDing calculated to patch " << patchI
//        //        << endl;
//        //    addCalculatedPatchFields(mesh);
//        //}
//        //Pout<< "** Added " << mesh.boundary().size()-nOldPatches
//        //    << " patches." << endl;
//    }


    // Set patches to use for edges to be extruded into boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // In order of edgeFaces: per edge, per originating face the
    // patch to use for the side face (from the extruded edge).
    // If empty size create an internal face.
    labelListList extrudeEdgePatches(extrudePatch.nEdges());

    // Is edge a non-manifold edge
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
        else if (sidePatchID[edgeI] != -1)
        {
            // Coupled extrusion
            ePatches.setSize(eFaces.size());
            forAll(eFaces, i)
            {
                ePatches[i] = sidePatchID[edgeI];
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
                label newPatchI = findPatchID
                (
                    regionPatches,
                    patches[patches.whichPatch(faceI)].name()
                );
                ePatches.setSize(eFaces.size(), newPatchI);
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
    faceList pointGlobalRegions;
    faceList pointLocalRegions;
    labelList localToGlobalRegion;

    createShellMesh::calcPointRegions
    (
        mesh.globalData(),
        extrudePatch,
        nonManifoldEdge,

        pointGlobalRegions,
        pointLocalRegions,
        localToGlobalRegion
    );

    // Per local region an originating point
    labelList localRegionPoints(localToGlobalRegion.size());
    forAll(pointLocalRegions, faceI)
    {
        const face& f = extrudePatch.localFaces()[faceI];
        const face& pRegions = pointLocalRegions[faceI];
        forAll(pRegions, fp)
        {
            localRegionPoints[pRegions[fp]] = f[fp];
        }
    }

    // Calculate region normals by reducing local region normals
    pointField localRegionNormals(localToGlobalRegion.size());
    {
        pointField localSum(localToGlobalRegion.size(), vector::zero);
        labelList localNum(localToGlobalRegion.size(), 0);

        forAll(pointLocalRegions, faceI)
        {
            const face& pRegions = pointLocalRegions[faceI];
            forAll(pRegions, fp)
            {
                label localRegionI = pRegions[fp];
                localSum[localRegionI] += extrudePatch.faceNormals()[faceI];
                localNum[localRegionI]++;
            }
        }

        Map<point> globalSum(2*localToGlobalRegion.size());
        Map<label> globalNum(2*localToGlobalRegion.size());

        forAll(localSum, localRegionI)
        {
            label globalRegionI = localToGlobalRegion[localRegionI];
            globalSum.insert(globalRegionI, localSum[localRegionI]);
            globalNum.insert(globalRegionI, localNum[localRegionI]);
        }

        // Reduce
        Pstream::mapCombineGather(globalSum, plusEqOp<point>());
        Pstream::mapCombineScatter(globalSum);
        Pstream::mapCombineGather(globalNum, plusEqOp<label>());
        Pstream::mapCombineScatter(globalNum);

        forAll(localToGlobalRegion, localRegionI)
        {
            label globalRegionI = localToGlobalRegion[localRegionI];
            localRegionNormals[localRegionI] =
                globalSum[globalRegionI]
              / globalNum[globalRegionI];
        }
        localRegionNormals /= mag(localRegionNormals);
    }



    // For debugging: dump hedgehog plot of normals
    if (false)
    {
        OFstream str(runTime.path()/"localRegionNormals.obj");
        label vertI = 0;

        scalar thickness = model().sumThickness(1);

        forAll(pointLocalRegions, faceI)
        {
            const face& f = extrudeFaces[faceI];

            forAll(f, fp)
            {
                label region = pointLocalRegions[faceI][fp];
                const point& pt = extrudePoints[f[fp]];

                meshTools::writeOBJ(str, pt);
                vertI++;
                meshTools::writeOBJ
                (
                    str,
                    pt+thickness*localRegionNormals[region]
                );
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Use model to create displacements of first layer
    vectorField firstDisp(localRegionNormals.size());
    forAll(firstDisp, regionI)
    {
        //const point& regionPt = regionCentres[regionI];
        const point& regionPt = extrudePatch.points()
        [
            extrudePatch.meshPoints()
            [
                localRegionPoints[regionI]
            ]
        ];
        const vector& n = localRegionNormals[regionI];
        firstDisp[regionI] = model()(regionPt, n, 1) - regionPt;
    }



    // Create a new mesh
    // ~~~~~~~~~~~~~~~~~

    createShellMesh extruder
    (
        extrudePatch,
        pointLocalRegions,
        localRegionPoints
    );


    autoPtr<mapPolyMesh> shellMap;
    fvMesh regionMesh
    (
        IOobject
        (
            shellRegionName,
            meshInstance,
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        xferCopy(pointField()),
        xferCopy(faceList()),
        xferCopy(labelList()),
        xferCopy(labelList()),
        false
    );

    // Add the new patches
    forAll(regionPatches, patchI)
    {
        regionPatches[patchI] = regionPatches[patchI]->clone
        (
            regionMesh.boundaryMesh()
        ).ptr();
    }
    regionMesh.clearOut();
    regionMesh.removeFvBoundary();
    regionMesh.addFvPatches(regionPatches, true);

    {
        polyTopoChange meshMod(regionPatches.size());

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

        shellMap = meshMod.changeMesh
        (
            regionMesh,     // mesh to change
            false           // inflate
        );
    }

    // Necessary?
    regionMesh.setInstance(meshInstance);


    // Update numbering on extruder.
    extruder.updateMesh(shellMap);


    // Calculate offsets from shell mesh back to original mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<pointField> topOffsets(zoneIDs.size());
    List<pointField> bottomOffsets(zoneIDs.size());

    forAll(regionMesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = regionMesh.boundaryMesh()[patchI];

        if
        (
            isA<mappedWallPolyPatch>(pp)
         && (findIndex(interRegionTopPatch, patchI) != -1)
        )
        {
            label zoneI = findIndex(interRegionTopPatch, patchI);
            topOffsets[zoneI] = calcOffset(extrudePatch, extruder, pp);
        }
        else if
        (
            isA<mappedWallPolyPatch>(pp)
         && (findIndex(interRegionBottomPatch, patchI) != -1)
        )
        {
            label zoneI = findIndex(interRegionBottomPatch, patchI);

            bottomOffsets[zoneI] = calcOffset(extrudePatch, extruder, pp);
        }
    }


    // Change top and bottom boundary conditions on regionMesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Correct top patches for offset
        setCouplingInfo
        (
            regionMesh,
            interRegionTopPatch,
            regionName,                 // name of main mesh
            sampleMode,                 // sampleMode
            topOffsets
        );

        // Correct bottom patches for offset
        setCouplingInfo
        (
            regionMesh,
            interRegionBottomPatch,
            regionName,
            sampleMode,                 // sampleMode
            bottomOffsets
        );

        // Remove any unused patches
        deleteEmptyPatches(regionMesh);
    }

    // Change top and bottom boundary conditions on main mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (adaptMesh)
    {
        // Correct top patches for offset
        setCouplingInfo
        (
            mesh,
            interMeshTopPatch,
            shellRegionName,                        // name of shell mesh
            sampleMode,                             // sampleMode
            -topOffsets
        );

        // Correct bottom patches for offset
        setCouplingInfo
        (
            mesh,
            interMeshBottomPatch,
            shellRegionName,
            sampleMode,
            -bottomOffsets
        );
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


    Pout<< "Writing mesh " << regionMesh.name()
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
                    interMeshBottomPatch[zoneI],// patch for face
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
                    interMeshBottomPatch[zoneI],    // patch for face
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
                        interMeshTopPatch[zoneI],   // patch for face
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
                        interMeshTopPatch[zoneI],       // patch for face
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
                label zoneI = zoneID[zoneFaceI];
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
                            interMeshTopPatch[zoneI],       // patch for face
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
                        interMeshTopPatch[zoneI],       // patch for face
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


        Pout<< "Writing mesh " << mesh.name()
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
