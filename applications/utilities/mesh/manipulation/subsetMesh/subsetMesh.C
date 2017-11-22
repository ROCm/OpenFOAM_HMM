/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    subsetMesh

Group
    grpMeshManipulationUtilities

Description
    Selects a section of mesh based on a cellSet.

    The utility sub-sets the mesh to choose only a part of interest. Check
    the setSet/cellSet/topoSet utilities to see how to select cells based on
    various shapes.

    The mesh will subset all points, faces and cells needed to make a sub-mesh
    but will not preserve attached boundary types.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "argList.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "topoDistanceData.H"
#include "FaceCellWave.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "ReadFields.H"
#include "processorMeshes.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

labelList nearestPatch(const polyMesh& mesh, const labelList& patchIDs)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Count number of faces in exposedPatchIDs
    label nFaces = 0;
    forAll(patchIDs, i)
    {
        const polyPatch& pp = pbm[patchIDs[i]];
        nFaces += pp.size();
    }

    // Field on cells and faces.
    List<topoDistanceData> cellData(mesh.nCells());
    List<topoDistanceData> faceData(mesh.nFaces());

    // Start of changes
    labelList patchFaces(nFaces);
    List<topoDistanceData> patchData(nFaces);
    nFaces = 0;
    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];
        const polyPatch& pp = pbm[patchI];

        forAll(pp, i)
        {
            patchFaces[nFaces] = pp.start()+i;
            patchData[nFaces] = topoDistanceData(patchI, 0);
            nFaces++;
        }
    }

    // Propagate information inwards
    FaceCellWave<topoDistanceData> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        faceData,
        cellData,
        mesh.globalData().nTotalCells()+1
    );

    // And extract

    labelList nearest(mesh.nFaces());

    bool haveWarned = false;
    forAll(faceData, faceI)
    {
        if (!faceData[faceI].valid(deltaCalc.data()))
        {
            if (!haveWarned)
            {
                WarningInFunction
                    << "Did not visit some faces, e.g. face " << faceI
                    << " at " << mesh.faceCentres()[faceI] << nl
                    << "Using patch " << patchIDs[0] << " as nearest"
                    << endl;
                haveWarned = true;
            }
            nearest[faceI] = patchIDs[0];
        }
        else
        {
            nearest[faceI] = faceData[faceI].data();
        }
    }

    return nearest;
}


//
// Subset field-type, availability information cached
// in the availableFields hashtable.
//
template<class Type, template<class> class PatchField, class GeoMesh>
void subsetFields
(
    const fvMeshSubset& subsetter,
    HashTable<wordHashSet>& availableFields,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& subFields
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> FieldType;
    const word fieldType = FieldType::typeName;

    const wordList fieldNames = availableFields(fieldType).sortedToc();
    subFields.setSize(fieldNames.size());

    const fvMesh& baseMesh = subsetter.baseMesh();

    if (fieldNames.empty())
    {
        return;
    }
    Info<< "Subsetting " << fieldType << " (";
    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];
        if (i) Info<< ' ';
        Info<< fieldName;

        FieldType fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, subsetter.interpolate(fld));

        // Subsetting adds 'subset' prefix - rename to match original.
        subFields[i].rename(fieldName);
    }
    Info<< ")" << nl;
}


template<class Type>
void subsetPointFields
(
    const fvMeshSubset& subsetter,
    const pointMesh& pMesh,
    HashTable<wordHashSet>& availableFields,
    PtrList<GeometricField<Type, pointPatchField, pointMesh>>& subFields
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> FieldType;
    const word fieldType = FieldType::typeName;

    const wordList fieldNames = availableFields(fieldType).sortedToc();
    subFields.setSize(fieldNames.size());

    const fvMesh& baseMesh = subsetter.baseMesh();

    if (fieldNames.empty())
    {
        return;
    }
    Info<< "Subsetting " << fieldType << " (";
    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];
        if (i) Info<< ' ';
        Info<< fieldName;

        FieldType fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            pMesh
        );

        subFields.set(i, subsetter.interpolate(fld));

        // Subsetting adds 'subset' prefix - rename to match original.
        subFields[i].rename(fieldName);
    }
    Info<< ")" << nl;
}


template<class Type>
void subsetDimensionedFields
(
    const fvMeshSubset& subsetter,
    HashTable<wordHashSet>& availableFields,
    PtrList<DimensionedField<Type, volMesh>>& subFields
)
{
    typedef DimensionedField<Type, volMesh> FieldType;
    const word fieldType = FieldType::typeName;

    const wordList fieldNames = availableFields(fieldType).sortedToc();
    subFields.setSize(fieldNames.size());

    const fvMesh& baseMesh = subsetter.baseMesh();

    if (fieldNames.empty())
    {
        return;
    }
    Info<< "Subsetting " << fieldType << " (";
    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];
        if (i) Info<< ' ';
        Info<< fieldName;

        FieldType fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, subsetter.interpolate(fld));

        // Subsetting adds 'subset' prefix - rename to match original.
        subFields[i].rename(fieldName);
    }
    Info<< ")" << nl;
}


template<class TopoSet>
void subsetTopoSets
(
    const fvMesh& mesh,
    const IOobjectList& objectsList,
    const labelList& map,
    const fvMesh& subMesh,
    PtrList<TopoSet>& subSets
)
{
    // Read original sets
    PtrList<TopoSet> sets;
    ReadFields<TopoSet>(objectsList, sets);

    subSets.setSize(sets.size());
    forAll(sets, i)
    {
        TopoSet& set = sets[i];

        Info<< "Subsetting " << set.type() << " " << set.name() << endl;

        // Map the data
        PackedBoolList isSet(set.maxSize(mesh));
        forAllConstIters(set, iter)
        {
            isSet[iter.key()] = true;
        }
        label nSet = 0;
        forAll(map, i)
        {
            if (isSet[map[i]])
            {
                nSet++;
            }
        }

        subSets.set
        (
            i,
            new TopoSet(subMesh, set.name(), nSet, IOobject::AUTO_WRITE)
        );
        TopoSet& subSet = subSets[i];
        forAll(map, i)
        {
            if (isSet[map[i]])
            {
                subSet.insert(i);
            }
        }
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "select a mesh subset based on a cellSet"
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    argList::addArgument("cellSet");
    argList::addOption
    (
        "patch",
        "name",
        "add exposed internal faces to specified patch instead of to "
        "'oldInternalFaces'"
    );
    argList::addOption
    (
        "patches",
        "names",
        "add exposed internal faces to nearest of specified patches"
        " instead of to 'oldInternalFaces'"
    );
    argList::addOption
    (
        "resultTime",
        "time",
        "specify a time for the resulting mesh"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    #include "createNamedMesh.H"

    const word setName = args[1];

    word meshInstance = mesh.pointsInstance();
    word fieldsInstance = runTime.timeName();

    const bool overwrite = args.optionFound("overwrite");
    const bool specifiedInstance = args.optionReadIfPresent
    (
        "resultTime",
        fieldsInstance
    );
    if (specifiedInstance)
    {
        // Set both mesh and field to this time
        meshInstance = fieldsInstance;
    }


    Info<< "Reading cell set from " << setName << nl << endl;

    // Create mesh subsetting engine
    fvMeshSubset subsetter(mesh);

    labelList exposedPatchIDs;

    if (args.optionFound("patch"))
    {
        const word patchName = args["patch"];

        exposedPatchIDs = { mesh.boundaryMesh().findPatchID(patchName) };

        if (exposedPatchIDs[0] == -1)
        {
            FatalErrorInFunction
                << nl << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        Info<< "Adding exposed internal faces to patch " << patchName
            << nl << endl;
    }
    else if (args.optionFound("patches"))
    {
        const wordReList patchNames(args.optionRead<wordReList>("patches"));

        exposedPatchIDs = mesh.boundaryMesh().patchSet(patchNames).sortedToc();

        Info<< "Adding exposed internal faces to nearest of patches "
            << patchNames << nl << endl;
    }
    else
    {
        Info<< "Adding exposed internal faces to a patch called"
            << " \"oldInternalFaces\" (created if necessary)" << endl
            << endl;
        exposedPatchIDs = { -1 };
    }


    cellSet currentSet(mesh, setName);

    if (exposedPatchIDs.size() == 1)
    {
        subsetter.setLargeCellSubset(currentSet, exposedPatchIDs[0], true);
    }
    else
    {
        // Find per face the nearest patch
        labelList nearestExposedPatch(nearestPatch(mesh, exposedPatchIDs));

        labelList region(mesh.nCells(), 0);
        forAllConstIter(cellSet, currentSet, iter)
        {
            region[iter.key()] = 1;
        }

        labelList exposedFaces(subsetter.getExposedFaces(region, 1, true));
        subsetter.setLargeCellSubset
        (
            region,
            1,
            exposedFaces,
            labelUIndList(nearestExposedPatch, exposedFaces)(),
            true
        );
    }


    IOobjectList objects(mesh, runTime.timeName());
    HashTable<wordHashSet> availableFields = objects.classes();


    // Read vol fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<volScalarField> scalarFlds;
    subsetFields(subsetter, availableFields, scalarFlds);

    PtrList<volVectorField> vectorFlds;
    subsetFields(subsetter, availableFields, vectorFlds);

    PtrList<volSphericalTensorField> sphTensorFlds;
    subsetFields(subsetter, availableFields, sphTensorFlds);

    PtrList<volSymmTensorField> symmTensorFlds;
    subsetFields(subsetter, availableFields, symmTensorFlds);

    PtrList<volTensorField> tensorFlds;
    subsetFields(subsetter, availableFields, tensorFlds);


    // Read surface fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<surfaceScalarField> surfScalarFlds;
    subsetFields(subsetter, availableFields, surfScalarFlds);

    PtrList<surfaceVectorField> surfVectorFlds;
    subsetFields(subsetter, availableFields, surfVectorFlds);

    PtrList<surfaceSphericalTensorField> surfSphTensorFlds;
    subsetFields(subsetter, availableFields, surfSphTensorFlds);

    PtrList<surfaceSymmTensorField> surfSymmTensorFlds;
    subsetFields(subsetter, availableFields, surfSymmTensorFlds);

    PtrList<surfaceTensorField> surfTensorFlds;
    subsetFields(subsetter, availableFields, surfTensorFlds);


    // Read point fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const pointMesh& pMesh = pointMesh::New(mesh);

    PtrList<pointScalarField> pointScalarFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pointScalarFlds);

    PtrList<pointVectorField> pointVectorFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pointVectorFlds);

    PtrList<pointSphericalTensorField> pointSphTensorFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pointSphTensorFlds);

    PtrList<pointSymmTensorField> pointSymmTensorFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pointSymmTensorFlds);

    PtrList<pointTensorField> pointTensorFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pointTensorFlds);


    // Read dimensioned fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<volScalarField::Internal> scalarDimFlds;
    subsetDimensionedFields(subsetter, availableFields, scalarDimFlds);

    PtrList<volVectorField::Internal> vectorDimFlds;
    subsetDimensionedFields(subsetter, availableFields, vectorDimFlds);

    PtrList<volSphericalTensorField::Internal> sphTensorDimFlds;
    subsetDimensionedFields(subsetter, availableFields, sphTensorDimFlds);

    PtrList<volSymmTensorField::Internal> symmTensorDimFlds;
    subsetDimensionedFields(subsetter, availableFields, symmTensorDimFlds);

    PtrList<volTensorField::Internal> tensorDimFlds;
    subsetDimensionedFields(subsetter, availableFields, tensorDimFlds);


    // topoSets and subset
    // ~~~~~~~~~~~~~~~~~~~

    PtrList<cellSet> cellSets;
    PtrList<faceSet> faceSets;
    PtrList<pointSet> pointSets;

    {
        IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
        objects.remove(currentSet);
        subsetTopoSets
        (
            mesh,
            objects,
            subsetter.cellMap(),
            subsetter.subMesh(),
            cellSets
        );
        subsetTopoSets
        (
            mesh,
            objects,
            subsetter.faceMap(),
            subsetter.subMesh(),
            faceSets
        );
        subsetTopoSets
        (
            mesh,
            objects,
            subsetter.pointMap(),
            subsetter.subMesh(),
            pointSets
        );
    }


    // Write mesh and fields to new time
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (overwrite || specifiedInstance)
    {
        runTime.setTime(instant(fieldsInstance), 0);
        subsetter.subMesh().setInstance(meshInstance);
        topoSet::setInstance(meshInstance, cellSets);
        topoSet::setInstance(meshInstance, faceSets);
        topoSet::setInstance(meshInstance, pointSets);
    }
    else
    {
        runTime++;
        subsetter.subMesh().setInstance(runTime.timeName());
    }

    Info<< "Writing subsetted mesh and fields to time " << runTime.timeName()
        << endl;
    subsetter.subMesh().write();
    processorMeshes::removeFiles(subsetter.subMesh());


    // Volume fields
    forAll(scalarFlds, i)
    {
        scalarFlds[i].write();
    }
    forAll(vectorFlds, i)
    {
        vectorFlds[i].write();
    }
    forAll(sphTensorFlds, i)
    {
        sphTensorFlds[i].write();
    }
    forAll(symmTensorFlds, i)
    {
        symmTensorFlds[i].write();
    }
    forAll(tensorFlds, i)
    {
        tensorFlds[i].write();
    }

    // Surface fields.
    forAll(surfScalarFlds, i)
    {
        surfScalarFlds[i].write();
    }
    forAll(surfVectorFlds, i)
    {
        surfVectorFlds[i].write();
    }
    forAll(surfSphTensorFlds, i)
    {
        surfSphTensorFlds[i].write();
    }
    forAll(surfSymmTensorFlds, i)
    {
        surfSymmTensorFlds[i].write();
    }
    forAll(surfTensorFlds, i)
    {
        surfTensorFlds[i].write();
    }

    // Point fields
    forAll(pointScalarFlds, i)
    {
        pointScalarFlds[i].write();
    }
    forAll(pointVectorFlds, i)
    {
        pointVectorFlds[i].write();
    }
    forAll(pointSphTensorFlds, i)
    {
        pointSphTensorFlds[i].write();
    }
    forAll(pointSymmTensorFlds, i)
    {
        pointSymmTensorFlds[i].write();
    }
    forAll(pointTensorFlds, i)
    {
        pointTensorFlds[i].write();
    }

    // Dimensioned fields
    forAll(scalarDimFlds, i)
    {
        scalarDimFlds[i].write();
    }
    forAll(vectorDimFlds, i)
    {
        vectorDimFlds[i].write();
    }
    forAll(sphTensorDimFlds, i)
    {
        sphTensorDimFlds[i].write();
    }
    forAll(symmTensorDimFlds, i)
    {
        symmTensorDimFlds[i].write();
    }
    forAll(tensorDimFlds, i)
    {
        tensorDimFlds[i].write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
