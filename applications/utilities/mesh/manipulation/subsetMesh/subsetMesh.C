/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Create a mesh subset for a particular region of interest based on a
    cellSet or cellZone.

    See setSet/topoSet utilities on how to define select cells based on
    various shapes.

    Will subset all points, faces and cells needed to make a sub-mesh,
    but not preserve attached boundary types.

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

// Get the exposed patchId or define the exposedPatchName in fvMeshSubset
label getExposedPatchId(const polyMesh& mesh, const word& patchName)
{
    const label patchId = mesh.boundaryMesh().findPatchID(patchName);

    if (patchId == -1)
    {
        fvMeshSubset::exposedPatchName = patchName;
    }

    Info<< "Adding exposed internal faces to "
        << (patchId == -1 ? "new" : "existing")
        << " patch \"" << patchName << "\"" << nl << endl;

    return patchId;
}


labelList nearestPatch(const polyMesh& mesh, const labelList& patchIDs)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Count number of faces in exposedPatchIDs
    label nFaces = 0;
    for (const label patchi : patchIDs)
    {
        nFaces += pbm[patchi].size();
    }

    // Field on cells and faces.
    List<topoDistanceData<label>> cellData(mesh.nCells());
    List<topoDistanceData<label>> faceData(mesh.nFaces());

    // Start of changes
    labelList patchFaces(nFaces);
    List<topoDistanceData<label>> patchData(nFaces);
    nFaces = 0;
    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, i)
        {
            patchFaces[nFaces] = pp.start()+i;
            patchData[nFaces] = topoDistanceData<label>(0, patchi);
            ++nFaces;
        }
    }

    // Propagate information inwards
    FaceCellWave<topoDistanceData<label>> deltaCalc
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

    label nFields = 0;
    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Subsetting " << fieldType << " (";
        }
        else
        {
            Info<< ' ';
        }
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

        subFields.set(nFields, subsetter.interpolate(fld));

        // Subsetting adds 'subset' prefix - rename to match original.
        subFields[nFields].rename(fieldName);

        ++nFields;
    }

    if (nFields)
    {
        Info<< ')' << nl;
    }
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

    label nFields = 0;
    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Subsetting " << fieldType << " (";
        }
        else
        {
            Info<< ' ';
        }
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

        subFields.set(nFields, subsetter.interpolate(fld));

        // Subsetting adds 'subset' prefix - rename to match original.
        subFields[nFields].rename(fieldName);

        ++nFields;
    }

    if (nFields)
    {
        Info<< ')' << nl;
    }
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

    label nFields = 0;
    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Subsetting " << fieldType << " (";
        }
        else
        {
            Info<< ' ';
        }
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

        subFields.set(nFields, subsetter.interpolate(fld));

        // Subsetting adds 'subset' prefix - rename to match original.
        subFields[nFields].rename(fieldName);

        ++nFields;
    }

    if (nFields)
    {
        Info<< ')' << nl;
    }
}


template<class TopoSet>
void subsetTopoSets
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const labelList& map,
    const fvMesh& subMesh,
    PtrList<TopoSet>& subSets
)
{
    // Read original sets
    PtrList<TopoSet> sets;
    ReadFields<TopoSet>(objects, sets);

    subSets.setSize(sets.size());
    forAll(sets, seti)
    {
        const TopoSet& set = sets[seti];

        Info<< "Subsetting " << set.type() << " " << set.name() << endl;

        labelHashSet subset(2*min(set.size(), map.size()));

        forAll(map, i)
        {
            if (set.found(map[i]))
            {
                subset.insert(i);
            }
        }

        subSets.set
        (
            seti,
            new TopoSet
            (
                subMesh,
                set.name(),
                std::move(subset),
                IOobject::AUTO_WRITE
            )
        );
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create a mesh subset for a particular region of interest based on a"
        " cellSet or cellZone(s) specified as the first command argument.\n"
        "See setSet/topoSet utilities on how to select cells based on"
        " various shapes."
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    argList::addArgument
    (
        "cell-selection",
        "The cellSet name, but with the -zone option this is interpreted"
        " to be a cellZone selection by name(s) or regex.\n"
        "Eg 'mixer' or '( mixer \"moving.*\" )'"
    );

    argList::addOption
    (
        "patch",
        "name",
        "Add exposed internal faces to specified patch"
        " instead of \"oldInternalFaces\""
    );
    argList::addOption
    (
        "patches",
        "wordRes",
        "Add exposed internal faces to closest of specified patches"
        " instead of \"oldInternalFaces\""
    );
    argList::addBoolOption
    (
        "zone",
        "Subset with cellZone(s) instead of cellSet."
        " The command argument may be a list of words or regexs"
    );
    argList::addOption
    (
        "resultTime",
        "time",
        "Specify a time for the resulting mesh"
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"

    #include "createNamedMesh.H"

    // arg[1] = word (cellSet) or wordRes (cellZone)
    // const word selectionName = args[1];

    word meshInstance = mesh.pointsInstance();
    word fieldsInstance = runTime.timeName();

    const bool useCellZone = args.found("zone");
    const bool overwrite = args.found("overwrite");
    const bool specifiedInstance = args.readIfPresent
    (
        "resultTime",
        fieldsInstance
    );
    if (specifiedInstance)
    {
        // Set both mesh and field to this time
        meshInstance = fieldsInstance;
    }


    // Default exposed patch id
    labelList exposedPatchIDs(one{}, -1);

    if (args.found("patches"))
    {
        const wordRes patchNames(args.getList<wordRe>("patches"));

        if (patchNames.size() == 1 && patchNames.first().isLiteral())
        {
            exposedPatchIDs.first() =
                getExposedPatchId(mesh, patchNames.first());
        }
        else
        {
            exposedPatchIDs =
                mesh.boundaryMesh().patchSet(patchNames).sortedToc();

            Info<< "Adding exposed internal faces to nearest of patches "
                << flatOutput(patchNames) << nl << endl;

            if (exposedPatchIDs.empty())
            {
                FatalErrorInFunction
                    << nl << "No patches matched. Patches: "
                    << mesh.boundaryMesh().names() << nl
                    << exit(FatalError);
            }
        }
    }
    else if (args.found("patch"))
    {
        exposedPatchIDs.first() =
            getExposedPatchId(mesh, args.get<word>("patch"));
    }
    else
    {
        Info<< "Adding exposed internal faces to patch \""
            << fvMeshSubset::exposedPatchName
            << "\" (created if necessary)" << nl
            << nl;
    }


    autoPtr<cellSet> cellSetPtr;

    // arg[1] can be a word (for cellSet) or wordRes (for cellZone)

    wordRes zoneNames;
    if (useCellZone)
    {
        wordRes selectionNames(args.getList<wordRe>(1));
        zoneNames.transfer(selectionNames);

        Info<< "Using cellZone " << flatOutput(zoneNames) << nl << endl;

        if (mesh.cellZones().findIndex(zoneNames) == -1)
        {
            FatalErrorInFunction
                << "No cellZones found: " << flatOutput(zoneNames) << nl << nl
                << exit(FatalError);
        }
    }
    else
    {
        const word selectionName = args[1];

        Info<< "Using cellSet " << selectionName << nl << endl;

        cellSetPtr = autoPtr<cellSet>::New(mesh, selectionName);
    }


    // Mesh subsetting engine
    fvMeshSubset subsetter(mesh);

    {
        bitSet selectedCells =
        (
            cellSetPtr
          ? BitSetOps::create(mesh.nCells(), *cellSetPtr)
          : mesh.cellZones().selection(zoneNames)
        );

        if (exposedPatchIDs.size() == 1)
        {
            // Single patch for exposed faces
            subsetter.setCellSubset
            (
                selectedCells,
                exposedPatchIDs.first(),
                true
            );
        }
        else
        {
            // The nearest patch per face
            labelList nearestExposedPatch(nearestPatch(mesh, exposedPatchIDs));

            labelList exposedFaces
            (
                subsetter.getExposedFaces(selectedCells, true)
            );

            subsetter.setCellSubset
            (
                selectedCells,
                exposedFaces,
                labelUIndList(nearestExposedPatch, exposedFaces)(),
                true
            );
        }

        Info<< "Subset "
            << returnReduce(subsetter.subMesh().nCells(), sumOp<label>())
            << " of "
            << returnReduce(mesh.nCells(), sumOp<label>())
            << " cells" << nl << nl;
    }


    IOobjectList objects(mesh, runTime.timeName());
    HashTable<wordHashSet> availableFields = objects.classes();


    // Read vol fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<volScalarField> vScalarFlds;
    subsetFields(subsetter, availableFields, vScalarFlds);

    PtrList<volVectorField> vVectorFlds;
    subsetFields(subsetter, availableFields, vVectorFlds);

    PtrList<volSphericalTensorField> vSphTensorFlds;
    subsetFields(subsetter, availableFields, vSphTensorFlds);

    PtrList<volSymmTensorField> vSymmTensorFlds;
    subsetFields(subsetter, availableFields, vSymmTensorFlds);

    PtrList<volTensorField> vTensorFlds;
    subsetFields(subsetter, availableFields, vTensorFlds);


    // Read surface fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<surfaceScalarField> sScalarFlds;
    subsetFields(subsetter, availableFields, sScalarFlds);

    PtrList<surfaceVectorField> sVectorFlds;
    subsetFields(subsetter, availableFields, sVectorFlds);

    PtrList<surfaceSphericalTensorField> sSphTensorFlds;
    subsetFields(subsetter, availableFields, sSphTensorFlds);

    PtrList<surfaceSymmTensorField> sSymmTensorFlds;
    subsetFields(subsetter, availableFields, sSymmTensorFlds);

    PtrList<surfaceTensorField> sTensorFlds;
    subsetFields(subsetter, availableFields, sTensorFlds);


    // Read point fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const pointMesh& pMesh = pointMesh::New(mesh);

    PtrList<pointScalarField> pScalarFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pScalarFlds);

    PtrList<pointVectorField> pVectorFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pVectorFlds);

    PtrList<pointSphericalTensorField> pSphTensorFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pSphTensorFlds);

    PtrList<pointSymmTensorField> pSymmTensorFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pSymmTensorFlds);

    PtrList<pointTensorField> pTensorFlds;
    subsetPointFields(subsetter, pMesh, availableFields, pTensorFlds);


    // Read dimensioned fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<volScalarField::Internal> dScalarFlds;
    subsetDimensionedFields(subsetter, availableFields, dScalarFlds);

    PtrList<volVectorField::Internal> dVectorFlds;
    subsetDimensionedFields(subsetter, availableFields, dVectorFlds);

    PtrList<volSphericalTensorField::Internal> dSphTensorFlds;
    subsetDimensionedFields(subsetter, availableFields, dSphTensorFlds);

    PtrList<volSymmTensorField::Internal> dSymmTensorFlds;
    subsetDimensionedFields(subsetter, availableFields, dSymmTensorFlds);

    PtrList<volTensorField::Internal> dTensorFlds;
    subsetDimensionedFields(subsetter, availableFields, dTensorFlds);


    // Read topoSets and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<cellSet> cellSets;
    PtrList<faceSet> faceSets;
    PtrList<pointSet> pointSets;

    {
        IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
        if (cellSetPtr)
        {
            objects.remove(*cellSetPtr);
        }
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
        ++runTime;
        subsetter.subMesh().setInstance(runTime.timeName());
    }

    Info<< "Writing subsetted mesh and fields to time " << runTime.timeName()
        << endl;
    subsetter.subMesh().write();
    processorMeshes::removeFiles(subsetter.subMesh());


    // Volume fields
    for (const auto& fld : vScalarFlds)     { fld.write(); }
    for (const auto& fld : vVectorFlds)     { fld.write(); }
    for (const auto& fld : vSphTensorFlds)  { fld.write(); }
    for (const auto& fld : vSymmTensorFlds) { fld.write(); }
    for (const auto& fld : vTensorFlds)     { fld.write(); }

    // Surface fields
    for (const auto& fld : sScalarFlds)     { fld.write(); }
    for (const auto& fld : sVectorFlds)     { fld.write(); }
    for (const auto& fld : sSphTensorFlds)  { fld.write(); }
    for (const auto& fld : sSymmTensorFlds) { fld.write(); }
    for (const auto& fld : sTensorFlds)     { fld.write(); }

    // Point fields
    for (const auto& fld : pScalarFlds)     { fld.write(); }
    for (const auto& fld : pVectorFlds)     { fld.write(); }
    for (const auto& fld : pSphTensorFlds)  { fld.write(); }
    for (const auto& fld : pSymmTensorFlds) { fld.write(); }
    for (const auto& fld : pTensorFlds)     { fld.write(); }

    // Dimensioned fields
    for (const auto& fld : dScalarFlds)     { fld.write(); }
    for (const auto& fld : dVectorFlds)     { fld.write(); }
    for (const auto& fld : dSphTensorFlds)  { fld.write(); }
    for (const auto& fld : dSymmTensorFlds) { fld.write(); }
    for (const auto& fld : dTensorFlds)     { fld.write(); }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
