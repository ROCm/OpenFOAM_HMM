/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

InClass
    vtkPVFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoamFieldTemplates_C
#define vtkPVFoamFieldTemplates_C

// OpenFOAM includes
#include "emptyFvPatchField.H"
#include "wallPolyPatch.H"
#include "faceSet.H"
#include "volPointInterpolation.H"
#include "zeroGradientFvPatchField.H"
#include "interpolatePointToCell.H"
#include "foamPvFields.H"

// vtk includes
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// volume-fields
//

template<class Type>
void Foam::vtkPVFoam::convertVolField
(
    const PtrList<patchInterpolator>& patchInterpList,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    vtkMultiBlockDataSet* output
)
{
    const fvMesh& mesh = fld.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const bool interpField = !patchInterpList.empty();
    const bool extrapPatch = reader_->GetExtrapolatePatches();

    // Interpolated field (demand driven)
    autoPtr<GeometricField<Type, pointPatchField, pointMesh>> ptfPtr;
    if (interpField)
    {
        if (debug)
        {
            Info<< "convertVolField interpolating:" << fld.name() << endl;
        }

        ptfPtr.reset
        (
            volPointInterpolation::New(mesh).interpolate(fld).ptr()
        );
    }

    // Convert activated internalMesh regions
    convertVolFieldBlock
    (
        fld,
        ptfPtr,
        output,
        arrayRangeVolume_,
        regionPolyDecomp_
    );

    // Convert activated cellZones
    convertVolFieldBlock
    (
        fld,
        ptfPtr,
        output,
        arrayRangeCellZones_,
        zonePolyDecomp_
    );

    // Convert activated cellSets
    convertVolFieldBlock
    (
        fld,
        ptfPtr,
        output,
        arrayRangeCellSets_,
        csetPolyDecomp_
    );


    //
    // Convert patches - if activated
    //
    for
    (
        int partId = arrayRangePatches_.start();
        partId < arrayRangePatches_.end();
        ++partId
    )
    {
        const word patchName = getPartName(partId);
        const label datasetNo = partDataset_[partId];
        const label patchId = patches.findPatchID(patchName);

        if (!partStatus_[partId] || patchId < 0)
        {
            continue;
        }

        vtkPolyData* vtkmesh = getDataFromBlock<vtkPolyData>
        (
            output, arrayRangePatches_, datasetNo
        );

        if (!vtkmesh)
        {
            continue;
        }

        const fvPatchField<Type>& ptf = fld.boundaryField()[patchId];

        if
        (
            isType<emptyFvPatchField<Type>>(ptf)
         ||
            (
                extrapPatch
             && !polyPatch::constraintType(patches[patchId].type())
            )
        )
        {
            fvPatch p(ptf.patch().patch(), mesh.boundary());

            tmp<Field<Type>> tpptf
            (
                fvPatchField<Type>(p, fld).patchInternalField()
            );

            vtkFloatArray* cdata = convertFieldToVTK(fld.name(), tpptf());
            vtkmesh->GetCellData()->AddArray(cdata);
            cdata->Delete();

            if (patchId < patchInterpList.size())
            {
                vtkFloatArray* pdata = convertFieldToVTK
                (
                    fld.name(),
                    patchInterpList[patchId].faceToPointInterpolate(tpptf)()
                );

                vtkmesh->GetPointData()->AddArray(pdata);
                pdata->Delete();
            }
        }
        else
        {
            vtkFloatArray* cdata = convertFieldToVTK(fld.name(), ptf);
            vtkmesh->GetCellData()->AddArray(cdata);
            cdata->Delete();

            if (patchId < patchInterpList.size())
            {
                vtkFloatArray* pdata = convertFieldToVTK
                (
                    fld.name(),
                    patchInterpList[patchId].faceToPointInterpolate(ptf)()
                );

                vtkmesh->GetPointData()->AddArray(pdata);
                pdata->Delete();
            }
        }
    }

    //
    // Convert face zones - if activated
    //
    for
    (
        int partId = arrayRangeFaceZones_.start();
        partId < arrayRangeFaceZones_.end();
        ++partId
    )
    {
        const word zoneName = getPartName(partId);
        const label datasetNo = partDataset_[partId];

        if (!partStatus_[partId] || datasetNo < 0)
        {
            continue;
        }

        const faceZoneMesh& zMesh = mesh.faceZones();
        const label zoneId = zMesh.findZoneID(zoneName);

        if (zoneId < 0)
        {
            continue;
        }

        vtkPolyData* vtkmesh = getDataFromBlock<vtkPolyData>
        (
            output, arrayRangeFaceZones_, datasetNo
        );

        if (vtkmesh)
        {
            vtkFloatArray* cdata = convertFaceFieldToVTK
            (
                fld,
                zMesh[zoneId]
            );

            vtkmesh->GetCellData()->AddArray(cdata);
            cdata->Delete();
        }

        // TODO: points
    }

    //
    // Convert face sets - if activated
    //
    for
    (
        int partId = arrayRangeFaceSets_.start();
        partId < arrayRangeFaceSets_.end();
        ++partId
    )
    {
        const word selectName = getPartName(partId);
        const label datasetNo = partDataset_[partId];

        if (!partStatus_[partId])
        {
            continue;
        }

        vtkPolyData* vtkmesh = getDataFromBlock<vtkPolyData>
        (
            output, arrayRangeFaceSets_, datasetNo
        );

        if (!vtkmesh)
        {
            continue;
        }

        const faceSet fSet(mesh, selectName);

        vtkFloatArray* cdata = convertFaceFieldToVTK
        (
            fld,
            fSet.sortedToc()
        );

        vtkmesh->GetCellData()->AddArray(cdata);
        cdata->Delete();

        // TODO: points
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolFields
(
    const fvMesh& mesh,
    const PtrList<patchInterpolator>& patchInterpList,
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    forAllConstIter(IOobjectList, objects, iter)
    {
        // restrict to GeometricField<Type, ...>
        if
        (
            iter()->headerClassName()
         != GeometricField<Type, fvPatchField, volMesh>::typeName
        )
        {
            continue;
        }

        // Load field
        GeometricField<Type, fvPatchField, volMesh> fld
        (
            *iter(),
            mesh
        );

        // Convert
        convertVolField(patchInterpList, fld, output);
    }
}


template<class Type>
void Foam::vtkPVFoam::convertDimFields
(
    const fvMesh& mesh,
    const PtrList<patchInterpolator>& patchInterpList,
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    forAllConstIter(IOobjectList, objects, iter)
    {
        // restrict to DimensionedField<Type, ...>
        if
        (
            iter()->headerClassName()
         != DimensionedField<Type, volMesh>::typeName
        )
        {
            continue;
        }

        // Load field
        DimensionedField<Type, volMesh> dimFld(*iter(), mesh);

        // Construct volField with zero-gradient patch fields

        IOobject io(dimFld);
        io.readOpt() = IOobject::NO_READ;

        PtrList<fvPatchField<Type>> patchFields(mesh.boundary().size());
        forAll(patchFields, patchI)
        {
            patchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    zeroGradientFvPatchField<Type>::typeName,
                    mesh.boundary()[patchI],
                    dimFld
                )
            );
        }

        VolFieldType volFld
        (
            io,
            dimFld.mesh(),
            dimFld.dimensions(),
            dimFld,
            patchFields
        );
        volFld.correctBoundaryConditions();

        convertVolField(patchInterpList, volFld, output);
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolFieldBlock
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const List<polyDecomp>& decompLst
)
{
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const label datasetNo = partDataset_[partId];

        if (!partStatus_[partId])
        {
            continue;
        }

        vtkUnstructuredGrid* vtkmesh =
            getDataFromBlock<vtkUnstructuredGrid>(output, range, datasetNo);

        if (!vtkmesh)
        {
            continue;
        }

        vtkFloatArray* cdata = convertVolFieldToVTK
        (
            fld,
            decompLst[datasetNo]
        );

        vtkmesh->GetCellData()->AddArray(cdata);
        cdata->Delete();

        if (ptfPtr.valid())
        {
            convertPointField(vtkmesh, ptfPtr(), fld, decompLst[datasetNo]);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// point-fields
//

template<class Type>
void Foam::vtkPVFoam::convertPointFields
(
    const pointMesh& pMesh,
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    const polyMesh& mesh = pMesh.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAllConstIter(IOobjectList, objects, iter)
    {
        const word& fieldName = iter()->name();
        // restrict to this GeometricField<Type, ...>
        if
        (
            iter()->headerClassName()
         != GeometricField<Type, pointPatchField, pointMesh>::typeName
        )
        {
            continue;
        }

        if (debug)
        {
            Info<< "convertPointFields : " << fieldName << endl;
        }

        GeometricField<Type, pointPatchField, pointMesh> pfld(*iter(), pMesh);


        // Convert activated internalMesh regions
        convertPointFieldBlock
        (
            pfld,
            output,
            arrayRangeVolume_,
            regionPolyDecomp_
        );

        // Convert activated cellZones
        convertPointFieldBlock
        (
            pfld,
            output,
            arrayRangeCellZones_,
            zonePolyDecomp_
        );

        // Convert activated cellSets
        convertPointFieldBlock
        (
            pfld,
            output,
            arrayRangeCellSets_,
            csetPolyDecomp_
        );


        //
        // Convert patches - if activated
        //
        for
        (
            int partId = arrayRangePatches_.start();
            partId < arrayRangePatches_.end();
            ++partId
        )
        {
            const word patchName = getPartName(partId);
            const label datasetNo = partDataset_[partId];
            const label patchId = patches.findPatchID(patchName);

            if (!partStatus_[partId] || patchId < 0)
            {
                continue;
            }

            vtkPolyData* vtkmesh = getDataFromBlock<vtkPolyData>
            (
                output, arrayRangePatches_, datasetNo
            );

            if (vtkmesh)
            {
                vtkFloatArray* pdata = convertFieldToVTK
                (
                    fieldName,
                    pfld.boundaryField()[patchId].patchInternalField()()
                );

                vtkmesh->GetPointData()->AddArray(pdata);
                pdata->Delete();
            }
        }

        //
        // Convert faceZones - if activated
        //
        for
        (
            int partId = arrayRangeFaceZones_.start();
            partId < arrayRangeFaceZones_.end();
            ++partId
        )
        {
            const word zoneName = getPartName(partId);
            const label datasetNo = partDataset_[partId];
            const label zoneId = mesh.faceZones().findZoneID(zoneName);

            if (!partStatus_[partId] || zoneId < 0)
            {
                continue;
            }

            vtkPolyData* vtkmesh = getDataFromBlock<vtkPolyData>
            (
                output, arrayRangeFaceZones_, datasetNo
            );

            if (vtkmesh)
            {
                // Extract the field on the zone
                Field<Type> znfld
                (
                    pfld.primitiveField(),
                    mesh.faceZones()[zoneId]().meshPoints()
                );

                vtkFloatArray* pdata = convertFieldToVTK(fieldName, znfld);

                vtkmesh->GetPointData()->AddArray(pdata);
                pdata->Delete();
            }
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertPointFieldBlock
(
    const GeometricField<Type, pointPatchField, pointMesh>& pfld,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const List<polyDecomp>& decompLst
)
{
   for (int partId = range.start(); partId < range.end(); ++partId)
   {
       const label datasetNo = partDataset_[partId];

       if (!partStatus_[partId])
       {
           continue;
       }

       vtkUnstructuredGrid* vtkmesh = getDataFromBlock<vtkUnstructuredGrid>
       (
           output, range, datasetNo
       );

       if (vtkmesh)
       {
           convertPointField
           (
               vtkmesh,
               pfld,
               GeometricField<Type, fvPatchField, volMesh>::null(),
               decompLst[datasetNo]
           );
       }
   }
}


template<class Type>
void Foam::vtkPVFoam::convertPointField
(
    vtkUnstructuredGrid* vtkmesh,
    const GeometricField<Type, pointPatchField, pointMesh>& pfld,
    const GeometricField<Type, fvPatchField, volMesh>& vfld,
    const polyDecomp& decomp
)
{
    if (!vtkmesh)
    {
        return;
    }

    const label nComp = pTraits<Type>::nComponents;
    const labelUList& addPointCellLabels = decomp.addPointCellLabels();
    const labelUList& pointMap = decomp.pointMap();

    // use a pointMap or address directly into mesh
    const label nPoints = (pointMap.size() ? pointMap.size() : pfld.size());

    vtkFloatArray* fldData = vtkFloatArray::New();
    fldData->SetNumberOfTuples(nPoints + addPointCellLabels.size());
    fldData->SetNumberOfComponents(nComp);
    fldData->Allocate(nComp*(nPoints + addPointCellLabels.size()));

    // Note: using the name of the original volField
    // not the name generated by the interpolation "volPointInterpolate(<name>)"

    if (&vfld != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        fldData->SetName(vfld.name().c_str());
    }
    else
    {
        fldData->SetName(pfld.name().c_str());
    }

    if (debug)
    {
        Info<< "convert Point field: "
            << pfld.name()
            << " size="  << (nPoints + addPointCellLabels.size())
            << " (" << nPoints << " + " << addPointCellLabels.size()
            << ") nComp=" << nComp << endl;
    }

    float vec[nComp];

    if (pointMap.size())
    {
        forAll(pointMap, i)
        {
            const Type& t = pfld[pointMap[i]];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            foamPvFields::remapTuple<Type>(vec);

            fldData->InsertTuple(i, vec);
        }
    }
    else
    {
        forAll(pfld, i)
        {
            const Type& t = pfld[i];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            foamPvFields::remapTuple<Type>(vec);

            fldData->InsertTuple(i, vec);
        }
    }

    // continue insertion from here
    label i = nPoints;

    if (&vfld != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels, apI)
        {
            const Type& t = vfld[addPointCellLabels[apI]];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            foamPvFields::remapTuple<Type>(vec);

            fldData->InsertTuple(i++, vec);
        }
    }
    else
    {
        forAll(addPointCellLabels, apI)
        {
            Type t = interpolatePointToCell(pfld, addPointCellLabels[apI]);
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            foamPvFields::remapTuple<Type>(vec);

            fldData->InsertTuple(i++, vec);
        }
    }

    vtkmesh->GetPointData()->AddArray(fldData);
    fldData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// lagrangian-fields
//

template<class Type>
void Foam::vtkPVFoam::convertLagrangianFields
(
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output,
    const label datasetNo
)
{
    const arrayRange& range = arrayRangeLagrangian_;

    forAllConstIter(IOobjectList, objects, iter)
    {
        // restrict to this IOField<Type>
        if (iter()->headerClassName() == IOField<Type>::typeName)
        {
            vtkPolyData* vtkmesh =
                getDataFromBlock<vtkPolyData>(output, range, datasetNo);

            if (vtkmesh)
            {
                IOField<Type> fld(*iter());

                vtkFloatArray* fldData = convertFieldToVTK(fld.name(), fld);
                vtkmesh->GetPointData()->AddArray(fldData);
                fldData->Delete();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// low-level conversions
//

template<class Type>
vtkFloatArray* Foam::vtkPVFoam::convertFieldToVTK
(
    const word& name,
    const Field<Type>& fld
)
{
    if (debug)
    {
        Info<< "convert Field<" << pTraits<Type>::typeName << "> "
            << name
            << " size="  << fld.size()
            << " nComp=" << label(pTraits<Type>::nComponents) << endl;
    }

    const label nComp = pTraits<Type>::nComponents;

    vtkFloatArray* fldData = vtkFloatArray::New();
    fldData->SetNumberOfTuples(fld.size());
    fldData->SetNumberOfComponents(nComp);
    fldData->Allocate(nComp*fld.size());
    fldData->SetName(name.c_str());

    float vec[nComp];
    forAll(fld, i)
    {
        const Type& t = fld[i];
        for (direction d=0; d<nComp; ++d)
        {
            vec[d] = component(t, d);
        }
        foamPvFields::remapTuple<Type>(vec);

        fldData->InsertTuple(i, vec);
    }

    return fldData;
}


template<class Type>
vtkFloatArray* Foam::vtkPVFoam::convertFaceFieldToVTK
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const labelUList& faceLabels
)
{
    if (debug)
    {
        Info<< "convert face field: "
            << fld.name()
            << " size="  << faceLabels.size()
            << " nComp=" << label(pTraits<Type>::nComponents) << endl;
    }

    const fvMesh& mesh = fld.mesh();

    const label nComp = pTraits<Type>::nComponents;
    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray* fldData = vtkFloatArray::New();
    fldData->SetNumberOfTuples(faceLabels.size());
    fldData->SetNumberOfComponents(nComp);
    fldData->Allocate(nComp*faceLabels.size());
    fldData->SetName(fld.name().c_str());

    float vec[nComp];

    // for interior faces: average owner/neighbour
    // for boundary faces: owner
    forAll(faceLabels, facei)
    {
        const label faceNo = faceLabels[facei];
        if (faceNo < nInternalFaces)
        {
            Type t = 0.5*(fld[faceOwner[faceNo]] + fld[faceNeigh[faceNo]]);

            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
        }
        else
        {
            const Type& t = fld[faceOwner[faceNo]];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
        }
        foamPvFields::remapTuple<Type>(vec);

        fldData->InsertTuple(facei, vec);
    }

    return fldData;
}


template<class Type>
vtkFloatArray* Foam::vtkPVFoam::convertVolFieldToVTK
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const polyDecomp& decompInfo
)
{
    const label nComp = pTraits<Type>::nComponents;
    const labelList& superCells = decompInfo.superCells();

    vtkFloatArray* fldData = vtkFloatArray::New();
    fldData->SetNumberOfTuples(superCells.size());
    fldData->SetNumberOfComponents(nComp);
    fldData->Allocate(nComp*superCells.size());
    fldData->SetName(fld.name().c_str());

    if (debug)
    {
        Info<< "convert volField: "
            << fld.name()
            << " size=" << superCells.size()
            << " (" << fld.size() << " + "
            << (superCells.size() - fld.size())
            << ") nComp=" << nComp << endl;
    }

    float vec[nComp];
    forAll(superCells, i)
    {
        const Type& t = fld[superCells[i]];
        for (direction d=0; d<nComp; ++d)
        {
            vec[d] = component(t, d);
        }
        foamPvFields::remapTuple<Type>(vec);

        fldData->InsertTuple(i, vec);
    }

    return fldData;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
