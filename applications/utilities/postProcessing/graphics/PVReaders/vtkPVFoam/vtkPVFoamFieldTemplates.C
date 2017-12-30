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
#include "error.H"
#include "emptyFvPatchField.H"
#include "wallPolyPatch.H"
#include "faceSet.H"
#include "volPointInterpolation.H"
#include "zeroGradientFvPatchField.H"
#include "interpolatePointToCell.H"
#include "foamPvFields.H"
#include "areaFaMesh.H"
#include "areaFields.H"

// vtk includes
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// volume-fields
//

template<class Type>
void Foam::vtkPVFoam::convertVolField
(
    const PtrList<patchInterpolator>& patchInterpList,
    const GeometricField<Type, fvPatchField, volMesh>& fld
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
            Info<< "convertVolField interpolating:" << fld.name() << nl;
        }

        ptfPtr.reset
        (
            volPointInterpolation::New(mesh).interpolate(fld).ptr()
        );
    }

    convertVolFieldBlock(fld, ptfPtr, rangeVolume_);    // internalMesh
    convertVolFieldBlock(fld, ptfPtr, rangeCellZones_); // cellZones
    convertVolFieldBlock(fld, ptfPtr, rangeCellSets_);  // cellSets

    // Patches - currently skip field conversion for groups
    for
    (
        const auto partId
      : rangePatches_.intersection(selectedPartIds_)
    )
    {
        const auto& longName = selectedPartIds_[partId];

        auto iter = cachedVtp_.find(longName);
        if (!iter.found() || !iter.object().dataset)
        {
            // Should not happen, but for safety require a vtk geometry
            continue;
        }

        foamVtpData& vtpData = iter.object();
        auto dataset = vtpData.dataset;

        const labelList& patchIds = vtpData.additionalIds();

        if (patchIds.empty())
        {
            continue;
        }

        // This is slightly roundabout, but we deal with groups and with
        // single patches.
        // For groups (spanning several patches) it is fairly messy to
        // get interpolated point fields. We would need to create a indirect
        // patch each time to obtain the mesh points. We thus skip that.
        //
        // To improve code reuse, we allocate the CellData as a zeroed-field
        // ahead of time.

        vtkSmartPointer<vtkFloatArray> cdata = zeroVTKField<Type>
        (
            fld.name(),
            dataset->GetNumberOfPolys()
        );

        vtkSmartPointer<vtkFloatArray> pdata;
        const bool allowPdata = (interpField && patchIds.size() == 1);

        label coffset = 0; // the write offset into cell-data
        for (label patchId : patchIds)
        {
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

                coffset += transcribeFloatData(cdata, tpptf(), coffset);

                if (allowPdata && patchId < patchInterpList.size())
                {
                    pdata = convertFieldToVTK
                    (
                        fld.name(),
                        patchInterpList[patchId].faceToPointInterpolate(tpptf)()
                    );
                }
            }
            else
            {
                coffset += transcribeFloatData(cdata, ptf, coffset);

                if (allowPdata && patchId < patchInterpList.size())
                {
                    pdata = convertFieldToVTK
                    (
                        fld.name(),
                        patchInterpList[patchId].faceToPointInterpolate(ptf)()
                    );
                }
            }
        }

        if (cdata)
        {
            dataset->GetCellData()->AddArray(cdata);
        }
        if (pdata)
        {
            dataset->GetPointData()->AddArray(pdata);
        }
    }


    // Face Zones
    for
    (
        const auto partId
      : rangeFaceZones_.intersection(selectedPartIds_)
    )
    {
        const auto& longName = selectedPartIds_[partId];
        const word zoneName = getFoamName(longName);

        auto iter = cachedVtp_.find(longName);
        if (!iter.found() || !iter.object().dataset)
        {
            // Should not happen, but for safety require a vtk geometry
            continue;
        }

        foamVtpData& vtpData = iter.object();
        auto dataset = vtpData.dataset;

        const faceZoneMesh& zMesh = mesh.faceZones();
        const label zoneId = zMesh.findZoneID(zoneName);

        if (zoneId < 0)
        {
            continue;
        }

        vtkSmartPointer<vtkFloatArray> cdata =
            convertFaceFieldToVTK
            (
                fld,
                zMesh[zoneId]
            );

        dataset->GetCellData()->AddArray(cdata);

        // TODO: point data
    }


    // Face Sets
    for
    (
        const auto partId
      : rangeFaceSets_.intersection(selectedPartIds_)
    )
    {
        const auto& longName = selectedPartIds_[partId];
        const word selectName = getFoamName(longName);

        auto iter = cachedVtp_.find(longName);
        if (!iter.found() || !iter.object().dataset)
        {
            // Should not happen, but for safety require a vtk geometry
            continue;
        }
        foamVtpData& vtpData = iter.object();
        auto dataset = vtpData.dataset;

        vtkSmartPointer<vtkFloatArray> cdata = convertFaceFieldToVTK
        (
            fld,
            vtpData.cellMap()
        );

        dataset->GetCellData()->AddArray(cdata);

        // TODO: point data
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolFields
(
    const fvMesh& mesh,
    const PtrList<patchInterpolator>& patchInterpList,
    const IOobjectList& objects
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    forAllConstIters(objects, iter)
    {
        // Restrict to GeometricField<Type, ...>
        const auto& ioobj = *(iter.object());

        if (ioobj.headerClassName() == FieldType::typeName)
        {
            // Load field
            FieldType fld(ioobj, mesh);

            // Convert
            convertVolField(patchInterpList, fld);
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertDimFields
(
    const fvMesh& mesh,
    const PtrList<patchInterpolator>& patchInterpList,
    const IOobjectList& objects
)
{
    typedef DimensionedField<Type, volMesh> FieldType;
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    forAllConstIters(objects, iter)
    {
        // Restrict to DimensionedField<Type, ...>
        const auto& ioobj = *(iter.object());

        if (ioobj.headerClassName() != FieldType::typeName)
        {
            continue;
        }

        // Load field
        FieldType dimFld(ioobj, mesh);

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

        convertVolField(patchInterpList, volFld);
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolFieldBlock
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr,
    const arrayRange& range
)
{
    for
    (
        const auto partId
      : range.intersection(selectedPartIds_)
    )
    {
        const auto& longName = selectedPartIds_[partId];

        auto iter = cachedVtu_.find(longName);
        if (!iter.found() || !iter.object().dataset)
        {
            // Should not happen, but for safety require a vtk geometry
            continue;
        }

        foamVtuData& vtuData = iter.object();
        auto dataset = vtuData.dataset;

        vtkSmartPointer<vtkFloatArray> cdata = convertVolFieldToVTK
        (
            fld,
            vtuData
        );
        dataset->GetCellData()->AddArray(cdata);

        if (ptfPtr.valid())
        {
            vtkSmartPointer<vtkFloatArray> pdata = convertPointField
            (
                ptfPtr(),
                fld,
                vtuData
            );
            dataset->GetPointData()->AddArray(pdata);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// area-fields
//

template<class Type>
void Foam::vtkPVFoam::convertAreaFields
(
    const faMesh& mesh,
    const IOobjectList& objects
)
{
    typedef GeometricField<Type, faPatchField, areaMesh> FieldType;

    const List<label> partIds =
        rangeArea_.intersection(selectedPartIds_);

    if (partIds.empty())
    {
        return;
    }

    forAllConstIters(objects, iter)
    {
        // Restrict to GeometricField<Type, ...>
        const auto& ioobj = *(iter.object());

        if (ioobj.headerClassName() == FieldType::typeName)
        {
            // Load field
            FieldType fld(ioobj, mesh);

            // Convert
            for (const auto partId : partIds)
            {
                const auto& longName = selectedPartIds_[partId];

                auto iter = cachedVtp_.find(longName);
                if (!iter.found() || !iter.object().dataset)
                {
                    // Should not happen, but for safety require a vtk geometry
                    continue;
                }

                foamVtpData& vtpData = iter.object();
                auto dataset = vtpData.dataset;

                vtkSmartPointer<vtkFloatArray> cdata = convertFieldToVTK
                (
                    fld.name(),
                    fld
                );
                dataset->GetCellData()->AddArray(cdata);
            }
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
    const IOobjectList& objects
)
{
    const polyMesh& mesh = pMesh.mesh();

    typedef GeometricField<Type, pointPatchField, pointMesh> FieldType;

    forAllConstIters(objects, iter)
    {
        // Restrict to this GeometricField<Type, ...>
        const auto& ioobj = *(iter.object());

        const word& fieldName = ioobj.name();
        if (ioobj.headerClassName() != FieldType::typeName)
        {
            continue;
        }

        if (debug)
        {
            Info<< "convertPointFields : " << fieldName << nl;
        }

        FieldType pfld(ioobj, pMesh);

        convertPointFieldBlock(pfld, rangeVolume_);     // internalMesh
        convertPointFieldBlock(pfld, rangeCellZones_);  // cellZones
        convertPointFieldBlock(pfld, rangeCellSets_);   // cellSets

        // Patches - currently skip field conversion for groups
        for
        (
            const auto partId
          : rangePatches_.intersection(selectedPartIds_)
        )
        {
            const auto& longName = selectedPartIds_[partId];

            auto iter = cachedVtp_.find(longName);
            if (!iter.found() || !iter.object().dataset)
            {
                // Should not happen, but for safety require a vtk geometry
                continue;
            }

            foamVtpData& vtpData = iter.object();
            auto dataset = vtpData.dataset;

            const labelList& patchIds = vtpData.additionalIds();
            if (patchIds.size() != 1)
            {
                continue;
            }

            const label patchId = patchIds[0];

            vtkSmartPointer<vtkFloatArray> pdata = convertFieldToVTK
            (
                fieldName,
                pfld.boundaryField()[patchId].patchInternalField()()
            );

            dataset->GetPointData()->AddArray(pdata);
        }

        // Face Zones
        for
        (
            const auto partId
          : rangeFaceZones_.intersection(selectedPartIds_)
        )
        {
            const auto& longName = selectedPartIds_[partId];
            const word zoneName = getFoamName(longName);

            auto iter = cachedVtp_.find(longName);
            if (!iter.found() || !iter.object().dataset)
            {
                // Should not happen, but for safety require a vtk geometry
                continue;
            }

            foamVtpData& vtpData = iter.object();
            auto dataset = vtpData.dataset;

            const label zoneId = mesh.faceZones().findZoneID(zoneName);

            if (zoneId < 0)
            {
                continue;
            }

            // Extract the field on the zone
            Field<Type> znfld
            (
                pfld.primitiveField(),
                mesh.faceZones()[zoneId]().meshPoints()
            );

            vtkSmartPointer<vtkFloatArray> pdata =
                convertFieldToVTK
                (
                    fieldName,
                    znfld
                );

            dataset->GetPointData()->AddArray(pdata);
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertPointFieldBlock
(
    const GeometricField<Type, pointPatchField, pointMesh>& pfld,
    const arrayRange& range
)
{
    for
    (
        const auto partId
      : range.intersection(selectedPartIds_)
    )
    {
        const auto& longName = selectedPartIds_[partId];

        auto iter = cachedVtu_.find(longName);
        if (!iter.found() || !iter.object().dataset)
        {
            // Should not happen, but for safety require a vtk geometry
            continue;
        }

        foamVtuData& vtuData = iter.object();
        auto dataset = vtuData.dataset;

        vtkSmartPointer<vtkFloatArray> pdata = convertPointField
        (
            pfld,
            GeometricField<Type, fvPatchField, volMesh>::null(),
            vtuData
        );

        dataset->GetPointData()->AddArray(pdata);
    }
}


template<class Type>
vtkSmartPointer<vtkFloatArray> Foam::vtkPVFoam::convertPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& pfld,
    const GeometricField<Type, fvPatchField, volMesh>& vfld,
    const foamVtuData& vtuData
)
{
    const int nComp(pTraits<Type>::nComponents);
    const labelUList& addPointCellLabels = vtuData.additionalIds();
    const labelUList& pointMap = vtuData.pointMap();

    // Use a pointMap or address directly into mesh
    const label nPoints = (pointMap.size() ? pointMap.size() : pfld.size());

    auto data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetNumberOfComponents(nComp);
    data->SetNumberOfTuples(nPoints + addPointCellLabels.size());

    // Note: using the name of the original volField
    // not the name generated by the interpolation "volPointInterpolate(<name>)"

    if (notNull(vfld))
    {
        data->SetName(vfld.name().c_str());
    }
    else
    {
        data->SetName(pfld.name().c_str());
    }

    if (debug)
    {
        Info<< "convert Point field: "
            << pfld.name()
            << " size="  << (nPoints + addPointCellLabels.size())
            << " (" << nPoints << " + " << addPointCellLabels.size()
            << ") nComp=" << nComp << nl;
    }

    float vec[nComp];

    label pointi = 0;
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

            data->SetTuple(pointi++, vec);
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

            data->SetTuple(pointi++, vec);
        }
    }

    // Continue with additional points

    if (notNull(vfld))
    {
        forAll(addPointCellLabels, apI)
        {
            const Type& t = vfld[addPointCellLabels[apI]];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            foamPvFields::remapTuple<Type>(vec);

            data->SetTuple(pointi++, vec);
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

            data->SetTuple(pointi++, vec);
        }
    }

    return data;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// lagrangian-fields
//

template<class Type>
void Foam::vtkPVFoam::convertLagrangianFields
(
    const IOobjectList& objects,
    vtkPolyData* vtkmesh
)
{
    forAllConstIters(objects, iter)
    {
        // Restrict to IOField<Type>
        const auto& ioobj = *(iter.object());

        if (ioobj.headerClassName() == IOField<Type>::typeName)
        {
            IOField<Type> fld(ioobj);

            vtkSmartPointer<vtkFloatArray> data =
                convertFieldToVTK
                (
                    ioobj.name(),
                    fld
                );

            // Provide identical data as cell and as point data
            vtkmesh->GetCellData()->AddArray(data);
            vtkmesh->GetPointData()->AddArray(data);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// low-level conversions
//

template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtkPVFoam::zeroVTKField
(
    const word& name,
    const label size
)
{
    auto data = vtkSmartPointer<vtkFloatArray>::New();

    data->SetName(name.c_str());
    data->SetNumberOfComponents(int(pTraits<Type>::nComponents));
    data->SetNumberOfTuples(size);

    data->Fill(0);

    return data;
}


template<class Type>
Foam::label Foam::vtkPVFoam::transcribeFloatData
(
    vtkFloatArray* array,
    const UList<Type>& input,
    const label start
)
{
    const int nComp(pTraits<Type>::nComponents);

    if (array->GetNumberOfComponents() != nComp)
    {
        FatalErrorInFunction
            << "vtk array '" << array->GetName()
            << "' has mismatch in number of components for type '"
            << pTraits<Type>::typeName
            << "' : target array has " << array->GetNumberOfComponents()
            << " components instead of " << nComp
            << nl;
    }

    const vtkIdType maxSize = array->GetNumberOfTuples();
    const vtkIdType endPos = vtkIdType(start) + vtkIdType(input.size());

    if (start < 0 || vtkIdType(start) >= maxSize)
    {
        WarningInFunction
            << "vtk array '" << array->GetName()
            << "' copy with out-of-range (0-" << long(maxSize-1) << ")"
            << " starting location"
            << nl;

        return 0;
    }
    else if (endPos > maxSize)
    {
        WarningInFunction
            << "vtk array '" << array->GetName()
            << "' copy ends out-of-range (" << long(maxSize) << ")"
            << " using sizing (start,size) = ("
            << start << "," << input.size() << ")"
            << nl;

        return 0;
    }

    float scratch[nComp];
    forAll(input, idx)
    {
        const Type& t = input[idx];
        for (direction d=0; d<nComp; ++d)
        {
            scratch[d] = component(t, d);
        }
        foamPvFields::remapTuple<Type>(scratch);

        array->SetTuple(start+idx, scratch);
    }

    return input.size();
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtkPVFoam::convertFieldToVTK
(
    const word& name,
    const UList<Type>& fld
) const
{
    const int nComp(pTraits<Type>::nComponents);

    if (debug)
    {
        Info<< "convert UList<" << pTraits<Type>::typeName << "> "
            << name
            << " size="  << fld.size()
            << " nComp=" << nComp << nl;
    }

    auto data = vtkSmartPointer<vtkFloatArray>::New();

    data->SetName(name.c_str());
    data->SetNumberOfComponents(nComp);
    data->SetNumberOfTuples(fld.size());

    transcribeFloatData(data, fld);

    return data;
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtkPVFoam::convertFaceFieldToVTK
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const labelUList& faceLabels
) const
{
    if (debug)
    {
        Info<< "convert face field: "
            << fld.name()
            << " size="  << faceLabels.size()
            << " nComp=" << int(pTraits<Type>::nComponents) << nl;
    }

    const fvMesh& mesh = fld.mesh();

    const int nComp(pTraits<Type>::nComponents);
    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    auto data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetName(fld.name().c_str());
    data->SetNumberOfComponents(nComp);
    data->SetNumberOfTuples(faceLabels.size());

    float scratch[nComp];

    // Interior faces: average owner/neighbour
    // Boundary faces: the owner value
    forAll(faceLabels, idx)
    {
        const label faceNo = faceLabels[idx];
        if (faceNo < nInternalFaces)
        {
            Type t = 0.5*(fld[faceOwner[faceNo]] + fld[faceNeigh[faceNo]]);

            for (direction d=0; d<nComp; ++d)
            {
                scratch[d] = component(t, d);
            }
        }
        else
        {
            const Type& t = fld[faceOwner[faceNo]];
            for (direction d=0; d<nComp; ++d)
            {
                scratch[d] = component(t, d);
            }
        }
        foamPvFields::remapTuple<Type>(scratch);

        data->SetTuple(idx, scratch);
    }

    return data;
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtkPVFoam::convertVolFieldToVTK
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const foamVtuData& vtuData
) const
{
    const int nComp(pTraits<Type>::nComponents);
    const labelUList& cellMap = vtuData.cellMap();

    auto data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetName(fld.name().c_str());
    data->SetNumberOfComponents(nComp);
    data->SetNumberOfTuples(cellMap.size());

    if (debug)
    {
        Info<< "convert volField: "
            << fld.name()
            << " size=" << cellMap.size()
            << " (" << fld.size() << " + "
            << (cellMap.size() - fld.size())
            << ") nComp=" << nComp << nl;
    }

    float scratch[nComp];
    forAll(cellMap, idx)
    {
        const Type& t = fld[cellMap[idx]];
        for (direction d=0; d<nComp; ++d)
        {
            scratch[d] = component(t, d);
        }
        foamPvFields::remapTuple<Type>(scratch);

        data->SetTuple(idx, scratch);
    }

    return data;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
