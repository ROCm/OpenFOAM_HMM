/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamConvertPointFields_H
#define vtkPV3FoamConvertPointFields_H

// Foam includes
#include "interpolatePointToCell.H"

// VTK includes

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertPointFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    vtkDataArraySelection* fieldSelection,
    vtkMultiBlockDataSet* output
)
{
    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            GeometricField<Type, pointPatchField, pointMesh>::typeName
        )
    );

    label nSelectedFields = fieldSelection->GetNumberOfArrays();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (label i=0; i<nSelectedFields; i++)
    {
        if (fieldSelection->GetArraySetting(i))
        {
            word fieldName = fieldSelection->GetArrayName(i);

            if (fieldObjects.found(fieldName))
            {

                if (debug)
                {
                    Info<< "converting Foam point field: " << fieldName
                        << endl;
                }

                pointMesh pMesh(mesh);

                GeometricField<Type, pointPatchField, pointMesh> ptf
                (
                    IOobject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    pMesh
                );

                // Convert internal mesh
                for
                (
                    int regionId = selectInfoVolume_.start();
                    regionId < selectInfoVolume_.end();
                    ++regionId
                )
                {
                    if (selectedRegions_[regionId])
                    {
                        convertPointField
                        (
                            ptf,
                            GeometricField<Type, fvPatchField, volMesh>::null(),
                            output,
                            selectInfoVolume_,
                            selectedRegionDatasetIds_[regionId]
                        );
                    }
                }


                // Convert patches
                label regionId = selectInfoPatches_.start();
                forAll (patches, patchI)
                {
                    if (patches[patchI].size())
                    {
                        if (selectedRegions_[regionId])
                        {
                            convertPatchPointField
                            (
                                ptf.name(),
                                ptf.boundaryField()[patchI]
                                   .patchInternalField()(),
                                output,
                                selectInfoPatches_,
                                selectedRegionDatasetIds_[regionId]
                            );
                        }
                        regionId++;
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::vtkPV3Foam::convertPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptf,
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(ptf.size() + addPointCellLabels_.size());
    pointData->SetNumberOfComponents(Type::nComponents);
    pointData->Allocate(Type::nComponents*ptf.size());
    pointData->SetName(ptf.name().c_str());

    float vec[Type::nComponents];

    forAll(ptf, i)
    {
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = ptf[i][d];
        }

        pointData->InsertTuple(i, vec);
    }

    label i = ptf.size();

    if (&tf != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels_, api)
        {
            Type t = tf[addPointCellLabels_[api]];

            for (direction d=0; d<Type::nComponents; d++)
            {
                vec[d] = t[d];
            }

            pointData->InsertTuple(i++, vec);
        }
    }
    else
    {
        forAll(addPointCellLabels_, api)
        {
            Type t = interpolatePointToCell(ptf, addPointCellLabels_[api]);

            for (direction d=0; d<Type::nComponents; d++)
            {
                vec[d] = t[d];
            }

            pointData->InsertTuple(i++, vec);
        }
    }

    internalMesh->GetPointData()->AddArray(pointData);
    pointData->Delete();
}


template<>
void Foam::vtkPV3Foam::convertPointField
(
    const GeometricField<scalar, pointPatchField, pointMesh>& psf,
    const GeometricField<scalar, fvPatchField, volMesh>& sf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(psf.size() + addPointCellLabels_.size());
    pointData->SetNumberOfComponents(1);
    pointData->Allocate(psf.size());
    pointData->SetName(psf.name().c_str());

    for (int i=0; i<psf.size(); i++)
    {
        pointData->InsertComponent(i, 0, psf[i]);
    }

    label i = psf.size();

    if (&sf != &GeometricField<scalar, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels_, api)
        {
            pointData->InsertComponent
            (
                i++,
                0,
                sf[addPointCellLabels_[api]]
            );
        }
    }
    else
    {
        forAll(addPointCellLabels_, api)
        {
            pointData->InsertComponent
            (
                i++,
                0,
                interpolatePointToCell(psf, addPointCellLabels_[api])
            );
        }
    }

    internalMesh->GetPointData()->AddArray(pointData);
    if (!internalMesh->GetPointData()->GetScalars())
    {
        internalMesh->GetPointData()->SetScalars(pointData);
    }

    pointData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
