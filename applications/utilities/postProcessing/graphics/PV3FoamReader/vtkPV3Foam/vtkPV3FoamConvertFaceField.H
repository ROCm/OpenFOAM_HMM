/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#ifndef vtkPV3FoamConvertFaceField_H
#define vtkPV3FoamConvertFaceField_H

// VTK includes
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertFaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo,
    const fvMesh& mesh,
    const labelList& faceLabels
)
{
    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(faceLabels.size());
    cellData->SetNumberOfComponents(Type::nComponents);
    cellData->Allocate(Type::nComponents*faceLabels.size());
    cellData->SetName(tf.name().c_str());

    float vec[Type::nComponents];

    forAll(faceLabels, faceI)
    {
        const label faceNo = faceLabels[faceI];
        if (faceNo < nInternalFaces)
        {
            Type t = 0.5 *
            (
                tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]]
            );

            for (direction d=0; d<Type::nComponents; d++)
            {
                vec[d] = t[d];
            }
        }
        else
        {
            const Type& t = tf[faceOwner[faceNo]];
            for (direction d=0; d<Type::nComponents; d++)
            {
                vec[d] = t[d];
            }
        }

        cellData->InsertTuple(faceI, vec);
    }

    vtkmesh->GetCellData()->AddArray(cellData);
    cellData->Delete();
}


template<class Type>
void Foam::vtkPV3Foam::convertFaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo,
    const fvMesh& mesh,
    const faceSet& fSet
)
{
    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(fSet.size());
    cellData->SetNumberOfComponents(Type::nComponents);
    cellData->Allocate(Type::nComponents*fSet.size());
    cellData->SetName(tf.name().c_str());

    float vec[Type::nComponents];

    label faceI = 0;
    forAllConstIter(faceSet, fSet, iter)
    {
        const label faceNo = iter.key();

        if (faceNo < nInternalFaces)
        {
            Type t = 0.5 *
            (
                tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]]
            );

            for (direction d=0; d<Type::nComponents; d++)
            {
                vec[d] = t[d];
            }
        }
        else
        {
            const Type& t = tf[faceOwner[faceNo]];
            for (direction d=0; d<Type::nComponents; d++)
            {
                vec[d] = t[d];
            }
        }

        cellData->InsertTuple(faceI, vec);
        ++faceI;
    }

    vtkmesh->GetCellData()->AddArray(cellData);
    cellData->Delete();
}


template<>
void Foam::vtkPV3Foam::convertFaceField
(
    const GeometricField<scalar, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo,
    const fvMesh& mesh,
    const labelList& faceLabels
)
{
    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(faceLabels.size());
    cellData->SetNumberOfComponents(1);
    cellData->Allocate(faceLabels.size());
    cellData->SetName(tf.name().c_str());

    forAll(faceLabels, faceI)
    {
        const label faceNo = faceLabels[faceI];
        if (faceNo < nInternalFaces)
        {
            cellData->InsertComponent
            (
                faceI,
                0,
                0.5 * (tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]])
            );
        }
        else
        {
            cellData->InsertComponent
            (
                faceI, 0, tf[faceOwner[faceNo]]
            );
        }

    }

    vtkmesh->GetCellData()->AddArray(cellData);

    if (!vtkmesh->GetCellData()->GetScalars())
    {
        vtkmesh->GetCellData()->SetScalars(cellData);
    }

    cellData->Delete();
}


template<>
void Foam::vtkPV3Foam::convertFaceField
(
    const GeometricField<scalar, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo,
    const fvMesh& mesh,
    const faceSet& fSet
)
{
    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(fSet.size());
    cellData->SetNumberOfComponents(1);
    cellData->Allocate(fSet.size());
    cellData->SetName(tf.name().c_str());

    label faceI = 0;
    forAllConstIter(faceSet, fSet, iter)
    {
        const label faceNo = iter.key();
        if (faceNo < nInternalFaces)
        {
            cellData->InsertComponent
            (
                faceI,
                0,
                0.5 * (tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]])
            );
        }
        else
        {
            cellData->InsertComponent
            (
                faceI, 0, tf[faceOwner[faceNo]]
            );
        }
        ++faceI;
    }

    vtkmesh->GetCellData()->AddArray(cellData);

    if (!vtkmesh->GetCellData()->GetScalars())
    {
        vtkmesh->GetCellData()->SetScalars(cellData);
    }

    cellData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
