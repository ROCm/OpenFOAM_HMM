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

#ifndef vtkPV3FoamConvertPatchFaceField_H
#define vtkPV3FoamConvertPatchFaceField_H

// VTK includes
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertPatchFaceField
(
    const word& name,
    const Field<Type>& ptf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(ptf.size());
    cellData->SetNumberOfComponents(Type::nComponents);
    cellData->Allocate(Type::nComponents*ptf.size());
    cellData->SetName(name.c_str());

    float vec[Type::nComponents];

    forAll(ptf, i)
    {
        const Type& t = ptf[i];
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = t[d];
        }

        cellData->InsertTuple(i, vec);
    }

    vtkmesh->GetCellData()->AddArray(cellData);
    cellData->Delete();
}


template<>
void Foam::vtkPV3Foam::convertPatchFaceField
(
    const word& name,
    const Field<scalar>& psf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(psf.size());
    cellData->SetNumberOfComponents(1);
    cellData->Allocate(psf.size());
    cellData->SetName(name.c_str());

    forAll(psf, i)
    {
        cellData->InsertComponent(i, 0, psf[i]);
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
