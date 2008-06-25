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

#ifndef vtkPV3FoamConvertPatchPointField_H
#define vtkPV3FoamConvertPatchPointField_H

// VTK includes
#include "vtkPointData.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertPatchPointField
(
    const word& name,
    const Field<Type>& pptf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(pptf.size());
    pointData->SetNumberOfComponents(Type::nComponents);
    pointData->Allocate(Type::nComponents*pptf.size());
    pointData->SetName(name.c_str());

    float vec[Type::nComponents];

    forAll(pptf, i)
    {
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = pptf[i][d];
        }

        pointData->InsertTuple(i, vec);
    }

    vtkmesh->GetPointData()->AddArray(pointData);
    pointData->Delete();
}


template<>
void Foam::vtkPV3Foam::convertPatchPointField
(
    const word& name,
    const Field<scalar>& ppsf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(ppsf.size());
    pointData->SetNumberOfComponents(1);
    pointData->Allocate(ppsf.size());
    pointData->SetName(name.c_str());

    for (int i=0; i<ppsf.size(); i++)
    {
        pointData->InsertComponent(i, 0, ppsf[i]);
    }

    vtkmesh->GetPointData()->AddArray(pointData);
    if (!vtkmesh->GetPointData()->GetScalars())
    {
        vtkmesh->GetPointData()->SetScalars(pointData);
    }

    pointData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
