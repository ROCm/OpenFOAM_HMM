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

#ifndef vtkPV3FoamPointFields_H
#define vtkPV3FoamPointFields_H

// Foam includes
#include "interpolatePointToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertPointFields
(
    const fvMesh& mesh,
    const pointMesh& pMesh,
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    // field subset based on type
    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            GeometricField<Type, pointPatchField, pointMesh>::typeName
        )
    );

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    forAllIter(IOobjectList, fieldObjects, iter)
    {
        word fieldName = iter()->name();

        if (debug)
        {
            Info<< "Foam::vtkPV3Foam::convertPointFields : "
                << fieldName << endl;
        }

        GeometricField<Type, pointPatchField, pointMesh> ptf
        (
            *iter(),
            pMesh
        );


        //
        // Convert internal mesh - if activated
        //
        for
        (
            int regionId = regionInfoVolume_.start();
            regionId < regionInfoVolume_.end();
            ++regionId
        )
        {
            const label datasetNo = regionDataset_[regionId];

            if (!regionStatus_[regionId] || datasetNo < 0)
            {
                continue;
            }

            convertPointField
            (
                ptf,
                GeometricField<Type, fvPatchField, volMesh>::null(),
                output,
                regionInfoVolume_,
                datasetNo
            );
        }

        //
        // Convert patches - if activated
        //
        for
        (
            int regionId = regionInfoPatches_.start();
            regionId < regionInfoPatches_.end();
            ++regionId
        )
        {
            word patchName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            const label datasetNo = regionDataset_[regionId];
            const label patchId = patches.findPatchID(patchName);

            if (!regionStatus_[regionId] || datasetNo < 0 || patchId < 0)
            {
                continue;
            }

            convertPatchPointField
            (
                fieldName,
                ptf.boundaryField()[patchId].patchInternalField()(),
                output,
                regionInfoPatches_,
                datasetNo
            );
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
    const label nComp = pTraits<Type>::nComponents;

    vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples( ptf.size() + addPointCellLabels_.size() );
    pointData->SetNumberOfComponents( nComp );
    pointData->Allocate( nComp*(ptf.size() + addPointCellLabels_.size()) );
    pointData->SetName( tf.name().c_str() );


    if (debug)
    {
        Info<< "convert convertPointField: "
            << tf.name()
            << " size = " << ptf.size()
            << " nComp=" << nComp
            << " nTuples = " << (ptf.size() + addPointCellLabels_.size())
            <<  endl;
    }

    float vec[nComp];
    forAll(ptf, i)
    {
        const Type& t = ptf[i];
        for (direction d=0; d<nComp; d++)
        {
            vec[d] = component(t, d);
        }

        pointData->InsertTuple(i, vec);
    }

    // continue insertion from here
    label i = ptf.size();

    if (&tf != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels_, api)
        {
            const Type& t = tf[addPointCellLabels_[api]];
            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }

            pointData->InsertTuple(i++, vec);
        }
    }
    else
    {
        forAll(addPointCellLabels_, api)
        {
            Type t = interpolatePointToCell(ptf, addPointCellLabels_[api]);

            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }

            pointData->InsertTuple(i++, vec);
        }
    }

    vtkmesh->GetPointData()->AddArray(pointData);
    pointData->Delete();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
