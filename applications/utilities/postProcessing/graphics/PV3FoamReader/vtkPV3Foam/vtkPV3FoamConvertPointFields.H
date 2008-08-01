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

#ifndef vtkPV3FoamConvertPointFields_H
#define vtkPV3FoamConvertPointFields_H

// Foam includes
#include "interpolatePointToCell.H"

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

    label nFields = fieldSelection->GetNumberOfArrays();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (label i=0; i<nFields; i++)
    {
        const word fieldName = fieldSelection->GetArrayName(i);

        if
        (
            !fieldSelection->GetArraySetting(i)
         || !fieldObjects.found(fieldName))
        {
            continue;
        }

        if (debug)
        {
            Info<< "Foam::vtkPV3Foam::convertPointFields : "
                << fieldName << endl;
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


        //
        // Convert internal mesh - if activated
        //
        for
        (
            int regionId = selectInfoVolume_.start();
            regionId < selectInfoVolume_.end();
            ++regionId
        )
        {
            const label datasetId = selectedRegionDatasetIds_[regionId];

            if (!selectedRegions_[regionId] || datasetId < 0)
            {
                continue;
            }

            convertPointField
            (
                ptf,
                GeometricField<Type, fvPatchField, volMesh>::null(),
                output,
                selectInfoVolume_,
                datasetId
            );
        }


        //
        // Convert patches - if activated
        //
        label regionId = selectInfoPatches_.start();
        forAll (patches, patchI)
        {
            if (patches[patchI].size())
            {
                const label datasetId = selectedRegionDatasetIds_[regionId];

                if (selectedRegions_[regionId] && datasetId >= 0)
                {

                    convertPatchPointField
                    (
                        fieldName,
                        ptf.boundaryField()[patchI].patchInternalField()(),
                        output,
                        selectInfoPatches_,
                        datasetId
                    );
                }
                regionId++;
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
    const label nComp = pTraits<Type>::nComponents;

    vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(ptf.size() + addPointCellLabels_.size());
    pointData->SetNumberOfComponents(nComp);
    pointData->Allocate(nComp*ptf.size());
    pointData->SetName(tf.name().c_str());

    float vec[nComp];

    forAll(ptf, i)
    {
        for (direction d=0; d<nComp; d++)
        {
            vec[d] = component(ptf[i], d);
        }

        pointData->InsertTuple(i, vec);
    }

    label i = ptf.size();

    if (&tf != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels_, api)
        {
            Type t = tf[addPointCellLabels_[api]];

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

    internalMesh->GetPointData()->AddArray(pointData);
    pointData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
