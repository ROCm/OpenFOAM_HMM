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
    Foam::vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamConvertVolFields_H
#define vtkPV3FoamConvertVolFields_H

// Foam includes
#include "emptyFvPatchField.H"
#include "vtkPV3FoamConvertPatchFaceField.H"
#include "vtkPV3FoamConvertPatchPointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertVolFields
(
    const fvMesh& mesh,
    const volPointInterpolation& pInterp,
    const IOobjectList& objects,
    vtkDataArraySelection *fieldSelection,
    vtkMultiBlockDataSet* output
)
{
    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            GeometricField<Type, fvPatchField, volMesh>::typeName
        )
    );

    label nFields = fieldSelection->GetNumberOfArrays();

    for (label i=0; i<nFields; i++)
    {
        if (fieldSelection->GetArraySetting(i))
        {
            const word fieldName = fieldSelection->GetArrayName(i);

            if (fieldObjects.found(fieldName))
            {
                if (debug)
                {
                    Info<< "converting Foam volume field: " << fieldName
                        << endl;
                }

                GeometricField<Type, fvPatchField, volMesh> tf
                (
                    IOobject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );

                tmp<GeometricField<Type, pointPatchField, pointMesh> > tptf
                (
                    pInterp.interpolate(tf)
                );

                // Convert internal mesh
                if (selectedRegions_[0])
                {
                    convertVolField(tf, output);
                    convertPointField(tptf(), tf, output);
                }

                // Convert patches
                const label regionStartId = idRegionPatches_;
                const label regionEndId = idRegionPatches_
                    + patchDataSize_ - 1;

                for (label i=regionStartId; i<=regionEndId; i++)
                {
                    if (selectedRegions_[i])
                    {

                        const word regionName = reader_->GetRegionSelection()
                            ->GetArrayName(i);
                        const label patchId = mesh.boundaryMesh()
                            .findPatchID(regionName);

                        const fvPatchField<Type>& ptf
                        (
                            tf.boundaryField()[patchId]
                        );

                        if (!isType<emptyFvPatchField<Type> >(ptf))
                        {
                            convertPatchFaceField
                            (
                                tf.name(),
                                ptf,
                                selectedRegionDatasetIds_[i],
                                output
                            );

                            convertPatchPointField
                            (
                                tptf().name(),
                                tptf().boundaryField()[patchId]
                                    .patchInternalField()(),
                                selectedRegionDatasetIds_[i],
                                output
                            );
                        }
                        else
                        {
                            fvPatch p
                            (
                                ptf.patch().patch(),
                                tf.mesh().boundary()
                            );

                            convertPatchFaceField
                            (
                                tf.name(),
                                fvPatchField<Type>(p, tf)
                                    .patchInternalField()(),
                                selectedRegionDatasetIds_[i],
                                output
                            );

                            convertPatchPointField
                            (
                                tptf().name(),
                                tptf().boundaryField()[patchId]
                                   .patchInternalField()(),
                                selectedRegionDatasetIds_[i],
                                output
                            );
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::vtkPV3Foam::convertVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output
)
{
    vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::SafeDownCast
    (
        output->GetDataSet(VOLUME, 0)
    );

    vtkFloatArray *cellTypes = vtkFloatArray::New();
    cellTypes->SetNumberOfTuples(superCells_.size());
    cellTypes->SetNumberOfComponents(Type::nComponents);
    cellTypes->Allocate(Type::nComponents*superCells_.size());
    cellTypes->SetName(tf.name().c_str());

    if (debug)
    {
        Info<< "converting vol<Type>Field: " << tf.name() << nl
            << "nTuples = " << superCells_.size() << nl
            << "nComps  = " << Type::nComponents << endl;
    }

    float vec[Type::nComponents];

    forAll(superCells_, sci)
    {
        const Type& t = tf[superCells_[sci]];
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = t[d];
        }

        cellTypes->InsertTuple(sci, vec);
    }

    internalMesh->GetCellData()->AddArray(cellTypes);
    cellTypes->Delete();
}


template<>
void Foam::vtkPV3Foam::convertVolField
(
    const GeometricField<scalar, fvPatchField, volMesh>& sf,
    vtkMultiBlockDataSet* output
)
{
    vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::SafeDownCast
    (
        output->GetDataSet(VOLUME, 0)
    );

    vtkFloatArray *cellScalars = vtkFloatArray::New();
    cellScalars->SetNumberOfTuples(superCells_.size());
    cellScalars->SetNumberOfComponents(1);
    cellScalars->Allocate(superCells_.size());
    cellScalars->SetName(sf.name().c_str());

    if (debug)
    {
        Info<< "converting volScalarField: " << sf.name() << nl
            << "sf.size() = " << sf.size() << nl
            << "nTuples = " << superCells_.size() << nl
            << "nComps  = 1" << endl;
    }

    forAll(superCells_, sci)
    {
        cellScalars->InsertComponent(sci, 0, sf[superCells_[sci]]);
    }

    internalMesh->GetCellData()->AddArray(cellScalars);
    if (!internalMesh->GetCellData()->GetScalars())
    {
        internalMesh->GetCellData()->SetScalars(cellScalars);
    }

    cellScalars->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
