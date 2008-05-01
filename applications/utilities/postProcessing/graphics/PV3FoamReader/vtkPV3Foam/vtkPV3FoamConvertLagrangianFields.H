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

Class
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamConvertLagrangianFields_H
#define vtkPV3FoamConvertLagrangianFields_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertLagrangianFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    vtkDataArraySelection* fieldSelection,
    vtkMultiBlockDataSet* output
)
{
    if (!lagrangianDataSize_)
    {
        return;
    }
    if (!selectedRegions_[LAGRANGIAN])
    {
        return;
    }

    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            IOField<Type>::typeName
        )
    );

    label nSelectedFields = fieldSelection->GetNumberOfArrays();

    for (label i=0; i<nSelectedFields; i++)
    {
        if (fieldSelection->GetArraySetting(i))
        {
            const word fieldName = fieldSelection->GetArrayName(i);

            if (fieldObjects.found(fieldName))
            {
                if (debug)
                {
                    Info<< "converting Foam Lagrangian field: " << fieldName
                        << endl;
                }

                IOField<Type> iotf
                (
                    IOobject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        "lagrangian"/cloudName_,
                        mesh,
                        IOobject::MUST_READ
                    )
                );

                convertLagrangianField(iotf, output);
            }
        }
    }
}


template<class Type>
void Foam::vtkPV3Foam::convertLagrangianField
(
    const IOField<Type>& tf,
    vtkMultiBlockDataSet* output
)
{
    vtkUnstructuredGrid* lagrangianMesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, LAGRANGIAN, 0)
    );

    vtkFloatArray *pointTypes = vtkFloatArray::New();
    pointTypes->SetNumberOfTuples(tf.size());
    pointTypes->SetNumberOfComponents(Type::nComponents);
    pointTypes->Allocate(Type::nComponents*tf.size());
    pointTypes->SetName(tf.name().c_str());

    if (debug)
    {
        Info<< "converting Lagrangian <Type>Field: " << tf.name() << nl
            << "tf.size() = " << tf.size() << nl
            << "nComps  = " << Type::nComponents << endl;
    }

    float vec[Type::nComponents];

    forAll(tf, i)
    {
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = tf[i][d];
        }

        pointTypes->InsertTuple(i, vec);
    }

    lagrangianMesh->GetPointData()->AddArray(pointTypes);
    pointTypes->Delete();
}


template<>
void Foam::vtkPV3Foam::convertLagrangianField
(
    const IOField<scalar>& sf,
    vtkMultiBlockDataSet* output
)
{
    vtkUnstructuredGrid* lagrangianMesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, LAGRANGIAN, 0)
    );

    vtkFloatArray *pointScalars = vtkFloatArray::New();
    pointScalars->SetNumberOfTuples(sf.size());
    pointScalars->SetNumberOfComponents(1);
    pointScalars->Allocate(sf.size());
    pointScalars->SetName(sf.name().c_str());

    if (debug)
    {
        Info<< "converting Lagrangian scalarField: " << sf.name() << nl
            << "sf.size() = " << sf.size() << nl
            << "nComps  = 1" << endl;
    }

    for (int i=0; i<sf.size(); i++)
    {
        pointScalars->InsertComponent(i, 0, sf[i]);
    }

    lagrangianMesh->GetPointData()->AddArray(pointScalars);
    if (!lagrangianMesh->GetPointData()->GetScalars())
    {
        lagrangianMesh->GetPointData()->SetScalars(pointScalars);
    }

    pointScalars->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
