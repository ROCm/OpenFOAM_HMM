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

Description

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "IOobjectList.H"
#include "vtkPV3FoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "vtkPV3FoamConvertVolFields.H"
#include "vtkPV3FoamConvertPointFields.H"
#include "vtkPV3FoamConvertLagrangianFields.H"

void Foam::vtkPV3Foam::updateFoamMesh()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateFoamMesh" << endl;
    }

    if
    (
        !reader_->GetCacheMesh()
     || reader_->GetTimeSelection()->GetArraySetting(0)
    )
    {
        delete meshPtr_;
        meshPtr_ = NULL;
    }

    // Check to see if the FOAM mesh has been created
    if (!meshPtr_)
    {
        if (debug)
        {
            Info<< "Creating Foam mesh" << endl;
        }
        meshPtr_ = new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                dbPtr_().timeName(),
                dbPtr_()
            )
        );
    }
    else
    {
        if (debug)
        {
            Info<< "Using existing Foam mesh" << endl;
        }
    }
}


void Foam::vtkPV3Foam::updateVolFields
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateVolFields" << endl;
    }

    const fvMesh& mesh = *meshPtr_;

    // Construct interpolation on the raw mesh
    Foam::pointMesh pMesh(mesh);

    // Search for list of objects for this time
    Foam::IOobjectList objects(mesh, dbPtr_().timeName());

    // Convert volume fields
    if (debug)
    {
        Info<< "converting Foam volume fields" << endl;
    }
    Foam::volPointInterpolation pInterp(mesh, pMesh);
    convertVolFields<Foam::scalar>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection(), output
    );
    convertVolFields<Foam::vector>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection(), output
    );
    convertVolFields<Foam::sphericalTensor>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection(), output
    );
    convertVolFields<Foam::symmTensor>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection(), output
    );
    convertVolFields<Foam::tensor>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection(), output
    );
}


void Foam::vtkPV3Foam::updatePointFields
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updatePointFields" << endl;
    }

    const fvMesh& mesh = *meshPtr_;

    // Search for list of objects for this time
    Foam::IOobjectList objects(mesh, dbPtr_().timeName());

    // Convert point fields
    if (debug)
    {
        Info<< "converting Foam point fields" << endl;
    }

    convertPointFields<Foam::scalar>
    (
        mesh, objects, reader_->GetPointFieldSelection(), output
    );
    convertPointFields<Foam::vector>
    (
        mesh, objects, reader_->GetPointFieldSelection(), output
    );
    convertPointFields<Foam::sphericalTensor>
    (
        mesh, objects, reader_->GetPointFieldSelection(), output
    );
    convertPointFields<Foam::symmTensor>
    (
        mesh, objects, reader_->GetPointFieldSelection(), output
    );
    convertPointFields<Foam::tensor>
    (
        mesh, objects, reader_->GetPointFieldSelection(), output
    );
}


void Foam::vtkPV3Foam::updateLagrangianFields
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateLagrangianFields" << endl;
    }

    const fvMesh& mesh = *meshPtr_;

    // Search for list of objects for this time
    //- TODO - currently hard-coded to ONE cloud
    Foam::IOobjectList lagrangianObjects
    (
        mesh,
        dbPtr_().timeName(),
        "lagrangian"/cloudName_
    );

    // Convert Lagrangian fields
    if (debug)
    {
        Info<< "converting Foam Lagrangian fields" << endl;
    }

    convertLagrangianFields<Foam::scalar>
    (
        mesh, lagrangianObjects, reader_->GetLagrangianFieldSelection(), output
    );
    convertLagrangianFields<Foam::vector>
    (
        mesh, lagrangianObjects, reader_->GetLagrangianFieldSelection(), output
    );
    convertLagrangianFields<Foam::sphericalTensor>
    (
        mesh, lagrangianObjects, reader_->GetLagrangianFieldSelection(), output
    );
    convertLagrangianFields<Foam::symmTensor>
    (
        mesh, lagrangianObjects, reader_->GetLagrangianFieldSelection(), output
    );
    convertLagrangianFields<Foam::tensor>
    (
        mesh, lagrangianObjects, reader_->GetLagrangianFieldSelection(), output
    );
}


// ************************************************************************* //
