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
#include "cellSet.H"
#include "faceSet.H"
#include "IOobjectList.H"
#include "pointSet.H"
#include "polyBoundaryMeshEntries.H"
#include "vtkPV3FoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "vtkPV3FoamAddFields.H"
#include "vtkPV3FoamUpdateInformationFields.H"

void Foam::vtkPV3Foam::updateInformationInternalMesh()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationInternalMesh" << endl;
    }

    // Determine number of meshes available
    HashTable<const fvMesh*> meshObjects = dbPtr_().lookupClass<const fvMesh>();
    nMesh_ = meshObjects.size();

    // Determine regions (internal mesh and patches...)
    //- Add internal mesh as first entry
    idRegionVolume_ = reader_->GetRegionSelection()->GetNumberOfArrays();
    reader_->GetRegionSelection()->AddArray("internalMesh");
}


void Foam::vtkPV3Foam::updateInformationLagrangian()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationLagrangian" << endl;
    }

    // Search for list of lagrangian objects for this time
//    IOobjectList lagrangianObjects(dbPtr(), dbPtr_().timeName(), "lagrangian");
//    lagrangianDataSize_ = lagrangianObjects.size();

    fileNameList cloudDirs
    (
        readDir(dbPtr_->timePath()/"lagrangian", fileName::DIRECTORY)
    );
    lagrangianDataSize_ = cloudDirs.size();

    if (lagrangianDataSize_)
    {
        idRegionLagrangian_ = reader_->GetRegionSelection()
            ->GetNumberOfArrays();
        reader_->GetRegionSelection()->AddArray("lagrangian");

        if (lagrangianDataSize_ > 1)
        {
            WarningIn("void Foam::vtkPV3Foam::updateInformationLagrangian()")
                << "Multiple lagrangian clouds identified. Currently only able "
                << "to process ONE cloud: " << cloudDirs[0]
                << endl;
        }
        // Set cloud name to first cloud found
        // TODO - multiple clouds
        cloudName_ = cloudDirs[0];
    }
}


void Foam::vtkPV3Foam::updateInformationPatches()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationPatches" << endl;
    }

    //- Read patches
    polyBoundaryMeshEntries patchEntries
    (
        IOobject
        (
            "boundary",
            dbPtr_().findInstance(polyMesh::meshSubDir, "boundary"),
            polyMesh::meshSubDir,
            dbPtr_(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Start regions at patches
    forAll (patchEntries, entryi)
    {
        label nFaces(readLabel(patchEntries[entryi].dict().lookup("nFaces")));

        //- Valid patch if nFace > 0
        if (nFaces)
        {
            // Add patch to GUI region list
            if (idRegionPatches_ < 0)
            {
                idRegionPatches_ = reader_->GetRegionSelection()
                    ->GetNumberOfArrays();
            }
            reader_->GetRegionSelection()->AddArray
            (
                patchEntries[entryi].keyword().c_str()
            );

            patchDataSize_++;
        }
    }
}


void Foam::vtkPV3Foam::updateInformationSets()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationSets" << endl;
    }

    // Add sets
    IOobjectList setObjects
    (
        dbPtr_(),
        dbPtr_().findInstance(polyMesh::meshSubDir, "faces"),
        polyMesh::meshSubDir/"sets"
    );
    IOobjectList cellSetObjects(setObjects.lookupClass(cellSet::typeName));
    cellSetDataSize_ = cellSetObjects.size();
    superCellSetCells_.setSize(cellSetDataSize_);
    idRegionCellSets_ = reader_->GetRegionSelection()->GetNumberOfArrays();
    addFields<cellSet>
    (
        reader_->GetRegionSelection(),
        cellSetObjects
    );

    IOobjectList faceSetObjects(setObjects.lookupClass(faceSet::typeName));
    faceSetDataSize_ = faceSetObjects.size();
    idRegionFaceSets_ = reader_->GetRegionSelection()->GetNumberOfArrays();
    addFields<faceSet>
    (
        reader_->GetRegionSelection(),
        faceSetObjects
    );

    IOobjectList pointSetObjects(setObjects.lookupClass(pointSet::typeName));
    pointSetDataSize_ = pointSetObjects.size();
    idRegionPointSets_ = reader_->GetRegionSelection()->GetNumberOfArrays();
    addFields<pointSet>
    (
        reader_->GetRegionSelection(),
        pointSetObjects
    );
}


void Foam::vtkPV3Foam::updateInformationLagrangianFields()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationLagrangianFields"
            << endl;
    }

    const wordList selectedArrayEntries = getSelectedArrayEntries
    (
        reader_->GetLagrangianFieldSelection()
    );
    reader_->GetLagrangianFieldSelection()->RemoveAllArrays();

    // TODO - currently hard-coded to ONE cloud
    IOobjectList lagrangianObjects(dbPtr_(), dbPtr_().timeName(), "lagrangian"/cloudName_);

    addFields<IOField<scalar> >
    (
        reader_->GetLagrangianFieldSelection(),
        lagrangianObjects
    );
    addFields<IOField<vector> >
    (
        reader_->GetLagrangianFieldSelection(),
        lagrangianObjects
    );
    addFields<IOField<sphericalTensor> >
    (
        reader_->GetLagrangianFieldSelection(),
        lagrangianObjects
    );
    addFields<IOField<symmTensor> >
    (
        reader_->GetLagrangianFieldSelection(),
        lagrangianObjects
    );
    addFields<IOField<tensor> >
    (
        reader_->GetLagrangianFieldSelection(),
        lagrangianObjects
    );

    setSelectedArrayEntries
    (
        reader_->GetLagrangianFieldSelection(),
        selectedArrayEntries
    );

    if (debug)
    {
        Info<< "lagrangianObjects.size() = " << lagrangianObjects.size()
            << endl;
    }
}


// ************************************************************************* //
