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

Description

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "cellSet.H"
#include "faceSet.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "pointSet.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "vtkPV3FoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"


// * * * * * * * * * * * * * * * Private Classes * * * * * * * * * * * * * * //

namespace Foam
{

//- A class for reading zone information without requiring a mesh
class zonesEntries
:
    public regIOobject,
    public PtrList<entry>
{

public:

    // Constructors

        explicit zonesEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream("regIOobject"))
        {
            close();
        }

   // Member functions

        bool writeData(Ostream&) const
        {
            notImplemented("zonesEntries::writeData(Ostream&) const");
            return false;
        }
};

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "vtkPV3FoamAddFields.H"
#include "vtkPV3FoamUpdateInformationFields.H"

void Foam::vtkPV3Foam::updateInformationInternalMesh()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationInternalMesh"
            << endl;
    }

    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    // Determine number of meshes available
    HashTable<const fvMesh*> meshObjects = dbPtr_().lookupClass<const fvMesh>();
    nMesh_ = meshObjects.size();

    // Determine regions (internal mesh and patches...)
    //- Add internal mesh as first entry
    selectInfoVolume_ = arraySelection->GetNumberOfArrays();
    arraySelection->AddArray("internalMesh");
    selectInfoVolume_ += 1;
}


void Foam::vtkPV3Foam::updateInformationLagrangian()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationLagrangian "
            << "at timePath " << dbPtr_->timePath()/"lagrangian" << endl;
    }

    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    // Search for list of lagrangian objects for this time
    fileNameList cloudDirs
    (
        readDir(dbPtr_->timePath()/"lagrangian", fileName::DIRECTORY)
    );

    selectInfoLagrangian_ = arraySelection->GetNumberOfArrays();

    if (cloudDirs.size())
    {
        arraySelection->AddArray("lagrangian");
        selectInfoLagrangian_ += 1;

        Info<<"added cloudDirs\n";

        if (cloudDirs.size() > 1)
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
    else
    {
        Info<<"no clouds identified in "
            << dbPtr_->timePath()/"lagrangian" << endl;
    }

}


void Foam::vtkPV3Foam::updateInformationPatches()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationPatches" << endl;
    }

    vtkDataArraySelection *arraySelection = reader_->GetRegionSelection();

    // Read patches
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

    selectInfoPatches_ = arraySelection->GetNumberOfArrays();
    int nPatches = 0;

    // Start regions at patches
    forAll (patchEntries, entryI)
    {
        label nFaces(readLabel(patchEntries[entryI].dict().lookup("nFaces")));

        // Valid patch if nFace > 0
        if (nFaces)
        {
            // Add patch to GUI region list
            arraySelection->AddArray
            (
                (patchEntries[entryI].keyword() + " - patch").c_str()
            );

            ++nPatches;
        }
    }
    selectInfoPatches_ += nPatches;

}


void Foam::vtkPV3Foam::updateInformationZones()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationZones" << endl;
    }

    vtkDataArraySelection *arraySelection = reader_->GetRegionSelection();

    // Read cell zone information
    {
        zonesEntries zones
        (
            IOobject
            (
                "cellZones",
                dbPtr_().findInstance(polyMesh::meshSubDir, "cellZones"),
                polyMesh::meshSubDir,
                dbPtr_(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        selectInfoCellZones_ = arraySelection->GetNumberOfArrays();
        if (zones.headerOk())
        {
            forAll(zones, zoneI)
            {
                arraySelection->AddArray
                (
                    (zones[zoneI].keyword() + " - cellZone").c_str()
                );
            }
            selectInfoCellZones_ += zones.size();
        }

        superCellZonesCells_.setSize(selectInfoCellZones_.size());
    }

    // Read face zone information
    {
        zonesEntries zones
        (
            IOobject
            (
                "faceZones",
                dbPtr_().findInstance(polyMesh::meshSubDir, "faceZones"),
                polyMesh::meshSubDir,
                dbPtr_(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        selectInfoFaceZones_ = arraySelection->GetNumberOfArrays();
        if (zones.headerOk())
        {
            forAll(zones, zoneI)
            {
                arraySelection->AddArray
                (
                    (zones[zoneI].keyword() + " - faceZone").c_str()
                );
            }
            selectInfoFaceZones_ += zones.size();
        }
    }

    // Read point zone information
    {
        zonesEntries zones
        (
            IOobject
            (
                "pointZones",
                dbPtr_().findInstance(polyMesh::meshSubDir, "pointZones"),
                polyMesh::meshSubDir,
                dbPtr_(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        selectInfoPointZones_ = arraySelection->GetNumberOfArrays();
        if (zones.headerOk())
        {
            forAll(zones, zoneI)
            {
                arraySelection->AddArray
                (
                    (zones[zoneI].keyword() + " - pointZone").c_str()
                );
            }
            selectInfoPointZones_ += zones.size();
        }
    }
}


void Foam::vtkPV3Foam::updateInformationSets()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationSets" << endl;
    }

    vtkDataArraySelection *arraySelection = reader_->GetRegionSelection();

    // Add sets
    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().findInstance(polyMesh::meshSubDir, "faces"),
        polyMesh::meshSubDir/"sets"
    );


    selectInfoCellSets_ = arraySelection->GetNumberOfArrays();
    selectInfoCellSets_ += addFields<cellSet>
    (
        arraySelection,
        objects,
        " - cellSet"
    );
    superCellSetCells_.setSize(selectInfoCellSets_.size());

    selectInfoFaceSets_ = arraySelection->GetNumberOfArrays();
    selectInfoFaceSets_ += addFields<faceSet>
    (
        arraySelection,
        objects,
        " - faceSet"
    );

    selectInfoPointSets_ = arraySelection->GetNumberOfArrays();
    selectInfoPointSets_ += addFields<pointSet>
    (
        arraySelection,
        objects,
        " - pointSet"
    );
}


void Foam::vtkPV3Foam::updateInformationLagrangianFields()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateInformationLagrangianFields"
            << endl;
    }

    vtkDataArraySelection *arraySelection =
        reader_->GetLagrangianFieldSelection();

    // preserve the currently selected values
    const stringList selectedEntries = getSelectedArrayEntries
    (
        arraySelection
    );
    arraySelection->RemoveAllArrays();

    // TODO - currently hard-coded to ONE cloud
    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().timeName(),
        "lagrangian"/cloudName_
    );

    addFields<IOField<label> >
    (
        arraySelection,
        objects
    );
    addFields<IOField<scalar> >
    (
        arraySelection,
        objects
    );
    addFields<IOField<vector> >
    (
        arraySelection,
        objects
    );
    addFields<IOField<sphericalTensor> >
    (
        arraySelection,
        objects
    );
    addFields<IOField<symmTensor> >
    (
        arraySelection,
        objects
    );
    addFields<IOField<tensor> >
    (
        arraySelection,
        objects
    );

    // restore the currently enabled values
    setSelectedArrayEntries
    (
        arraySelection,
        selectedEntries
    );

    if (debug)
    {
        Info<< "lagrangian objects.size() = " << objects.size()
            << endl;
    }
}


// ************************************************************************* //
