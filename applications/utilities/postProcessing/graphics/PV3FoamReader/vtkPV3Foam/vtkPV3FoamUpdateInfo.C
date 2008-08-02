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

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "vtkPV3FoamReader.h"

// local headers
#include "vtkPV3FoamAddToSelection.H"
#include "vtkPV3FoamUpdateInfoFields.H"

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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::vtkPV3Foam::readZoneNames(const word& zoneType)
{
    wordList zoneNames;

    // mesh not loaded - read from file
    IOobject ioObj
    (
        zoneType,
        dbPtr_().findInstance
        (
            polyMesh::meshSubDir,
            zoneType,
            IOobject::READ_IF_PRESENT
        ),
        polyMesh::meshSubDir,
        dbPtr_(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (ioObj.headerOk())
    {
        zonesEntries zones(ioObj);

        zoneNames.setSize(zones.size());
        forAll(zones, zoneI)
        {
            zoneNames[zoneI] = zones[zoneI].keyword();
        }
    }

    return zoneNames;
}


void Foam::vtkPV3Foam::updateInfoInternalMesh()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoInternalMesh" << endl;
    }

    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // Determine number of meshes available
    HashTable<const fvMesh*> meshObjects = dbPtr_().lookupClass<const fvMesh>();
    nMesh_ = meshObjects.size();

    // Determine regions (internal mesh and patches...)
    //- Add internal mesh as first entry
    regionInfoVolume_ = regionSelection->GetNumberOfArrays();
    regionSelection->AddArray("internalMesh");
    regionInfoVolume_ += 1;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(regionSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoInternalMesh" << endl;
    }

}


void Foam::vtkPV3Foam::updateInfoLagrangian()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoLagrangian" << nl
            << "    " << dbPtr_->timePath()/"lagrangian" << endl;
    }

    // Search for list of lagrangian objects for this time
    fileNameList cloudDirs
    (
        readDir(dbPtr_->timePath()/"lagrangian", fileName::DIRECTORY)
    );

    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();
    regionInfoLagrangian_ = regionSelection->GetNumberOfArrays();

    int nClouds = 0;
    forAll(cloudDirs, cloudI)
    {
        // Add cloud to GUI list
        regionSelection->AddArray
        (
            (cloudDirs[cloudI] + " - lagrangian").c_str()
        );

        ++nClouds;
    }

    regionInfoLagrangian_ += nClouds;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(regionSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoLagrangian" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoPatches()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoPatches"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }

    vtkDataArraySelection *regionSelection = reader_->GetRegionSelection();
    regionInfoPatches_ = regionSelection->GetNumberOfArrays();

    int nPatches = 0;
    if (meshPtr_)
    {
        const polyBoundaryMesh& patches = meshPtr_->boundaryMesh();
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.size())
            {
                // Add patch to GUI list
                regionSelection->AddArray
                (
                    (pp.name() + " - patch").c_str()
                );

                ++nPatches;
            }
        }
    }
    else
    {
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

        // Start regions at patches
        forAll(patchEntries, entryI)
        {
            label nFaces
            (
                readLabel(patchEntries[entryI].dict().lookup("nFaces"))
            );

            // Valid patch if nFace > 0
            if (nFaces)
            {
                // Add patch to GUI list
                regionSelection->AddArray
                (
                    (patchEntries[entryI].keyword() + " - patch").c_str()
                );

                ++nPatches;
            }
        }
    }

    regionInfoPatches_ += nPatches;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(regionSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoPatches" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoZones()
{
    if (!reader_->GetIncludeZones())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoZones"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }

    vtkDataArraySelection *regionSelection = reader_->GetRegionSelection();

    wordList namesLst;

    //
    // cellZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->cellZones().names();
    }
    else
    {
        namesLst = readZoneNames("cellZones");
    }

    regionInfoCellZones_ = regionSelection->GetNumberOfArrays();
    forAll(namesLst, elemI)
    {
        regionSelection->AddArray((namesLst[elemI] + " - cellZone").c_str());
    }
    regionInfoCellZones_ += namesLst.size();
    zoneSuperCells_.setSize(regionInfoCellZones_.size());


    //
    // faceZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->faceZones().names();
    }
    else
    {
        namesLst = readZoneNames("faceZones");
    }

    regionInfoFaceZones_ = regionSelection->GetNumberOfArrays();
    forAll(namesLst, elemI)
    {
        regionSelection->AddArray
        (
            (namesLst[elemI] + " - faceZone").c_str()
        );
    }
    regionInfoFaceZones_ += namesLst.size();


    //
    // pointZones information
    // ~~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->pointZones().names();
    }
    else
    {
        namesLst = readZoneNames("pointZones");
    }

    regionInfoPointZones_ = regionSelection->GetNumberOfArrays();
    forAll(namesLst, elemI)
    {
        regionSelection->AddArray
        (
            (namesLst[elemI] + " - pointZone").c_str()
        );
    }
    regionInfoPointZones_ += namesLst.size();


    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(regionSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoZones" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoSets()
{
    if (!reader_->GetIncludeSets())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoSets" << endl;
    }

    vtkDataArraySelection *regionSelection = reader_->GetRegionSelection();

    // Add names of sets
    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().findInstance(polyMesh::meshSubDir, "faces"),
        polyMesh::meshSubDir/"sets"
    );


    regionInfoCellSets_ = regionSelection->GetNumberOfArrays();
    regionInfoCellSets_ += addToSelection<cellSet>
    (
        regionSelection,
        objects,
        " - cellSet"
    );
    csetSuperCells_.setSize(regionInfoCellSets_.size());

    regionInfoFaceSets_ = regionSelection->GetNumberOfArrays();
    regionInfoFaceSets_ += addToSelection<faceSet>
    (
        regionSelection,
        objects,
        " - faceSet"
    );

    regionInfoPointSets_ = regionSelection->GetNumberOfArrays();
    regionInfoPointSets_ += addToSelection<pointSet>
    (
        regionSelection,
        objects,
        " - pointSet"
    );

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(regionSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoSets" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoLagrangianFields()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoLagrangianFields"
            << endl;
    }

    vtkDataArraySelection *fieldSelection =
        reader_->GetLagrangianFieldSelection();

    vtkDataArraySelection *regionSelection =
        reader_->GetRegionSelection();

    // preserve the enabled selections
    stringList selectedEntries = getSelectedArrayEntries
    (
        fieldSelection,
        true
    );

    fieldSelection->RemoveAllArrays();


    //
    // TODO - can currently only get fields from ONE cloud
    //

    const selectionInfo& selector = regionInfoLagrangian_;
    int regionId = selector.start();

    if (!selector.size() || regionId < 0)
    {
        return;
    }

    word cloudName = getFirstWord
    (
        regionSelection->GetArrayName(regionId)
    );

    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().timeName(),
        "lagrangian"/cloudName
    );

    addToSelection<IOField<label> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<scalar> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<vector> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<sphericalTensor> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<symmTensor> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<tensor> >
    (
        fieldSelection,
        objects
    );

    // restore the enabled selections
    setSelectedArrayEntries
    (
        fieldSelection,
        selectedEntries
    );

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateInfoLagrangianFields - "
            << "lagrangian objects.size() = " << objects.size() << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
