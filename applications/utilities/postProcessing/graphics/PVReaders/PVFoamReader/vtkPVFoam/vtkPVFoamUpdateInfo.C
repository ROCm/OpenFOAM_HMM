/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "Cloud.H"
#include "vtkPVFoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"

// Templates (only needed here)
#include "vtkPVFoamUpdateTemplates.C"


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
            NotImplemented;
            return false;
        }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType>
Foam::wordList Foam::vtkPVFoam::getZoneNames
(
    const ZoneMesh<ZoneType, polyMesh>& zmesh
)
{
    wordList names(zmesh.size());
    label nZone = 0;

    forAll(zmesh, zonei)
    {
        if (!zmesh[zonei].empty())
        {
            names[nZone++] = zmesh[zonei].name();
        }
    }
    names.setSize(nZone);

    return names;
}


Foam::wordList Foam::vtkPVFoam::getZoneNames(const word& zoneType) const
{
    wordList names;

    // mesh not loaded - read from file
    IOobject ioObj
    (
        zoneType,
        dbPtr_().findInstance
        (
            meshDir_,
            zoneType,
            IOobject::READ_IF_PRESENT
        ),
        meshDir_,
        dbPtr_(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (ioObj.typeHeaderOk<cellZoneMesh>(false, false))
    {
        zonesEntries zones(ioObj);

        names.setSize(zones.size());
        forAll(zones, zonei)
        {
            names[zonei] = zones[zonei].keyword();
        }
    }

    return names;
}


void Foam::vtkPVFoam::updateInfoInternalMesh
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        Info<< "<beg> updateInfoInternalMesh" << endl;
    }

    // Determine mesh parts (internalMesh, patches...)
    //- Add internal mesh as first entry
    rangeVolume_.reset(arraySelection->GetNumberOfArrays(), 1);
    arraySelection->AddArray("internalMesh");

    if (debug)
    {
        Info<< "<end> updateInfoInternalMesh" << endl;
    }
}


void Foam::vtkPVFoam::updateInfoLagrangian
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        Info<< "<beg> updateInfoLagrangian" << nl
            << "    " << dbPtr_->timePath()/cloud::prefix << endl;
    }


    // use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/cloud::prefix;
    }

    // Search for list of lagrangian objects for this time
    fileNameList cloudDirs
    (
        readDir(dbPtr_->timePath()/lagrangianPrefix, fileName::DIRECTORY)
    );

    rangeLagrangian_.reset(arraySelection->GetNumberOfArrays());
    forAll(cloudDirs, cloudi)
    {
        // Add cloud to GUI list
        arraySelection->AddArray
        (
            ("lagrangian/" + cloudDirs[cloudi]).c_str()
        );
        ++rangeLagrangian_;
    }

    if (debug)
    {
        Info<< "<end> updateInfoLagrangian" << endl;
    }
}


void Foam::vtkPVFoam::updateInfoPatches
(
    vtkDataArraySelection* arraySelection,
    stringList& enabledEntries
)
{
    if (debug)
    {
        Info<< "<beg> updateInfoPatches"
            << " [meshPtr=" << (meshPtr_ ? "set" : "null") << "]" << endl;
    }

    HashSet<string> enabledEntriesSet(enabledEntries);

    rangePatches_.reset(arraySelection->GetNumberOfArrays());

    if (meshPtr_)
    {
        const polyBoundaryMesh& patches = meshPtr_->boundaryMesh();
        const HashTable<labelList>& groups = patches.groupPatchIDs();

        // Add (non-zero) patch groups to the list of mesh parts
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAllConstIters(groups, iter)
        {
            const auto& groupName = iter.key();
            const auto& patchIDs  = iter.object();

            label nFaces = 0;
            for (auto patchId : patchIDs)
            {
                nFaces += patches[patchId].size();
            }

            if (!nFaces)
            {
                // Skip if group has no faces
                continue;
            }

            // Valid patch if nFace > 0 - add patch to GUI list
            const string dpyName = "group/" + groupName;
            arraySelection->AddArray(dpyName.c_str());
            ++rangePatches_;

            if (enabledEntriesSet.found(dpyName))
            {
                if (!reader_->GetShowGroupsOnly())
                {
                    enabledEntriesSet.erase(dpyName);
                    for (auto patchId : patchIDs)
                    {
                        const polyPatch& pp = patches[patchId];
                        if (pp.size())
                        {
                            enabledEntriesSet.insert
                            (
                                "patch/" + pp.name()
                            );
                        }
                    }
                }
            }
        }


        // Add (non-zero) patches to the list of mesh parts
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (!reader_->GetShowGroupsOnly())
        {
            forAll(patches, patchi)
            {
                const polyPatch& pp = patches[patchi];

                if (pp.size())
                {
                    // Add patch to GUI list
                    arraySelection->AddArray
                    (
                        ("patch/" + pp.name()).c_str()
                    );
                    ++rangePatches_;
                }
            }
        }
    }
    else
    {
        // mesh not loaded - read from file
        // but this could fail if we've supplied a bad region name
        IOobject ioObj
        (
            "boundary",
            dbPtr_().findInstance
            (
                meshDir_,
                "boundary",
                IOobject::READ_IF_PRESENT
            ),
            meshDir_,
            dbPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        // this should only ever fail if the mesh region doesn't exist
        if (ioObj.typeHeaderOk<polyBoundaryMesh>(true, false))
        {
            polyBoundaryMeshEntries patchEntries(ioObj);

            // Read patches and determine sizes
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            wordList names(patchEntries.size());
            labelList sizes(patchEntries.size());

            forAll(patchEntries, patchi)
            {
                const dictionary& patchDict = patchEntries[patchi].dict();

                sizes[patchi] = readLabel(patchDict.lookup("nFaces"));
                names[patchi] = patchEntries[patchi].keyword();
            }


            // Add (non-zero) patch groups to the list of mesh parts
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            HashTable<labelHashSet> groups(patchEntries.size());

            forAll(patchEntries, patchi)
            {
                const dictionary& patchDict = patchEntries[patchi].dict();
                wordList groupNames;

                if
                (
                    sizes[patchi]  // Valid patch if nFace > 0
                 && patchDict.readIfPresent("inGroups", groupNames)
                )
                {
                    forAll(groupNames, groupI)
                    {
                        groups(groupNames[groupI]).insert(patchi);
                    }
                }
            }

            forAllConstIters(groups, iter)
            {
                const auto& groupName = iter.key();
                const auto& patchIDs  = iter.object();

                const string dpyName = "group/" + groupName;
                arraySelection->AddArray(dpyName.c_str());
                ++rangePatches_;

                // Optionally replace with patch name selection
                if
                (
                    enabledEntriesSet.found(dpyName)
                 && !reader_->GetShowGroupsOnly()
                )
                {
                    enabledEntriesSet.erase(dpyName);
                    for (auto patchId : patchIDs)
                    {
                        enabledEntriesSet.insert
                        (
                            "patch/" + names[patchId]
                        );
                    }
                }
            }


            // Add (non-zero) patches to the list of mesh parts
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (!reader_->GetShowGroupsOnly())
            {
                forAll(names, patchi)
                {
                    // Valid patch if nFace > 0 - add patch to GUI list
                    if (sizes[patchi])
                    {
                        arraySelection->AddArray
                        (
                            ("patch/" + names[patchi]).c_str()
                        );
                        ++rangePatches_;
                    }
                }
            }
        }
    }

    // Update enabled entries in case of group selection
    enabledEntries = enabledEntriesSet.toc();

    if (debug)
    {
        Info<< "<end> updateInfoPatches" << endl;
    }
}


void Foam::vtkPVFoam::updateInfoZones
(
    vtkDataArraySelection* arraySelection
)
{
    if (!reader_->GetIncludeZones())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> updateInfoZones"
            << " [meshPtr=" << (meshPtr_ ? "set" : "null") << "]" << endl;
    }

    wordList namesLst;

    //
    // cellZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->cellZones());
    }
    else
    {
        namesLst = getZoneNames("cellZones");
    }

    rangeCellZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            ("cellZone/" + namesLst[elemI]).c_str()
        );
        ++rangeCellZones_;
    }


    //
    // faceZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->faceZones());
    }
    else
    {
        namesLst = getZoneNames("faceZones");
    }

    rangeFaceZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            ("faceZone/" + namesLst[elemI]).c_str()
        );
        ++rangeFaceZones_;
    }


    //
    // pointZones information
    // ~~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->pointZones());
    }
    else
    {
        namesLst = getZoneNames("pointZones");
    }

    rangePointZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            ("pointZone/" + namesLst[elemI]).c_str()
        );
        ++rangePointZones_;
    }

    if (debug)
    {
        Info<< "<end> updateInfoZones" << endl;
    }
}


void Foam::vtkPVFoam::updateInfoSets
(
    vtkDataArraySelection* arraySelection
)
{
    if (!reader_->GetIncludeSets())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> updateInfoSets" << endl;
    }

    // Add names of sets. Search for last time directory with a sets
    // subdirectory. Take care not to search beyond the last mesh.

    word facesInstance = dbPtr_().findInstance
    (
        meshDir_,
        "faces",
        IOobject::READ_IF_PRESENT
    );

    word setsInstance = dbPtr_().findInstance
    (
        meshDir_/"sets",
        word::null,
        IOobject::READ_IF_PRESENT,
        facesInstance
    );

    IOobjectList objects(dbPtr_(), setsInstance, meshDir_/"sets");

    if (debug)
    {
        Info<< "     updateInfoSets read "
            << objects.names() << " from " << setsInstance << endl;
    }


    rangeCellSets_.reset(arraySelection->GetNumberOfArrays());
    rangeCellSets_ += addToSelection<cellSet>
    (
        arraySelection,
        objects,
        "cellSet/"
    );

    rangeFaceSets_.reset(arraySelection->GetNumberOfArrays());
    rangeFaceSets_ += addToSelection<faceSet>
    (
        arraySelection,
        objects,
        "faceSet/"
    );

    rangePointSets_.reset(arraySelection->GetNumberOfArrays());
    rangePointSets_ += addToSelection<pointSet>
    (
        arraySelection,
        objects,
        "pointSet/"
    );

    if (debug)
    {
        Info<< "<end> updateInfoSets" << endl;
    }
}


void Foam::vtkPVFoam::updateInfoFields()
{
    updateInfoFields<fvPatchField, volMesh>
    (
        reader_->GetVolFieldSelection()
    );
    updateInfoFields<pointPatchField, pointMesh>
    (
        reader_->GetPointFieldSelection()
    );
    updateInfoLagrangianFields
    (
        reader_->GetLagrangianFieldSelection()
    );
}


void Foam::vtkPVFoam::updateInfoLagrangianFields
(
    vtkDataArraySelection* select
)
{
    if (debug)
    {
        Info<< "<beg> updateInfoLagrangianFields" << endl;
    }

    // preserve the enabled selections
    stringList enabledEntries = getSelectedArrayEntries(select);
    select->RemoveAllArrays();

    // TODO - currently only get fields from ONE cloud
    // have to decide if the second set of fields get mixed in
    // or dealt with separately

    const arrayRange& range = rangeLagrangian_;
    if (range.empty())
    {
        return;
    }

    int partId = range.start();
    word cloudName = getPartName(partId);

    // use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/cloud::prefix;
    }

    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().timeName(),
        lagrangianPrefix/cloudName
    );

    addToSelection<IOField<label>>(select, objects);
    addToSelection<IOField<scalar>>(select, objects);
    addToSelection<IOField<vector>>(select, objects);
    addToSelection<IOField<sphericalTensor>>(select, objects);
    addToSelection<IOField<symmTensor>>(select, objects);
    addToSelection<IOField<tensor>>(select, objects);

    // restore the enabled selections
    setSelectedArrayEntries(select, enabledEntries);

    if (debug > 1)
    {
        boolList status;
        const label nElem = getSelected(status, select);

        forAll(status, i)
        {
            Info<< "  lagrangian[" << i << "] = "
                << status[i]
                << " : " << select->GetArrayName(i) << nl;
        }

        if (!nElem)
        {
            Info<< "  lagrangian[none]" << nl;
        }
    }

    if (debug)
    {
        Info<< "<end> updateInfoLagrangianFields - "
            << "lagrangian objects.size() = " << objects.size() << endl;
    }
}


// ************************************************************************* //
