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
    vtkDataArraySelection* select
)
{
    if (debug)
    {
        Info<< "<beg> updateInfoInternalMesh" << endl;
    }

    // Determine mesh parts (internalMesh, patches...)
    //- Add internal mesh as first entry
    rangeVolume_.reset(select->GetNumberOfArrays(), 1);
    select->AddArray("internalMesh");

    if (debug)
    {
        Info<< "<end> updateInfoInternalMesh" << endl;
    }
}


void Foam::vtkPVFoam::updateInfoLagrangian
(
    vtkDataArraySelection* select
)
{
    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl
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

    rangeLagrangian_.reset(select->GetNumberOfArrays());
    forAll(cloudDirs, cloudi)
    {
        // Add cloud to GUI list
        select->AddArray
        (
            ("lagrangian/" + cloudDirs[cloudi]).c_str()
        );
        ++rangeLagrangian_;
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
    }
}


void Foam::vtkPVFoam::updateInfoPatches
(
    vtkDataArraySelection* select,
    HashSet<string>& enabledEntries
)
{
    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME
            << " [meshPtr=" << (meshPtr_ ? "set" : "null") << "]" << endl;
    }

    rangePatches_.reset(select->GetNumberOfArrays());

    if (meshPtr_)
    {
        const polyBoundaryMesh& patches = meshPtr_->boundaryMesh();
        const HashTable<labelList>& groups = patches.groupPatchIDs();
        DynamicList<string> displayNames(groups.size());

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
            displayNames.append(dpyName);

            // Optionally replace group with patch name selections
            // - must remove the group from the select itself, otherwise
            //   it can toggle on, but not toggle off very well
            if
            (
                !reader_->GetShowGroupsOnly()
             && enabledEntries.erase(dpyName)
            )
            {
                for (auto patchId : patchIDs)
                {
                    const polyPatch& pp = patches[patchId];
                    if (pp.size())
                    {
                        enabledEntries.insert
                        (
                            "patch/" + pp.name()
                        );
                    }
                }
            }
        }

        // Sort group names
        Foam::sort(displayNames);
        for (const auto& name : displayNames)
        {
            select->AddArray(name.c_str());
            ++rangePatches_;
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
                    select->AddArray
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

        // This should only ever fail if the mesh region doesn't exist
        if (ioObj.typeHeaderOk<polyBoundaryMesh>(true, false))
        {
            polyBoundaryMeshEntries patchEntries(ioObj);

            // Read patches, determine sizes and patch groups
            wordList names(patchEntries.size());
            labelList sizes(patchEntries.size());
            HashTable<labelHashSet> groups(2*patchEntries.size());

            forAll(patchEntries, patchi)
            {
                const dictionary& patchDict = patchEntries[patchi].dict();
                wordList groupNames;

                sizes[patchi] = readLabel(patchDict.lookup("nFaces"));
                names[patchi] = patchEntries[patchi].keyword();

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


            // Add (non-zero) patch groups to the list of mesh parts
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            DynamicList<string> displayNames(groups.size());

            forAllConstIters(groups, iter)
            {
                const auto& groupName = iter.key();
                const auto& patchIDs  = iter.object();

                const string dpyName = "group/" + groupName;
                displayNames.append(dpyName);

                // Optionally replace group with patch name selections
                // - must remove the group from the select itself, otherwise
                //   it can toggle on, but not toggle off very well
                if
                (
                    !reader_->GetShowGroupsOnly()
                 && enabledEntries.erase(dpyName)
                )
                {
                    for (auto patchId : patchIDs)
                    {
                        enabledEntries.insert
                        (
                            "patch/" + names[patchId]
                        );
                    }
                }
            }

            // Sort group names
            Foam::sort(displayNames);
            for (const auto& name : displayNames)
            {
                select->AddArray(name.c_str());
                ++rangePatches_;
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
                        select->AddArray
                        (
                            ("patch/" + names[patchi]).c_str()
                        );
                        ++rangePatches_;
                    }
                }
            }
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
    }
}


void Foam::vtkPVFoam::updateInfoZones
(
    vtkDataArraySelection* select
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

    rangeCellZones_.reset(select->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        select->AddArray
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

    rangeFaceZones_.reset(select->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        select->AddArray
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

    rangePointZones_.reset(select->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        select->AddArray
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
    vtkDataArraySelection* select
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


    rangeCellSets_.reset(select->GetNumberOfArrays());
    rangeCellSets_ += addToSelection<cellSet>
    (
        select,
        objects,
        "cellSet/"
    );

    rangeFaceSets_.reset(select->GetNumberOfArrays());
    rangeFaceSets_ += addToSelection<faceSet>
    (
        select,
        objects,
        "faceSet/"
    );

    rangePointSets_.reset(select->GetNumberOfArrays());
    rangePointSets_ += addToSelection<pointSet>
    (
        select,
        objects,
        "pointSet/"
    );

    if (debug)
    {
        Info<< "<end> updateInfoSets" << endl;
    }
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

    // Preserve the enabled selections
    HashSet<string> enabled = getSelectedArraySet(select);
    select->RemoveAllArrays();

    // TODO - currently only get fields from ONE cloud
    // have to decide if the second set of fields get mixed in
    // or dealt with separately

    const arrayRange& range = rangeLagrangian_;
    if (range.empty())
    {
        return;
    }

    // Add Lagrangian fields even if particles are not enabled?
    const int partId = range.start();
    const word cloudName = getReaderPartName(partId);

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

    // Restore the enabled selections
    setSelectedArrayEntries(select, enabled);

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
    }
}


// ************************************************************************* //
