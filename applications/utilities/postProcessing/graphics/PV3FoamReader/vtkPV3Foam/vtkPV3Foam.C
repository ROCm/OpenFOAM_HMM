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
#include "vtkPV3FoamReader.h"

// Foam includes
#include "fvMesh.H"
#include "Time.H"
#include "patchZones.H"
#include "IFstream.H"

// VTK includes
#include "vtkCharArray.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkInformation.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::vtkPV3Foam, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "vtkPV3FoamAddToSelection.H"
#include "vtkPV3FoamUpdateInfoFields.H"

void Foam::vtkPV3Foam::AddToBlock
(
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo,
    vtkDataSet* dataset,
    const string& blockName
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);
    if (blockDO && !block)
    {
        FatalErrorIn("Foam::vtkPV3Foam::AddToBlock")
            << "Block already has a vtkDataSet assigned to it" << nl << endl;
        return;
    }

    if (!block)
    {
        block = vtkMultiBlockDataSet::New();
        output->SetBlock(blockNo, block);
        block->Delete();
    }

    if (block)
    {
        if (debug)
        {
            Info<< "block[" << blockNo << "] has "
                << block->GetNumberOfBlocks()
                <<  " datasets prior to adding set " << datasetNo
                <<  " with name: " << blockName << endl;
        }

        // when assigning dataset 0, also name the parent block
        if (!datasetNo && selector.name())
        {
            output->GetMetaData(blockNo)->Set
            (
                vtkCompositeDataSet::NAME(),
                selector.name()
            );
        }
    }


    block->SetBlock(datasetNo, dataset);

    if (blockName.size())
    {
        block->GetMetaData(datasetNo)->Set
        (
            vtkCompositeDataSet::NAME(), blockName.c_str()
        );
    }

}


vtkDataSet* Foam::vtkPV3Foam::GetDataSetFromBlock
(
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);
    if (block)
    {
        return vtkDataSet::SafeDownCast(block->GetBlock(datasetNo));
    }

    return 0;
}


Foam::label Foam::vtkPV3Foam::GetNumberOfDataSets
(
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);
    if (block)
    {
        return block->GetNumberOfBlocks();
    }

    return 0;
}


void Foam::vtkPV3Foam::resetCounters()
{
    // Reset region ids and sizes
    regionInfoVolume_.reset();
    regionInfoPatches_.reset();
    regionInfoLagrangian_.reset();
    regionInfoCellZones_.reset();
    regionInfoFaceZones_.reset();
    regionInfoPointZones_.reset();
    regionInfoCellSets_.reset();
    regionInfoFaceSets_.reset();
    regionInfoPointSets_.reset();
}


int Foam::vtkPV3Foam::setTime(const double& requestedTime)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::setTime(" << requestedTime << ")"
            << endl;
    }

    Time& runTime = dbPtr_();

    // Get times list
    instantList Times = runTime.times();

    int foundIndex = Time::findClosestTimeIndex(Times, requestedTime);
    int nearestIndex = foundIndex;

    if (foundIndex < 0)
    {
        nearestIndex = 0;
    }

    // see what has changed
    if (timeIndex_ != nearestIndex)
    {
        timeIndex_ = nearestIndex;
        runTime.setTime(Times[nearestIndex], nearestIndex);

        // the fields change each time
        fieldsChanged_ = true;

        if (meshPtr_)
        {
            if (meshPtr_->readUpdate() != polyMesh::UNCHANGED)
            {
                meshChanged_ = true;
            }
        }
        else
        {
            meshChanged_ = true;
        }

        reader_->UpdateProgress(0.05);

        // this seems to be needed for catching Lagrangian fields
        updateInfo();
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::setTime() - selected time "
            << Times[nearestIndex].name() << " index=" << nearestIndex
            << " meshChanged=" << meshChanged_
            << " fieldsChanged=" << fieldsChanged_ << endl;
    }

    return foundIndex;
}


void Foam::vtkPV3Foam::updateRegionStatus()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateRegionStatus" << endl;
    }

    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();
    const label nSelect = regionSelection->GetNumberOfArrays();
    if (regionStatus_.size() != nSelect)
    {
        regionStatus_.setSize(nSelect);
        regionStatus_ = false;
        meshChanged_ = true;
    }

    // this needs fixing if we wish to re-use the datasets
    regionDataset_.setSize(nSelect);
    regionDataset_ = -1;

    // Read the selected cell regions, zones, patches and add to region list
    forAll(regionStatus_, regionId)
    {
        int setting = regionSelection->GetArraySetting(regionId);

        if (regionStatus_[regionId] != setting)
        {
            regionStatus_[regionId] = setting;
            meshChanged_ = true;
        }

        if (debug)
        {
            Info<< "  region[" << regionId << "] = "
                << regionStatus_[regionId]
                << " : " << regionSelection->GetArrayName(regionId) << endl;
        }
    }
    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateRegionStatus" << endl;
    }
}

Foam::wordHashSet Foam::vtkPV3Foam::getSelected
(
    vtkDataArraySelection* select
)
{
    int nElem = select->GetNumberOfArrays();
    wordHashSet selections(2*nElem);

    for (int elemI=0; elemI < nElem; ++elemI)
    {
        if (select->GetArraySetting(elemI))
        {
            word name(getFirstWord(select->GetArrayName(elemI)));
            selections.insert(name);
        }
    }

    return selections;
}


Foam::stringList Foam::vtkPV3Foam::getSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const bool firstWord
)
{
    stringList selections(select->GetNumberOfArrays());
    label nElem = 0;

    if (debug)
    {
        Info<< "available(";
        forAll(selections, elemI)
        {
            Info<< " \"" << select->GetArrayName(elemI) << "\"";
        }
        Info<< " )\n"
            << "selected(";
    }

    forAll(selections, elemI)
    {
        if (select->GetArraySetting(elemI))
        {
            if (firstWord)
            {
                selections[nElem] = getFirstWord
                (
                    select->GetArrayName(elemI)
                );
            }
            else
            {
                selections[nElem] = select->GetArrayName(elemI);
            }

            if (debug)
            {
                Info<< " " << selections[nElem];
            }

            ++nElem;
        }
    }

    if (debug)
    {
        Info<< " )" << endl;
    }

    selections.setSize(nElem);
    return selections;
}


Foam::stringList Foam::vtkPV3Foam::getSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const selectionInfo& selector,
    const bool firstWord
)
{
    stringList selections(selector.size());
    label nElem = 0;

    if (debug)
    {
        Info<< "available(";
        for
        (
            int elemI = selector.start();
            elemI < selector.end();
            ++elemI
        )
        {
            Info<< " \"" << select->GetArrayName(elemI) << "\"";
        }

        Info<< " )\n"
            << "selected(";
    }

    for
    (
        int elemI = selector.start();
        elemI < selector.end();
        ++elemI
    )
    {
        if (select->GetArraySetting(elemI))
        {
            if (firstWord)
            {
                selections[nElem] = getFirstWord
                (
                    select->GetArrayName(elemI)
                );
            }
            else
            {
                selections[nElem] = select->GetArrayName(elemI);
            }

            if (debug)
            {
                Info<< " " << selections[nElem];
            }

            ++nElem;
        }
    }

    if (debug)
    {
        Info<< " )" << endl;
    }

    selections.setSize(nElem);
    return selections;
}


void Foam::vtkPV3Foam::setSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const stringList& selections
)
{
    if (debug > 1)
    {
        Info<< "<beg> Foam::vtkPV3Foam::setSelectedArrayEntries" << endl;
    }
    const label nEntries = select->GetNumberOfArrays();

    // Reset all current entries to 'not selected'
    select->DisableAllArrays();

    // Loop through entries, setting values from selectedEntries
    forAll(selections, elemI)
    {
        if (debug > 1)
        {
            Info<< "selections[" << elemI << "] = " << selections[elemI]
                << endl;
        }

        for (label i=0; i<nEntries; i++)
        {
            string arrayName = select->GetArrayName(i);

            if (arrayName == selections[elemI])
            {
                if (debug > 1)
                {
                    Info<< "enabling array: " << arrayName << " Index = "
                        << i
                        << endl;
                }

                select->EnableArray(arrayName.c_str());
                break;
            }
        }
    }
    if (debug > 1)
    {
        Info<< "<end> Foam::vtkPV3Foam::setSelectedArrayEntries" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPV3Foam::vtkPV3Foam
(
    const char* const FileName,
    vtkPV3FoamReader* reader
)
:
    reader_(reader),
    dbPtr_(NULL),
    meshPtr_(NULL),
    nMesh_(0),
    timeIndex_(-1),
    meshChanged_(true),
    fieldsChanged_(true),
    regionInfoVolume_("unzoned"),
    regionInfoPatches_("patches"),
    regionInfoLagrangian_("lagrangian"),
    regionInfoCellZones_("cellZone"),
    regionInfoFaceZones_("faceZone"),
    regionInfoPointZones_("pointZone"),
    regionInfoCellSets_("cellSet"),
    regionInfoFaceSets_("faceSet"),
    regionInfoPointSets_("pointSet")
{
    if (debug)
    {
        Info<< "Foam::vtkPV3Foam::vtkPV3Foam - " << FileName << endl;
        printMemory();
    }

    // avoid argList and get rootPath/caseName directly from the file
    fileName fullCasePath(fileName(FileName).path());

    if (!dir(fullCasePath))
    {
        return;
    }
    if (fullCasePath == ".")
    {
        fullCasePath = cwd();
    }

    // Set the case as an environment variable - some BCs might use this
    if (fullCasePath.name().find("processor", 0) == 0)
    {
        setEnv("FOAM_CASE", fullCasePath.path(), true);
    }
    else
    {
        setEnv("FOAM_CASE", fullCasePath, true);
    }

    if (debug)
    {
        Info<< "fullCasePath=" << fullCasePath << nl
            << "FOAM_CASE=" << getEnv("FOAM_CASE") << endl;
    }

    // Create time object
    dbPtr_.reset
    (
        new Time
        (
            Time::controlDictName,
            fileName(fullCasePath.path()),
            fileName(fullCasePath.name())
        )
    );

    dbPtr_().functionObjects().off();

    updateInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkPV3Foam::~vtkPV3Foam()
{
    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::~vtkPV3Foam" << endl;
    }

    delete meshPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::updateInfo()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfo"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "] timeIndex="
            << timeIndex_ << endl;
    }

    resetCounters();

    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    stringList selectedEntries;
    // enable 'internalMesh' on the first call
    if (regionSelection->GetNumberOfArrays() == 0 && !meshPtr_)
    {
        selectedEntries.setSize(1);
        selectedEntries[0] = "internalMesh";
    }
    else
    {
        // preserve the enabled selections
        selectedEntries = getSelectedArrayEntries
        (
            regionSelection,
            false
        );
    }

    // Clear current region list/array
    regionSelection->RemoveAllArrays();

    // Update region array - add Lagrangian at the bottom
    updateInfoInternalMesh();
    updateInfoPatches();
    updateInfoSets();
    updateInfoZones();
    updateInfoLagrangian();

    // restore the enabled selections
    setSelectedArrayEntries
    (
        regionSelection,
        selectedEntries
    );

    if (meshChanged_)
    {
        fieldsChanged_ = true;
    }

    // Update volume, point and lagrangian fields
    updateInfoFields<fvPatchField, volMesh>
    (
        reader_->GetVolFieldSelection()
    );
    updateInfoFields<pointPatchField, pointMesh>
    (
        reader_->GetPointFieldSelection()
    );
    updateInfoLagrangianFields();

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateInfo" << endl;
    }

}


void Foam::vtkPV3Foam::updateFoamMesh()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateFoamMesh" << endl;
        printMemory();
    }

    if (!reader_->GetCacheMesh())
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

        meshChanged_ = true;
    }
    else
    {
        if (debug)
        {
            Info<< "Using existing Foam mesh" << endl;
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateFoamMesh" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::Update
(
    vtkMultiBlockDataSet* output,
    vtkMultiBlockDataSet* lagrangianOutput
)
{
    if (debug)
    {
        cout<< "<beg> Foam::vtkPV3Foam::Update - output with "
            << output->GetNumberOfBlocks() << " and "
            << lagrangianOutput->GetNumberOfBlocks() << " blocks\n";

        output->Print(cout);
        lagrangianOutput->Print(cout);
        printMemory();
    }
    reader_->UpdateProgress(0.1);

    // Set up region selection(s)
    updateRegionStatus();

    // Update the Foam mesh
    updateFoamMesh();
    reader_->UpdateProgress(0.2);

    // Convert meshes
    int blockNo = 0;

    convertMeshVolume(output, blockNo);
    convertMeshPatches(output, blockNo);
    reader_->UpdateProgress(0.4);

    if (reader_->GetIncludeZones())
    {
        convertMeshCellZones(output, blockNo);
        convertMeshFaceZones(output, blockNo);
        convertMeshPointZones(output, blockNo);
    }

    if (reader_->GetIncludeSets())
    {
        convertMeshCellSets(output, blockNo);
        convertMeshFaceSets(output, blockNo);
        convertMeshPointSets(output, blockNo);
    }

    blockNo = 0;
    convertMeshLagrangian(lagrangianOutput, blockNo);

    reader_->UpdateProgress(0.8);

    // Update fields
    convertVolFields(output);
    convertPointFields(output);
    convertLagrangianFields(lagrangianOutput);
    reader_->UpdateProgress(1.0);

    meshChanged_ = fieldsChanged_ = false;
}


double* Foam::vtkPV3Foam::findTimes(int& nTimeSteps)
{
    int nTimes = 0;
    double* tsteps = NULL;

    if (dbPtr_.valid())
    {
        Time& runTime = dbPtr_();
        instantList timeLst = runTime.times();

        // always skip "constant" time, unless there are no other times
        nTimes = timeLst.size();
        label timeI = 0;

        if (nTimes > 1)
        {
            timeI = 1;
            --nTimes;
        }

        if (nTimes)
        {
            tsteps = new double[nTimes];
            for (label stepI = 0; stepI < nTimes; ++stepI, ++timeI)
            {
                tsteps[stepI] = timeLst[timeI].value();
            }
        }
    }
    else
    {
        if (debug)
        {
            cout<< "no valid dbPtr:\n";
        }
    }

    // vector length returned via the parameter
    nTimeSteps = nTimes;

    return tsteps;
}


void Foam::vtkPV3Foam::addPatchNames(vtkRenderer* renderer)
{
    // Remove any patch names previously added to the renderer
    removePatchNames(renderer);

    // get the display patches, strip off any suffix
    const stringList selectedPatches = getSelectedArrayEntries
    (
        reader_->GetRegionSelection(),
        regionInfoPatches_,
        true
    );

    if (!selectedPatches.size())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::addPatchNames" << nl
            <<"... add patches: " << selectedPatches << endl;
    }

    const polyBoundaryMesh& pbMesh = meshPtr_->boundaryMesh();

    // Find the total number of zones
    // Each zone will take the patch name
    // Number of zones per patch ... zero zones should be skipped
    labelList nZones(pbMesh.size(), 0);

    // Per global zone number the average face centre position
    DynamicList<point> zoneCentre(pbMesh.size());

    if (debug)
    {
        Info<< "... determining patch zones" << endl;
    }

    // Loop through all patches to determine zones, and centre of each zone
    forAll(pbMesh, patchI)
    {
        const polyPatch& pp = pbMesh[patchI];

        // Only include the patch if it is selected
        bool isSelected = false;
        forAll(selectedPatches, elemI)
        {
            if (pp.name() == selectedPatches[elemI])
            {
                isSelected = true;
                break;
            }
        }

        if (isSelected)
        {
            const labelListList& edgeFaces = pp.edgeFaces();
            const vectorField& n = pp.faceNormals();

            boolList featEdge(pp.nEdges(), false);

            forAll(edgeFaces, edgeI)
            {
                const labelList& eFaces = edgeFaces[edgeI];

                if (eFaces.size() != 2)
                {
                    featEdge[edgeI] = true;
                }
                else if (mag(n[eFaces[0]] & n[eFaces[1]]) < 0.5)
                {
                    featEdge[edgeI] = true;
                }
            }

            // Do topological analysis of patch. Determine disconnected regions
            patchZones pZones(pp, featEdge);

            nZones[patchI] = pZones.nZones();

            labelList zoneNFaces(pZones.nZones(), 0);

            // Save start of information for current patch
            label patchStart = zoneCentre.size();

            // Create storage for additional zone centres
            forAll(zoneNFaces, zoneI)
            {
                zoneCentre.append(vector::zero);
            }

            // Do averaging per individual zone

            forAll(pp, faceI)
            {
                label zoneI = pZones[faceI];
                zoneCentre[patchStart+zoneI] += pp[faceI].centre(pp.points());
                zoneNFaces[zoneI]++;
            }

            for (label i=0; i<nZones[patchI]; i++)
            {
                zoneCentre[patchStart + i] /= zoneNFaces[i];
            }
        }
    }
    zoneCentre.shrink();

    if (debug)
    {
        Info<< "patch zone centres = " << zoneCentre << nl
            << "zones per patch = " << nZones << endl;
    }

    // Set the size of the patch labels to max number of zones
    patchTextActorsPtrs_.setSize(zoneCentre.size());

    if (debug)
    {
        Info<< "constructing patch labels" << endl;
    }

    label globalZoneI = 0;
    forAll(pbMesh, patchI)
    {
        const polyPatch& pp = pbMesh[patchI];

        // Only selected patches will have a non-zero number of zones
        for (label i=0; i<nZones[patchI]; i++)
        {
            if (debug)
            {
                Info<< "patch name = " << pp.name() << nl
                    << "anchor = " << zoneCentre[globalZoneI] << nl
                    << "globalZoneI = " << globalZoneI << endl;
            }

            vtkTextActor* txt = vtkTextActor::New();

            txt->SetInput(pp.name().c_str());

            // Set text properties
            vtkTextProperty* tprop = txt->GetTextProperty();
            tprop->SetFontFamilyToArial();
            tprop->BoldOff();
            tprop->ShadowOff();
            tprop->SetLineSpacing(1.0);
            tprop->SetFontSize(12);
            tprop->SetColor(1.0, 0.0, 0.0);
            tprop->SetJustificationToCentered();

            // Set text to use 3-D world co-ordinates
            txt->GetPositionCoordinate()->SetCoordinateSystemToWorld();

            txt->GetPositionCoordinate()->SetValue
            (
                zoneCentre[globalZoneI].x(),
                zoneCentre[globalZoneI].y(),
                zoneCentre[globalZoneI].z()
            );

            // Add text to each renderer
            renderer->AddViewProp(txt);

            // Maintain a list of text labels added so that they can be
            // removed later
            patchTextActorsPtrs_[globalZoneI] = txt;

            globalZoneI++;
        }
    }

    // Resize the patch names list to the actual number of patch names added
    patchTextActorsPtrs_.setSize(globalZoneI);

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::addPatchNames" << endl;
    }
}


void Foam::vtkPV3Foam::removePatchNames(vtkRenderer* renderer)
{
    forAll(patchTextActorsPtrs_, patchI)
    {
        renderer->RemoveViewProp(patchTextActorsPtrs_[patchI]);
        patchTextActorsPtrs_[patchI]->Delete();
    }
    patchTextActorsPtrs_.clear();
}


void Foam::vtkPV3Foam::PrintSelf(ostream& os, vtkIndent indent) const
{
    os  << indent << "Number of meshes: " << nMesh_ << "\n";
    os  << indent << "Number of nodes: "
        << (meshPtr_ ? meshPtr_->nPoints() : 0) << "\n";

    os  << indent << "Number of cells: "
        << (meshPtr_ ? meshPtr_->nCells() : 0) << "\n";

    os  << indent << "Number of available time steps: "
        << (dbPtr_.valid() ? dbPtr_().times().size() : 0) << endl;
}


// parse these bits of info from /proc/meminfo (Linux)
//
// MemTotal:      2062660 kB
// MemFree:       1124400 kB
//
// used = MemTotal - MemFree is what the free(1) uses.
//
void Foam::vtkPV3Foam::printMemory()
{
    const char* meminfo = "/proc/meminfo";

    if (exists(meminfo))
    {
        IFstream is(meminfo);
        label memTotal = 0;
        label memFree = 0;

        string line;

        while (is.getLine(line).good())
        {
            char tag[32];
            int value;

            if (sscanf(line.c_str(), "%30s %d", tag, &value) == 2)
            {
                if (!strcmp(tag, "MemTotal:"))
                {
                    memTotal = value;
                }
                else if (!strcmp(tag, "MemFree:"))
                {
                    memFree = value;
                }
            }
        }

        Info << "memUsed: " << (memTotal - memFree) << " kB\n";
    }
}

// ************************************************************************* //
