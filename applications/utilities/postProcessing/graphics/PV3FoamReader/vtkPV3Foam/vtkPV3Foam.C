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
#include "argList.H"
#include "JobInfo.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOobjectList.H"
#include "patchZones.H"
#include "vtkPV3FoamReader.h"

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

#include "vtkPV3FoamAddFields.H"
#include "vtkPV3FoamUpdateInformationFields.H"


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
    selectInfoVolume_.reset();
    selectInfoPatches_.reset();
    selectInfoLagrangian_.reset();
    selectInfoCellZones_.reset();
    selectInfoFaceZones_.reset();
    selectInfoPointZones_.reset();
    selectInfoCellSets_.reset();
    selectInfoFaceSets_.reset();
    selectInfoPointSets_.reset();
}


void Foam::vtkPV3Foam::initializeTime()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::initializeTime" << endl;
    }

    Time& runTime = dbPtr_();

    // Get times list
    instantList times = runTime.times();

    vtkDataArraySelection* arraySelection = reader_->GetTimeSelection();

    // only execute this if there is a mismatch between
    // the times available according to FOAM and according to VTK
    int nArrays = arraySelection->GetNumberOfArrays();

    if (nArrays && nArrays == times.size() - 1)
    {
        return;
    }

    // "constant" is implicit - skip it
    // All the time selections are enabled by default
    for (label timeI = 1; timeI < times.size(); ++timeI)
    {
        arraySelection->AddArray
        (
            times[timeI].name().c_str()
        );
    }
}


bool Foam::vtkPV3Foam::setTime(const double& requestedTime)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::setTime("
            << requestedTime << ")" << endl;
    }

    Time& runTime = dbPtr_();

    // Get times list
    instantList times = runTime.times();

    // logic as per "checkTimeOption.H"
    bool found = false;
    int nearestIndex = -1;
    scalar nearestDiff = Foam::GREAT;

    forAll (times, timeIndex)
    {
        if (times[timeIndex].name() == "constant") continue;

        scalar diff = fabs(times[timeIndex].value() - requestedTime);
        if (diff < nearestDiff)
        {
            nearestDiff = diff;
            nearestIndex = timeIndex;
        }
    }

    if (nearestIndex == -1)
    {
        nearestIndex = 0;
        found = false;
    }
    else
    {
        found = true;
    }

    if (debug)
    {
        Info<< "Selecting time " << times[nearestIndex].name() << endl;
    }

    runTime.setTime(times[nearestIndex], nearestIndex);
    return found;
}


void Foam::vtkPV3Foam::updateSelectedRegions()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateSelectedRegions" << endl;
    }

    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    const label nRegions = arraySelection->GetNumberOfArrays();

    selectedRegions_.setSize(nRegions);
    selectedRegionDatasetIds_.setSize(nRegions);

    // Read the selected patches and add to the region list
    for (int regionId=0; regionId < nRegions; ++regionId)
    {
        selectedRegions_[regionId] = arraySelection->GetArraySetting(regionId);
        selectedRegionDatasetIds_[regionId] = -1;

        if (debug)
        {
            Info<< "region " << regionId
                << " = " << selectedRegions_[regionId] << endl;
        }
    }
}


Foam::stringList Foam::vtkPV3Foam::getSelectedArrayEntries
(
    vtkDataArraySelection* arraySelection,
    const bool firstWord
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::getSelectedArrayEntries" << endl;
    }

    stringList selections(arraySelection->GetNumberOfArrays());
    label nElem = 0;
    forAll (selections, elemI)
    {
        if (arraySelection->GetArraySetting(elemI))
        {
            if (firstWord)
            {
                selections[nElem] = getFirstWord
                (
                    arraySelection->GetArrayName(elemI)
                );
            }
            else
            {
                selections[nElem] = arraySelection->GetArrayName(elemI);
            }
            ++nElem;
        }
    }

    selections.setSize(nElem);
    if (debug)
    {
        Info<< "Active array: " << selections << endl;
    }

    return selections;
}


Foam::stringList Foam::vtkPV3Foam::getSelectedArrayEntries
(
    vtkDataArraySelection* arraySelection,
    const selectionInfo& selector,
    const bool firstWord
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::getSelectedArrayEntries" << endl;
    }

    stringList selections(selector.size());
    label nElem = 0;

    for
    (
        int regionId = selector.start();
        regionId < selector.end();
        ++regionId
    )
    {
        if (arraySelection->GetArraySetting(regionId))
        {
            if (firstWord)
            {
                selections[nElem] = getFirstWord
                (
                    arraySelection->GetArrayName(regionId)
                );
            }
            else
            {
                selections[nElem] = arraySelection->GetArrayName(regionId);
            }

            ++nElem;
        }
    }

    selections.setSize(nElem);
    if (debug)
    {
        Info<< "Active array: " << selections << endl;
    }

    return selections;
}


void Foam::vtkPV3Foam::setSelectedArrayEntries
(
    vtkDataArraySelection* arraySelection,
    const stringList& selections
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::setSelectedArrayEntries" << endl;
    }
    const label nEntries = arraySelection->GetNumberOfArrays();

    // Reset all current entries to 'not selected'
    arraySelection->DisableAllArrays();

    // Loop through entries, setting values from selectedEntries
    forAll (selections, elemI)
    {
        if (debug)
        {
            Info<< "selections[" << elemI << "] = " << selections[elemI]
                << endl;
        }

        for (label i=0; i<nEntries; i++)
        {
            string arrayName = arraySelection->GetArrayName(i);

            if (arrayName == selections[elemI])
            {
                if (debug)
                {
                    Info<< "enabling array: " << arrayName << " Index = "
                        << i
                        << endl;
                }

                arraySelection->EnableArray
                (
                    arrayName.c_str()
                );
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPV3Foam::vtkPV3Foam
(
    const char* const FileName,
    vtkPV3FoamReader* reader,
    vtkMultiBlockDataSet* output
)
:
    reader_(reader),
    output_(output),
    dbPtr_(NULL),
    meshPtr_(NULL),
    selectInfoVolume_(VOLUME, "unzoned"),
    selectInfoPatches_(PATCHES, "patches"),
    selectInfoLagrangian_(LAGRANGIAN, "lagrangian"),
    selectInfoCellZones_(CELLZONE, "cellZone"),
    selectInfoFaceZones_(FACEZONE, "faceZone"),
    selectInfoPointZones_(POINTZONE, "pointZone"),
    selectInfoCellSets_(CELLSET, "cellSet"),
    selectInfoFaceSets_(FACESET, "faceSet"),
    selectInfoPointSets_(POINTSET, "pointSet"),
    patchTextActorsPtrs_(0),
    nMesh_(0)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::vtkPV3Foam with "
            << FileName << endl;
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
    setEnv("FOAM_CASE", fullCasePath, true);

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

    if (debug)
    {
        cout<< "constructed with output: ";
        output_->Print(cout);
    }

    resetCounters();

    // Set initial cloud name
    // TODO - TEMPORARY MEASURE UNTIL CAN PROCESS MULTIPLE CLOUDS
    cloudName_ = "";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkPV3Foam::~vtkPV3Foam()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::~vtkPV3Foam" << endl;
    }

    if (meshPtr_)
    {
        delete meshPtr_;
        meshPtr_ = NULL;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::UpdateInformation()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::UpdateInformation" << nl
            << "TimeStep = " << reader_->GetTimeStep() << endl;
    }

    resetCounters();

    // preserve the currently selected values
    const stringList selectedEntries = getSelectedArrayEntries
    (
        reader_->GetRegionSelection()
    );
    // Clear current region list/array
    reader_->GetRegionSelection()->RemoveAllArrays();

    initializeTime();

    // Update region array
    updateInformationInternalMesh();

    updateInformationLagrangian();

    updateInformationPatches();

    if (reader_->GetIncludeSets())
    {
        updateInformationSets();
    }

    if (reader_->GetIncludeZones())
    {
        updateInformationZones();
    }

    // Update region selection with the data just read in
    setSelectedArrayEntries
    (
        reader_->GetRegionSelection(),
        selectedEntries
    );

    // Update volField array
    updateInformationFields<fvPatchField, volMesh>
    (
        reader_->GetVolFieldSelection()
    );

    // Update pointField array
    updateInformationFields<pointPatchField, pointMesh>
    (
        reader_->GetPointFieldSelection()
    );

    // Update lagrangian field array
    updateInformationLagrangianFields();
}


void Foam::vtkPV3Foam::Update
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        cout<< "entered Foam::vtkPV3Foam::Update" << nl
             <<"Update\n";
        output->Print(cout);

        cout<<"Internally:\n";
        output_->Print(cout);

        cout<< " has " << output_->GetNumberOfBlocks() << " blocks\n";
    }



    // Set up region selection(s)
    updateSelectedRegions();

    // Update the Foam mesh
    updateFoamMesh();

    // Convert meshes
    convertMeshVolume(output);

    convertMeshLagrangian(output);

    convertMeshPatches(output);

    if (reader_->GetIncludeZones())
    {
        convertMeshCellZones(output);
        convertMeshFaceZones(output);
        convertMeshPointZones(output);
    }

    if (reader_->GetIncludeSets())
    {
        convertMeshCellSet(output);
        convertMeshFaceSet(output);
        convertMeshPointSet(output);
    }

    // Update fields
    updateVolFields(output);

    updatePointFields(output);

    updateLagrangianFields(output);

    if (debug)
    {
        Info<< "Number of data sets after update" << nl
            << "  VOLUME = "
            << GetNumberOfDataSets(output, selectInfoVolume_) << nl
            << "  PATCHES = "
            << GetNumberOfDataSets(output, selectInfoPatches_) << nl
            << "  LAGRANGIAN = "
            << GetNumberOfDataSets(output, selectInfoLagrangian_) << nl
            << "  CELLZONE = "
            << GetNumberOfDataSets(output, selectInfoCellZones_) << nl
            << "  FACEZONE = "
            << GetNumberOfDataSets(output, selectInfoFaceZones_) << nl
            << "  POINTZONE = "
            << GetNumberOfDataSets(output, selectInfoPointZones_) << nl
            << "  CELLSET = "
            << GetNumberOfDataSets(output, selectInfoCellSets_) << nl
            << "  FACESET = "
            << GetNumberOfDataSets(output, selectInfoFaceSets_) << nl
            << "  POINTSET = "
            << GetNumberOfDataSets(output, selectInfoPointSets_) << nl;

        // traverse blocks:

        int nBlocks = output->GetNumberOfBlocks();
        Info << "nBlocks = " << nBlocks << endl;

        cout<< "done Update\n";
        output_->Print(cout);
        cout<< " has " << output_->GetNumberOfBlocks() << " blocks\n";
        output_->GetInformation()->Print(cout);

        cout <<"ShouldIReleaseData :" << output_->ShouldIReleaseData() << "\n";
    }
}


double* Foam::vtkPV3Foam::timeSteps(int& nTimeSteps)
{
    int nTimes = 0;
    double* ts = NULL;

    vtkDataArraySelection* arraySelection = reader_->GetTimeSelection();

    if (dbPtr_.valid())
    {
        Time& runTime = dbPtr_();

        instantList times = runTime.times();
        List<bool> selected = List<bool>(times.size(), false);

        const label nSelectedTimes = arraySelection->GetNumberOfArrays();

        for (int i = 0; i < nSelectedTimes; ++i)
        {
            // always skip "constant" time
            const int timeI = i + 1;
            if
            (
                arraySelection->GetArraySetting(i)
             && timeI < times.size()
            )
            {
#if 0
                Info<<"timeSelection["
                    << i
                    <<"] = "
                    << arraySelection->GetArraySetting(i)
                    << " is "
                    << arraySelection->GetArrayName(i) << endl;
#endif
                selected[timeI] = true;
                ++nTimes;
            }
        }

        if (debug)
        {
            Info<< "selected " << nTimes << " times ";
            Info<< "found " << times.size() << " times: (";
            forAll(times, timeI)
            {
                Info<< " " << times[timeI].value();
            }
            Info<< " )" << endl;
        }

        if (nTimes)
        {
            ts = new double[nTimes];
            int stepI = 0;

            forAll(selected, selectI)
            {
                if (selected[selectI])
                {
                    ts[stepI] = times[selectI].value();
                    stepI++;
                }
            }
        }
    }
    else
    {
        if (debug)
        {
            Info<< "no valid dbPtr:" <<endl;
        }
    }

    // return length via the parameter
    nTimeSteps = nTimes;

    return ts;
}


void Foam::vtkPV3Foam::addPatchNames(vtkRenderer* renderer)
{
    // Remove any patch names previously added to the renderer
    removePatchNames(renderer);

    if (debug)
    {
        Info<< "addPatchNames()" << endl;
    }

    const fvMesh& mesh = *meshPtr_;
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    const selectionInfo& selector = selectInfoPatches_;

    // the currently selected patches, strip off any suffix
    const stringList selectedPatches = getSelectedArrayEntries
    (
        reader_->GetRegionSelection(),
        selector,
        true
    );

    if (debug)
    {
        Info<<"patches: " << selectedPatches <<endl;
    }

    // Find the total number of zones
    // Each zone will take the patch name

    // Number of zones per patch ... zero zones should be skipped
    labelList nZones(pbMesh.size(), 0);

    // Per global zone number the average face centre position
    DynamicList<point> zoneCentre(pbMesh.size());

    if (debug)
    {
        Info<< "determining patch zones" << endl;
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
        Info<< "patch zone centres = " << zoneCentre << endl;
        Info<< "zones per patch = " << nZones << endl;
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
                Info<< "patch name = " << pp.name() << endl;
                Info<< "anchor = " << zoneCentre[globalZoneI] << endl;
                Info<< "globalZoneI = " << globalZoneI << endl;
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
}


void Foam::vtkPV3Foam::removePatchNames(vtkRenderer* renderer)
{
    if (debug)
    {
        Info<< "removePatchNames()" << endl;
    }

    forAll(patchTextActorsPtrs_, patchI)
    {
        renderer->RemoveViewProp(patchTextActorsPtrs_[patchI]);
        patchTextActorsPtrs_[patchI]->Delete();
    }
    patchTextActorsPtrs_.setSize(0);
}


int Foam::vtkPV3Foam::numberOfAvailableTimes() const
{
    if (dbPtr_.valid())
    {
        return dbPtr_().times().size();
    }
    else
    {
        return 0;
    }
}


int Foam::vtkPV3Foam::numberOfPoints()
{
    if (meshPtr_)
    {
        return meshPtr_->nPoints();
    }
    else
    {
        return 0;
    }
}


int Foam::vtkPV3Foam::numberOfCells()
{
    if (meshPtr_)
    {
        return meshPtr_->nCells();
    }
    {
        return 0;
    }
}


int Foam::vtkPV3Foam::numberOfMeshes() const
{
    return nMesh_;
}


// ************************************************************************* //
