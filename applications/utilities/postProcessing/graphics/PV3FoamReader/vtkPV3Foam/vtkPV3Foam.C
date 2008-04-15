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

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "argList.H"
#include "Time.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "vtkPV3FoamReader.h"
#include "patchZones.H"

// VTK includes
#include "vtkCharArray.h"
#include "vtkDataArraySelection.h"
#include "vtkFieldData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::vtkPV3Foam, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "vtkPV3FoamAddFields.H"

#include "vtkPV3FoamUpdateInformationFields.H"


void Foam::vtkPV3Foam::resetCounters()
{
    // Reset data size counts
    lagrangianDataSize_ = 0;
    patchDataSize_ = 0;
    cellSetDataSize_ = 0;
    faceSetDataSize_ = 0;
    pointSetDataSize_ = 0;

    // Reset region ids
    idRegionVolume_ = -1;
    idRegionLagrangian_ = -1;
    idRegionPatches_ = -1;
    idRegionCellSets_ = -1;
    idRegionFaceSets_ = -1;
    idRegionPointSets_ = -1;
}


void Foam::vtkPV3Foam::SetName
(
    vtkUnstructuredGrid* vtkMesh,
    const char* name
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::setName" << endl;
    }
    vtkCharArray* nmArray =  vtkCharArray::New();
    nmArray->SetName("Name");
    size_t len = strlen(name);
    nmArray->SetNumberOfTuples(static_cast<vtkIdType>(len) + 1);
    char* copy = nmArray->GetPointer(0);
    memcpy(copy, name, len);
    copy[len] = '\0';
    vtkMesh->GetFieldData()->AddArray(nmArray);
    nmArray->Delete();
}


void Foam::vtkPV3Foam::setSelectedTime
(
    Time& runTime,
    vtkPV3FoamReader* reader
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::setSelectedTime" << endl;
    }

    // Get times list
    instantList times = runTime.times();
    int timeIndex = min(max(reader->GetTimeStep() + 1, 0), times.size() - 1);

    // If this is the first call timeIndex will be 0 ("constant")
    // so reset to the first time step if one exists and deselect every
    // element of the selection array
    if (timeIndex == 0)
    {
        timeIndex = min(1, times.size() - 1);
        reader->GetTimeSelection()->DisableAllArrays();
    }

    label selectedTimeIndex = -1;
    const label nSelectedTimes =
        reader->GetTimeSelection()->GetNumberOfArrays();

    // Choose latest time
    for (label i=nSelectedTimes - 1; i>=0; i--)
    {
        if (reader->GetTimeSelection()->GetArraySetting(i))
        {
            word timeName = string::validate<word>
            (
                reader->GetTimeSelection()->GetArrayName(i)
            );

            forAll(times, j)
            {
                if (times[j].name() == timeName)
                {
                    selectedTimeIndex = j;
                    break;
                }
            }
            break;
        }
    }

    if (selectedTimeIndex != -1)
    {
        timeIndex = min(selectedTimeIndex, times.size() - 1);
    }

    if (debug)
    {
        Info<< "Selecting time " << times[timeIndex].name() << endl;
    }

    runTime.setTime(times[timeIndex], timeIndex);

    times = runTime.times();

    // set time step slider range
    reader->SetTimeStepRange(0, max(times.size() - 2, 0));
    reader->SetTimeStep(timeIndex); // TODO - does nothing???

    // reset the time steps ...
    reader->GetTimeSelection()->RemoveAllArrays();

    forAll (times, timeI)
    {
        reader->GetTimeSelection()->AddArray(times[timeI].name().c_str());
    }


    // Disable all the time selections (which are all selected by default) ...
    reader->GetTimeSelection()->DisableAllArrays();

    // But maintain the selection made previously
    if (selectedTimeIndex != -1 && selectedTimeIndex < times.size())
    {
        reader->GetTimeSelection()->EnableArray
            (times[selectedTimeIndex].name().c_str());
    }
}


void Foam::vtkPV3Foam::updateSelectedRegions()
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::updateSelectedRegions" << endl;
    }

    const label nRegions = reader_->GetRegionSelection()->GetNumberOfArrays();

    selectedRegions_.setSize(nRegions);
    selectedRegionDatasetIds_.setSize(nRegions);

    // Read the selected patches and add to the region list
    for (int i=0; i<nRegions; i++)
    {
        selectedRegions_[i] = reader_->GetRegionSelection()
            ->GetArraySetting(i);
        selectedRegionDatasetIds_[i] = -1;
        if (debug)
        {
            Info<< "region" << i << " selected = " << selectedRegions_[i]
                << endl;
        }
    }
}


Foam::wordList Foam::vtkPV3Foam::getSelectedArrayEntries
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::getSelectedArrayEntries" << endl;
    }
    const label nArrays = arraySelection->GetNumberOfArrays();
    wordList selections(nArrays);
    label j = 0;
    for (label i=0; i<nArrays; i++)
    {
        if (arraySelection->GetArraySetting(i))
        {
            selections[j] = arraySelection->GetArrayName(i);
            if (debug)
            {
                Info<< "Active array: " << selections[j] << " index: " << j
                    << endl;
            }
            j++;
        }
    }

    // Fill non-selected arrays with 'default' value
    for (label k=j; k<nArrays; k++)
    {
        selections[k] = "undefinedArrayEntry";
    }

    return selections;
}


void Foam::vtkPV3Foam::setSelectedArrayEntries
(
    vtkDataArraySelection* arraySelection,
    const wordList& selectedEntries
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
    forAll (selectedEntries, entryI)
    {
        if (debug)
        {
            Info<< "selectedEntries[entryI] = " << selectedEntries[entryI]
                << endl;
        }

        for (label i=0; i<nEntries; i++)
        {
            const word arrayName = string::validate<word>
            (
                arraySelection->GetArrayName(i)
            );
            if (arrayName == selectedEntries[entryI])
            {
                if (debug)
                {
                    Info<< "enabling array: " << arrayName << " Index = " << i
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
    argsPtr_(NULL),
    dbPtr_(NULL),
    meshPtr_(NULL),
    patchTextActorsPtrs_(0),
    nMesh_(0)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::vtkPV3Foam" << endl;
    }

    // Set paths
    fileName fullCasePath(fileName(FileName).path());

    if (!dir(fullCasePath))
    {
        return;
    }

    char* argvStrings[3];
    argvStrings[0] = new char[12];
    strcpy(argvStrings[0], "/vtkPV3Foam");
    argvStrings[1] = new char[6];
    strcpy(argvStrings[1], "-case");
    argvStrings[2] = new char[fullCasePath.size() + 1];
    strcpy(argvStrings[2], fullCasePath.c_str());

    int argc = 3;
    char** argv = &argvStrings[0];
    argsPtr_.reset(new argList(argc, argv));

    for (int i=0; i<argc; i++)
    {
        delete[] argvStrings[i];
    }

    // Create time object
    dbPtr_.reset
    (
        new Time
        (
            Time::controlDictName,
            argsPtr_().rootPath(),
            argsPtr_().caseName()
        )
    );
    dbPtr_().functionObjects().off();

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

    const wordList selectedArrayEntries = getSelectedArrayEntries
    (
        reader_->GetRegionSelection()
    );
    // Clear current region list/array
    reader_->GetRegionSelection()->RemoveAllArrays();

    // Set current time
    setSelectedTime(dbPtr_(), reader_);

    // Update region array
    updateInformationInternalMesh();

    updateInformationLagrangian();

    updateInformationPatches();

    if (reader_->GetIncludeSets())
    {
        updateInformationSets();
    }

    // Update region selection with the data just read in
    setSelectedArrayEntries
    (
        reader_->GetRegionSelection(),
        selectedArrayEntries
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
        Info<< "entered Foam::vtkPV3Foam::Update" << endl;
    }

    // Set up region selection(s)
    updateSelectedRegions();

    // Update the Foam mesh
    updateFoamMesh();

    // Convert meshes
    convertMeshVolume(output);

    convertMeshLagrangian(output);

    convertMeshPatches(output);

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
            << "    VOLUME = " << output->GetNumberOfDataSets(VOLUME) << nl
            << "    LAGRANGIAN = " << output->GetNumberOfDataSets(LAGRANGIAN)
            << nl << "    CELLSET = " << output->GetNumberOfDataSets(CELLSET)
            << nl << "    FACESET = " << output->GetNumberOfDataSets(FACESET)
            << nl << "    POINTSET = " << output->GetNumberOfDataSets(POINTSET)
            << endl;
    }
}


double* Foam::vtkPV3Foam::timeSteps()
{
    if (dbPtr_.valid())
    {
        const instantList foamTimes =
            dbPtr_().findTimes(dbPtr_().path().path());
        double* ts = new double[foamTimes.size()];

        for (int i = 0; i < foamTimes.size(); i++)
        {
            ts[i] = foamTimes[i].value();
            if (debug)
            {
                Info<< "found time = " << ts[i] << endl;
            }
        }
        return ts;
    }
    else
    {
        return NULL;
    }
}


void Foam::vtkPV3Foam::addPatchNames(vtkRenderer* renderer)
{
    if (debug)
    {
        Info<< "addPatchNames()" << endl;
    }

    // Remove any patch names previously added to the renderer
    removePatchNames(renderer);

    const fvMesh& mesh = *meshPtr_;
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    // Find the total number of zones
    // Each zone will take the patch name

    // Number of zones per patch
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
        if (reader_->GetRegionArrayStatus(pp.name().c_str()) == 1)
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

        if (reader_->GetRegionArrayStatus(pp.name().c_str()) == 1)
        {
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
                tprop->SetFontSize(20);
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
    }
    patchTextActorsPtrs_.setSize(0);
}


int Foam::vtkPV3Foam::numberOfTimeSteps()
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
