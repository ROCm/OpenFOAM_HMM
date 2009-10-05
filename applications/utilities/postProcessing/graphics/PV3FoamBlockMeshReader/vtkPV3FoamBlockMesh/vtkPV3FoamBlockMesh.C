/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "vtkPV3FoamBlockMesh.H"
#include "vtkPV3FoamBlockMeshReader.h"

// Foam includes
#include "blockMesh.H"
#include "Time.H"
#include "patchZones.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::vtkPV3FoamBlockMesh, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkPV3FoamBlockMesh::resetCounters()
{
    // Reset mesh part ids and sizes
    partInfoBlocks_.reset();
    partInfoCorners_.reset();
}


void Foam::vtkPV3FoamBlockMesh::updateMeshPartsStatus()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3FoamBlockMesh::updateMeshPartsStatus" << endl;
    }

    vtkDataArraySelection* selection = reader_->GetPartSelection();
    label nElem = selection->GetNumberOfArrays();

    if (partStatus_.size() != nElem)
    {
        partStatus_.setSize(nElem);
        partStatus_ = false;
    }

    // this needs fixing if we wish to re-use the datasets
    partDataset_.setSize(nElem);
    partDataset_ = -1;

    // Read the selected mesh parts (blocks only) and add to list
    forAll(partStatus_, partId)
    {
        const int setting = selection->GetArraySetting(partId);

        if (partStatus_[partId] != setting)
        {
            partStatus_[partId] = setting;
        }

        if (debug)
        {
            Info<< "  part[" << partId << "] = "
                << partStatus_[partId]
                << " : " << selection->GetArrayName(partId) << endl;
        }
    }
    if (debug)
    {
        Info<< "<end> Foam::vtkPV3FoamBlockMesh::updateMeshPartsStatus" << endl;
    }
}


void Foam::vtkPV3FoamBlockMesh::updateInfoBlocks()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3FoamBlockMesh::updateInfoBlocks"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();
    partInfoBlocks_ = partSelection->GetNumberOfArrays();

    int nBlocks = meshPtr_->size();

    for (int blockI = 0; blockI < nBlocks; ++blockI)
    {
        // Add blockId to GUI list
        partSelection->AddArray
        (
            Foam::name(blockI).c_str()
        );
    }

    partInfoBlocks_ += nBlocks;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(partSelection);

        Info<< "<end> Foam::vtkPV3FoamBlockMesh::updateInfoBlocks" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPV3FoamBlockMesh::vtkPV3FoamBlockMesh
(
    const char* const FileName,
    vtkPV3FoamBlockMeshReader* reader
)
:
    reader_(reader),
    dbPtr_(NULL),
    meshPtr_(NULL),
    partInfoBlocks_("block"),
    partInfoCorners_("corners")
{
    if (debug)
    {
        Info<< "Foam::vtkPV3FoamBlockMesh::vtkPV3FoamBlockMesh - "
            << FileName << endl;
    }

    // avoid argList and get rootPath/caseName directly from the file
    fileName fullCasePath(fileName(FileName).path());

    if (!isDir(fullCasePath))
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
        const fileName globalCase = fullCasePath.path();

        setEnv("FOAM_CASE", globalCase, true);
        setEnv("FOAM_CASENAME", globalCase.name(), true);
    }
    else
    {
        setEnv("FOAM_CASE", fullCasePath, true);
        setEnv("FOAM_CASENAME", fullCasePath.name(), true);
    }

    // look for 'case{region}.OpenFOAM'
    // could be stringent and insist the prefix match the directory name...
    // Note: cannot use fileName::name() due to the embedded '{}'
    string caseName(fileName(FileName).lessExt());

    if (debug)
    {
        Info<< "fullCasePath=" << fullCasePath << nl
            << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
            << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << endl;
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

Foam::vtkPV3FoamBlockMesh::~vtkPV3FoamBlockMesh()
{
    if (debug)
    {
        Info<< "<end> Foam::vtkPV3FoamBlockMesh::~vtkPV3FoamBlockMesh" << endl;
    }

    delete meshPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3FoamBlockMesh::updateInfo()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3FoamBlockMesh::updateInfo"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "] " << endl;
    }

    resetCounters();

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();

    // enable 'internalMesh' on the first call
    // or preserve the enabled selections
    stringList enabledEntries;
    bool firstTime = false;
    if (!partSelection->GetNumberOfArrays() && !meshPtr_)
    {
        firstTime = true;
    }
    else
    {
        enabledEntries = getSelectedArrayEntries(partSelection);
    }

    // Clear current mesh parts list
    partSelection->RemoveAllArrays();

    // need a blockMesh
    updateFoamMesh();

    // Update mesh parts list - add corrner points at the bottom
    updateInfoBlocks();

    // restore the enabled selections
    if (!firstTime)
    {
        setSelectedArrayEntries(partSelection, enabledEntries);
    }


    if (debug)
    {
        Info<< "<end> Foam::vtkPV3FoamBlockMesh::updateInfo" << endl;
    }
}


void Foam::vtkPV3FoamBlockMesh::updateFoamMesh()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3FoamBlockMesh::updateFoamMesh" << endl;
    }

    // Check to see if the FOAM mesh has been created
    if (!meshPtr_)
    {
        if (debug)
        {
            Info<< "Creating blockMesh at time=" << dbPtr_().timeName()
                << endl;
        }

        IOdictionary meshDict
        (
            IOobject
            (
                "blockMeshDict",
                dbPtr_().constant(),
                polyMesh::meshSubDir,
                dbPtr_(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        meshPtr_ = new blockMesh(meshDict);
    }


    if (debug)
    {
        Info<< "<end> Foam::vtkPV3FoamBlockMesh::updateFoamMesh" << endl;
    }
}


void Foam::vtkPV3FoamBlockMesh::Update
(
    vtkMultiBlockDataSet* output
)
{
    reader_->UpdateProgress(0.1);

    // Set up mesh parts selection(s)
    updateMeshPartsStatus();

    reader_->UpdateProgress(0.2);

    // Update the Foam mesh
    updateFoamMesh();
    reader_->UpdateProgress(0.5);

    // Convert meshes - start port0 at block=0
    int blockNo = 0;

    convertMeshBlocks(output, blockNo);
    convertMeshCorners(output, blockNo);

    reader_->UpdateProgress(0.8);

}


void Foam::vtkPV3FoamBlockMesh::CleanUp()
{
    reader_->UpdateProgress(1.0);
}


void Foam::vtkPV3FoamBlockMesh::renderPointNumbers
(
    vtkRenderer* renderer,
    const bool show
)
{
    // always remove old actors first

    forAll(pointNumberTextActorsPtrs_, pointI)
    {
        renderer->RemoveViewProp(pointNumberTextActorsPtrs_[pointI]);
        pointNumberTextActorsPtrs_[pointI]->Delete();
    }
    pointNumberTextActorsPtrs_.clear();

    if (show && meshPtr_)
    {
        const pointField& cornerPts = meshPtr_->blockPointField();

        pointNumberTextActorsPtrs_.setSize(cornerPts.size());
        forAll(cornerPts, pointI)
        {
            vtkTextActor* txt = vtkTextActor::New();

            txt->SetInput(Foam::name(pointI).c_str());

            // Set text properties
            vtkTextProperty* tprop = txt->GetTextProperty();
            tprop->SetFontFamilyToArial();
            tprop->BoldOn();
            tprop->ShadowOff();
            tprop->SetLineSpacing(1.0);
            tprop->SetFontSize(12);
            tprop->SetColor(1.0, 0.0, 0.0);
            tprop->SetJustificationToCentered();

            // Set text to use 3-D world co-ordinates
            txt->GetPositionCoordinate()->SetCoordinateSystemToWorld();

            txt->GetPositionCoordinate()->SetValue
            (
                cornerPts[pointI].x(),
                cornerPts[pointI].y(),
                cornerPts[pointI].z()
            );

            // Add text to each renderer
            renderer->AddViewProp(txt);

            // Maintain a list of text labels added so that they can be
            // removed later
            pointNumberTextActorsPtrs_[pointI] = txt;
        }
    }
}



void Foam::vtkPV3FoamBlockMesh::PrintSelf(ostream& os, vtkIndent indent) const
{
#if 0
    os  << indent << "Number of nodes: "
        << (meshPtr_ ? meshPtr_->nPoints() : 0) << "\n";

    os  << indent << "Number of cells: "
        << (meshPtr_ ? meshPtr_->nCells() : 0) << "\n";

    os  << indent << "Number of available time steps: "
        << (dbPtr_.valid() ? dbPtr_().times().size() : 0) << endl;
#endif
}

// ************************************************************************* //
