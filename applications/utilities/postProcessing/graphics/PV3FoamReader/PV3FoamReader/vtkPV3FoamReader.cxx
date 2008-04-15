/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPV3FoamReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPV3FoamReader.h"

#include "pqApplicationCore.h"
#include "pqRenderView.h"
#include "pqServerManagerModel.h"

// VTK includes
#include "vtkCallbackCommand.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkDirectory.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkRenderer.h"
#include "vtkSMRenderViewProxy.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridAlgorithm.h"

// Foam includes
#include "vtkPV3Foam.H"


vtkCxxRevisionMacro(vtkPV3FoamReader, "$Revision: 1.2$");
vtkStandardNewMacro(vtkPV3FoamReader);


vtkPV3FoamReader::vtkPV3FoamReader()
{
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName  = NULL;
    foamData_ = NULL;

    CacheMesh = 0;

    UpdateGUI = 1;
    UpdateGUIOld = 1;
    TimeStep = 0;
    TimeStepRange[0] = 0;
    TimeStepRange[1] = 0;

    ShowPatchNames = 0;

    TimeSelection = vtkDataArraySelection::New();
    PointFieldSelection = vtkDataArraySelection::New();
    RegionSelection = vtkDataArraySelection::New();
    VolFieldSelection = vtkDataArraySelection::New();
    LagrangianFieldSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPV3FoamReader::SelectionModifiedCallback
    );
    SelectionObserver->SetClientData(this);

    TimeSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    RegionSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    VolFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    PointFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    LagrangianFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
}

vtkPV3FoamReader::~vtkPV3FoamReader()
{
    vtkDebugMacro(<<"Deconstructor");

    if (foamData_)
    {
        delete foamData_;
    }

    if (FileName)
    {
        delete [] FileName;
    }

    TimeSelection->RemoveObserver(this->SelectionObserver);
    RegionSelection->RemoveObserver(this->SelectionObserver);
    VolFieldSelection->RemoveObserver(this->SelectionObserver);
    PointFieldSelection->RemoveObserver(this->SelectionObserver);
    SelectionObserver->Delete();

    TimeSelection->Delete();
    RegionSelection->Delete();
    VolFieldSelection->Delete();
    PointFieldSelection->Delete();
}


// Do everything except set the output info
int vtkPV3FoamReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"Request Information");

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    if (!foamData_)
    {
        vtkDebugMacro("RequestInformation: creating foamData_");

        vtkInformation* outInfo = outputVector->GetInformationObject(0);
        vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
        (
            outInfo->Get(vtkMultiBlockDataSet::DATA_OBJECT())
        );
        foamData_ = new Foam::vtkPV3Foam(FileName, this, output);
        foamData_->UpdateInformation();
    }
    else
    {
        vtkDebugMacro("RequestInformation: updating information");

        foamData_->UpdateInformation();
    }

    double* timeSteps = foamData_->timeSteps();
    outputVector->GetInformationObject(0)->Set
    (
        vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
        timeSteps,
        foamData_->numberOfTimeSteps()
    );
    delete timeSteps;

    return 1;
}


// Set the output info
int vtkPV3FoamReader::RequestData
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestData");

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
    (
        outInfo->Get(vtkMultiBlockDataSet::DATA_OBJECT())
    );

    if
    (
        (UpdateGUIOld == GetUpdateGUI())
     || (output->GetNumberOfDataSets(0) == 0)
    )
    {
        foamData_->Update(output);

        if (ShowPatchNames == 1)
        {
            addPatchNamesToView();
        }
        else
        {
            removePatchNamesFromView();
        }
    }
    UpdateGUIOld = GetUpdateGUI();

    return 1;
}


void vtkPV3FoamReader::addPatchNamesToView()
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();

    // Get all the pqRenderView instances
    QList<pqRenderView*> renderViews = smModel->findItems<pqRenderView*>();

    for (int viewI=0; viewI<renderViews.size(); viewI++)
    {
        foamData_->addPatchNames
        (
            renderViews[viewI]->getRenderViewProxy()->GetRenderer()
        );
    }
}


void vtkPV3FoamReader::removePatchNamesFromView()
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();

    // Get all the pqRenderView instances
    QList<pqRenderView*> renderViews = smModel->findItems<pqRenderView*>();

    for (int viewI=0; viewI<renderViews.size(); viewI++)
    {
        foamData_->removePatchNames
        (
            renderViews[viewI]->getRenderViewProxy()->GetRenderer()
        );
    }
}


void vtkPV3FoamReader::PrintSelf
(
    ostream& os,
    vtkIndent indent
)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os<< indent << "File name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";
    os<< indent << "Number of meshes: " << foamData_->numberOfMeshes() << "\n";
    os<< indent << "Number of nodes: " << foamData_->numberOfPoints() << "\n";
    os<< indent << "Number of cells: " << foamData_->numberOfCells() << "\n";
    os<< indent << "Number of time steps: " << foamData_->numberOfTimeSteps()
      << endl;
    os<< indent << "Time step range: "
      << this->TimeStepRange[0] << " - " << this->TimeStepRange[1]
      << endl;
    os<< indent << "Time step: " << this->TimeStep << endl;
    return;
}


vtkDataArraySelection* vtkPV3FoamReader::GetTimeSelection()
{
    vtkDebugMacro(<<"GetTimeSelection");

    return TimeSelection;
}


int vtkPV3FoamReader::GetNumberOfTimeArrays()
{
    vtkDebugMacro(<<"GetNumberOf TimeArrays");

    return TimeSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetTimeArrayName
(
    int index
)
{
    vtkDebugMacro(<<"GetTimeArrayName");

    return TimeSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetTimeArrayStatus
(
    const char* name
)
{
    vtkDebugMacro(<<"GetTimeArrayStatus");

    return TimeSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetTimeArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetTimeArrayStatus");

    if (status)
    {
        TimeSelection->EnableArray(name);
    }
    else
    {
        TimeSelection->DisableArray(name);
    }
}


vtkDataArraySelection* vtkPV3FoamReader::GetRegionSelection()
{
    vtkDebugMacro(<<"GetRegionSelection");

    return RegionSelection;
}


int vtkPV3FoamReader::GetNumberOfRegionArrays()
{
    vtkDebugMacro(<<"GetNumberOfRegionArrays");

    return RegionSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetRegionArrayName
(
    int index
)
{
    vtkDebugMacro(<<"GetRegionArrayName");

    return RegionSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetRegionArrayStatus
(
    const char* name
)
{
    vtkDebugMacro(<<"GetRegionArrayStatus");

    return RegionSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetRegionArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetRegionArrayStatus");

    if(status)
    {
        RegionSelection->EnableArray(name);
    }
    else
    {
        RegionSelection->DisableArray(name);
    }
}


vtkDataArraySelection* vtkPV3FoamReader::GetVolFieldSelection()
{
    vtkDebugMacro(<<"GetVolFieldSelection");

    return VolFieldSelection;
}


int vtkPV3FoamReader::GetNumberOfVolFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfVolFieldArrays");

    return VolFieldSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetVolFieldArrayName
(
    int index
)
{
    vtkDebugMacro(<<"GetVolFieldArrayName");

    return VolFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetVolFieldArrayStatus
(
    const char* name
)
{
    vtkDebugMacro(<<"GetVolFieldArrayStatus");

    return VolFieldSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetVolFieldArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetVolFieldArrayStatus");

    if (status)
    {
        VolFieldSelection->EnableArray(name);
    }
    else
    {
        VolFieldSelection->DisableArray(name);
    }
}


vtkDataArraySelection* vtkPV3FoamReader::GetPointFieldSelection()
{
    vtkDebugMacro(<<"GetPointFieldSelection");

    return PointFieldSelection;
}


int vtkPV3FoamReader::GetNumberOfPointFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfPointFieldArrays");

    return PointFieldSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetPointFieldArrayName
(
    int index
)
{
    vtkDebugMacro(<<"GetPointFieldArrayName");

    return PointFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetPointFieldArrayStatus
(
    const char* name
)
{
    vtkDebugMacro(<<"GetPointFieldArrayStatus");

    return PointFieldSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetPointFieldArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetPointFieldArrayStatus");

    if (status)
    {
        PointFieldSelection->EnableArray(name);
    }
    else
    {
        PointFieldSelection->DisableArray(name);
    }
}


vtkDataArraySelection* vtkPV3FoamReader::GetLagrangianFieldSelection()
{
    vtkDebugMacro(<<"GetLagrangianFieldSelection");

    return LagrangianFieldSelection;
}


int vtkPV3FoamReader::GetNumberOfLagrangianFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfLagrangianFieldArrays");

    return LagrangianFieldSelection->GetNumberOfArrays();
}


const char* vtkPV3FoamReader::GetLagrangianFieldArrayName
(
    int index
)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayName");

    return LagrangianFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetLagrangianFieldArrayStatus
(
    const char* name
)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayStatus");

    return LagrangianFieldSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetLagrangianFieldArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetLagrangianFieldArrayStatus");

    if (status)
    {
        LagrangianFieldSelection->EnableArray(name);
    }
    else
    {
        LagrangianFieldSelection->DisableArray(name);
    }
}


void vtkPV3FoamReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkPV3FoamReader*>(clientdata)->SelectionModified();
}


void vtkPV3FoamReader::SelectionModified()
{
    vtkDebugMacro(<<"SelectionModified");

    Modified();
}

