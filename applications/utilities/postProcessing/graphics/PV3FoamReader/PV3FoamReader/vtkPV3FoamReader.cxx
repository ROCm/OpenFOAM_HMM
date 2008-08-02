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
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkSMRenderViewProxy.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"

// Foam includes
#include "vtkPV3Foam.H"

vtkCxxRevisionMacro(vtkPV3FoamReader, "$Revision: 1.5$");
vtkStandardNewMacro(vtkPV3FoamReader);

#undef EXPERIMENTAL_TIME_CACHING

vtkPV3FoamReader::vtkPV3FoamReader()
{
    Debug = 0;
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName  = NULL;
    foamData_ = NULL;

    output0_  = NULL;

    // Add second output for the Lagrangian
    this->SetNumberOfOutputPorts(2);
    vtkMultiBlockDataSet *lagrangian;
    lagrangian = vtkMultiBlockDataSet::New();
    lagrangian->ReleaseData();

    this->GetExecutive()->SetOutputData(1, lagrangian);
    lagrangian->Delete();

    TimeStepRange[0] = 0;
    TimeStepRange[1] = 0;

    CacheMesh = 1;

    ExtrapolateWalls = 0;
    IncludeSets = 0;
    IncludeZones = 0;
    ShowPatchNames = 0;

    UpdateGUI = 0;

    RegionSelection = vtkDataArraySelection::New();
    VolFieldSelection = vtkDataArraySelection::New();
    PointFieldSelection = vtkDataArraySelection::New();
    LagrangianFieldSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPV3FoamReader::SelectionModifiedCallback
    );
    SelectionObserver->SetClientData(this);

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

    if (output0_)
    {
        output0_->Delete();
    }


    RegionSelection->RemoveObserver(this->SelectionObserver);
    VolFieldSelection->RemoveObserver(this->SelectionObserver);
    PointFieldSelection->RemoveObserver(this->SelectionObserver);
    LagrangianFieldSelection->RemoveObserver(this->SelectionObserver);

    SelectionObserver->Delete();

    RegionSelection->Delete();
    VolFieldSelection->Delete();
    PointFieldSelection->Delete();
    LagrangianFieldSelection->Delete();
}


// Do everything except set the output info
int vtkPV3FoamReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestInformation");

    if (Foam::vtkPV3Foam::debug)
    {
        cout<<"REQUEST_INFORMATION\n";
    }

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPV3Foam::debug)
    {
        cout<<"RequestInformation with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    if (!foamData_)
    {
        foamData_ = new Foam::vtkPV3Foam(FileName, this);
    }
    else
    {
        foamData_->updateInfo();
    }

    int nTimeSteps = 0;
    double* timeSteps = foamData_->findTimes(nTimeSteps);


    // set identical time steps for all ports
    for (int infoI = 0; infoI < nInfo; ++infoI)
    {
        outputVector->GetInformationObject(infoI)->Set
        (
            vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
            timeSteps,
            nTimeSteps
        );
    }

    double timeRange[2];
    if (nTimeSteps)
    {
        timeRange[0] = timeSteps[0];
        timeRange[1] = timeSteps[nTimeSteps-1];

        if (Foam::vtkPV3Foam::debug > 1)
        {
            cout<<"nTimeSteps " << nTimeSteps << "\n"
                <<"timeRange " << timeRange[0] << " to " << timeRange[1]
                << "\n";

            for (int timeI = 0; timeI < nTimeSteps; ++timeI)
            {
                cout<< "step[" << timeI << "] = " << timeSteps[timeI] << "\n";
            }
        }

        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Set
            (
                vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
                timeRange,
                2
            );
        }
    }

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

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPV3Foam::debug)
    {
        cout<<"RequestData with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    // take port0 as the lead for other outputs
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
    (
        outInfo->Get
        (
            vtkMultiBlockDataSet::DATA_OBJECT()
        )
    );


    vtkMultiBlockDataSet* lagrangianOutput = vtkMultiBlockDataSet::SafeDownCast
    (
        outputVector->GetInformationObject(1)->Get
        (
            vtkMultiBlockDataSet::DATA_OBJECT()
        )
    );

    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
        // Get the requested time step.
        // We only support requests for a single time step
        int nRequestedTimeSteps = outInfo->Length
        (
            vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()
        );
        if (nRequestedTimeSteps >= 1)
        {
            double *requestedTimeSteps = outInfo->Get
            (
                vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()
            );

            foamData_->setTime(requestedTimeSteps[0]);
        }
    }

    if (Foam::vtkPV3Foam::debug)
    {
        cout<< "update output with "
            << output->GetNumberOfBlocks() << " blocks\n";
    }

#ifdef EXPERIMENTAL_TIME_CACHING
    bool needsUpdate = false;

    if (!output0_)
    {
        output0_ = vtkMultiBlockDataSet::New();
        needsUpdate = true;
    }

    // This experimental bit of code seems to work for the geometry,
    // but trashes the fields and still triggers the GeometryFilter
    if (needsUpdate)
    {
        foamData_->Update(output);
        output0_->ShallowCopy(output);
    }
    else
    {
        output->ShallowCopy(output0_);
    }

    if (Foam::vtkPV3Foam::debug)
    {
        if (needsUpdate)
        {
            cout<< "full UPDATE ---------\n";
        }
        else
        {
            cout<< "cached UPDATE ---------\n";
        }

        cout<< "UPDATED output: ";
        output->Print(cout);

        cout<< "UPDATED output0_: ";
        output0_->Print(cout);
    }

#else

    foamData_->Update(output, lagrangianOutput);

    if (ShowPatchNames)
    {
        addPatchNamesToView();
    }
    else
    {
        removePatchNamesFromView();
    }

#endif

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


void vtkPV3FoamReader::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os<< indent << "File name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";

    foamData_->PrintSelf(os, indent);

    os<< indent << "Time step range: "
      << this->TimeStepRange[0] << " - " << this->TimeStepRange[1]
      << "\n";
    os<< indent << "Time step: " << this->GetTimeStep() << endl;
}


int vtkPV3FoamReader::GetTimeStep()
{
    return foamData_ ? foamData_->timeIndex() : -1;
}


// ----------------------------------------------------------------------
// Region selection list control

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


const char* vtkPV3FoamReader::GetRegionArrayName(int index)
{
    vtkDebugMacro(<<"GetRegionArrayName");
    return RegionSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetRegionArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetRegionArrayStatus");
    return RegionSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetRegionArrayStatus(const char* name, int status)
{
    vtkDebugMacro(<<"SetRegionArrayStatus");
    if (status)
    {
        RegionSelection->EnableArray(name);
    }
    else
    {
        RegionSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// volField selection list control

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


const char* vtkPV3FoamReader::GetVolFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetVolFieldArrayName");
    return VolFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetVolFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetVolFieldArrayStatus");
    return VolFieldSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetVolFieldArrayStatus(const char* name, int status)
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


// ----------------------------------------------------------------------
// pointField selection list control

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


const char* vtkPV3FoamReader::GetPointFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetPointFieldArrayName");
    return PointFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetPointFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetPointFieldArrayStatus");
    return PointFieldSelection->ArrayIsEnabled(name);
}


void vtkPV3FoamReader::SetPointFieldArrayStatus(const char* name, int status)
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


// ----------------------------------------------------------------------
// lagrangianField selection list control

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


const char* vtkPV3FoamReader::GetLagrangianFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayName");
    return LagrangianFieldSelection->GetArrayName(index);
}


int vtkPV3FoamReader::GetLagrangianFieldArrayStatus(const char* name)
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

// ----------------------------------------------------------------------

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


int vtkPV3FoamReader::FillOutputPortInformation
(
    int port,
    vtkInformation* info
)
{
    if (port == 0)
    {
        return this->Superclass::FillOutputPortInformation(port, info);
    }
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
}

// ************************************************************************* //
