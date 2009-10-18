/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPV3blockMeshReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPV3blockMeshReader.h"

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
#include "vtkPV3blockMesh.H"

vtkCxxRevisionMacro(vtkPV3blockMeshReader, "$Revision: 1.5$");
vtkStandardNewMacro(vtkPV3blockMeshReader);

vtkPV3blockMeshReader::vtkPV3blockMeshReader()
{
    Debug = 0;
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName  = NULL;
    foamData_ = NULL;

    ShowPointNumbers = 1;
    UpdateGUI = 0;

    PartSelection = vtkDataArraySelection::New();
    CurvedEdgesSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPV3blockMeshReader::SelectionModifiedCallback
    );
    SelectionObserver->SetClientData(this);


    PartSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );

    CurvedEdgesSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
}


vtkPV3blockMeshReader::~vtkPV3blockMeshReader()
{
    vtkDebugMacro(<<"Deconstructor");

    if (foamData_)
    {
        // remove point numbers
        updatePointNumbersView(false);
        delete foamData_;
    }

    if (FileName)
    {
        delete [] FileName;
    }

    PartSelection->RemoveObserver(this->SelectionObserver);
    CurvedEdgesSelection->RemoveObserver(this->SelectionObserver);

    SelectionObserver->Delete();
    PartSelection->Delete();
}


// Do everything except set the output info
int vtkPV3blockMeshReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestInformation");

    if (Foam::vtkPV3blockMesh::debug)
    {
        cout<<"REQUEST_INFORMATION\n";
    }

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPV3blockMesh::debug)
    {
        cout<<"RequestInformation with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    if (!foamData_)
    {
        foamData_ = new Foam::vtkPV3blockMesh(FileName, this);
    }
    else
    {
        foamData_->updateInfo();
    }

    // might need some other type of error handling

//    {
//        vtkErrorMacro("could not find valid OpenFOAM blockMesh");
//
//        // delete foamData and flag it as fatal error
//        delete foamData_;
//        foamData_ = NULL;
//        return 0;
//    }


    return 1;
}


// Set the output info
int vtkPV3blockMeshReader::RequestData
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

    // catch previous error
    if (!foamData_)
    {
        vtkErrorMacro("Reader failed - perhaps no mesh?");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPV3blockMesh::debug)
    {
        cout<<"RequestData with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
    (
        outputVector->GetInformationObject(0)->Get
        (
            vtkMultiBlockDataSet::DATA_OBJECT()
        )
    );

    if (Foam::vtkPV3blockMesh::debug)
    {
        cout<< "update output with "
            << output->GetNumberOfBlocks() << " blocks\n";
    }


    foamData_->Update(output);
    updatePointNumbersView(ShowPointNumbers);

    // Do any cleanup on the Foam side
    foamData_->CleanUp();

    return 1;
}


void vtkPV3blockMeshReader::updatePointNumbersView(const bool show)
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();

    // Get all the pqRenderView instances
    QList<pqRenderView*> renderViews = smModel->findItems<pqRenderView*>();

    for (int viewI=0; viewI<renderViews.size(); ++viewI)
    {
        foamData_->renderPointNumbers
        (
            renderViews[viewI]->getRenderViewProxy()->GetRenderer(),
            show
        );
    }
}


void vtkPV3blockMeshReader::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os<< indent << "File name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";

    foamData_->PrintSelf(os, indent);
}


// ----------------------------------------------------------------------
// Parts selection list control

vtkDataArraySelection* vtkPV3blockMeshReader::GetPartSelection()
{
    vtkDebugMacro(<<"GetPartSelection");
    return PartSelection;
}


int vtkPV3blockMeshReader::GetNumberOfPartArrays()
{
    vtkDebugMacro(<<"GetNumberOfPartArrays");
    return PartSelection->GetNumberOfArrays();
}


const char* vtkPV3blockMeshReader::GetPartArrayName(int index)
{
    vtkDebugMacro(<<"GetPartArrayName");
    return PartSelection->GetArrayName(index);
}


int vtkPV3blockMeshReader::GetPartArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetPartArrayStatus");
    return PartSelection->ArrayIsEnabled(name);
}


void vtkPV3blockMeshReader::SetPartArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetPartArrayStatus");
    if (status)
    {
        PartSelection->EnableArray(name);
    }
    else
    {
        PartSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// CurvedEdges selection list control

vtkDataArraySelection* vtkPV3blockMeshReader::GetCurvedEdgesSelection()
{
    vtkDebugMacro(<<"GetCurvedEdgesSelection");
    return CurvedEdgesSelection;
}


int vtkPV3blockMeshReader::GetNumberOfCurvedEdgesArrays()
{
    vtkDebugMacro(<<"GetNumberOfCurvedEdgesArrays");
    return CurvedEdgesSelection->GetNumberOfArrays();
}


const char* vtkPV3blockMeshReader::GetCurvedEdgesArrayName(int index)
{
    vtkDebugMacro(<<"GetCurvedEdgesArrayName");
    return CurvedEdgesSelection->GetArrayName(index);
}


int vtkPV3blockMeshReader::GetCurvedEdgesArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetCurvedEdgesArrayStatus");
    return CurvedEdgesSelection->ArrayIsEnabled(name);
}


void vtkPV3blockMeshReader::SetCurvedEdgesArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetCurvedEdgesArrayStatus");
    if (status)
    {
        CurvedEdgesSelection->EnableArray(name);
    }
    else
    {
        CurvedEdgesSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------

void vtkPV3blockMeshReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkPV3blockMeshReader*>(clientdata)->SelectionModified();
}


void vtkPV3blockMeshReader::SelectionModified()
{
    vtkDebugMacro(<<"SelectionModified");
    Modified();
}


int vtkPV3blockMeshReader::FillOutputPortInformation
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
