/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPV3FoamBlockMeshReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPV3FoamBlockMeshReader - reads a dataset in OpenFOAM bockMesh format
// .SECTION Description
// vtkPV3FoamBlockMeshReader creates an multiblock dataset.
// It uses the OpenFOAM infrastructure (blockMesh).

#ifndef __vtkPV3FoamBlockMeshReader_h
#define __vtkPV3FoamBlockMeshReader_h

// Foam forward declarations
namespace Foam
{
    class vtkPV3FoamBlockMesh;
}

// VTK includes
#include "vtkMultiBlockDataSetAlgorithm.h"

// VTK forward declarations
class vtkDataArraySelection;
class vtkCallbackCommand;


class VTK_IO_EXPORT vtkPV3FoamBlockMeshReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:
    vtkTypeRevisionMacro(vtkPV3FoamBlockMeshReader,vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream&, vtkIndent);

    static vtkPV3FoamBlockMeshReader* New();

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Description:
    // GUI update control
    vtkSetMacro(UpdateGUI, int);
    vtkGetMacro(UpdateGUI, int);

    // Description:
    // FOAM display patch names control
    vtkSetMacro(ShowPointNumbers, int);
    vtkGetMacro(ShowPointNumbers, int);

    // Description:
    // Parts (blocks) selection list control
    vtkDataArraySelection* GetPartSelection();
    int  GetNumberOfPartArrays();
    int  GetPartArrayStatus(const char* name);
    void SetPartArrayStatus(const char* name, int status);
    const char* GetPartArrayName(int index);

    // Description:
    // Callback registered with the SelectionObserver
    // for all the selection lists
    static void SelectionModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );

    void SelectionModified();


protected:

    //- Construct null
    vtkPV3FoamBlockMeshReader();

    //- Destructor
    ~vtkPV3FoamBlockMeshReader();

    //- Return information about mesh, times, etc without loading anything
    virtual int RequestInformation
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    //- Get the mesh/fields for a particular time
    virtual int RequestData
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    //- Fill in additional port information
    virtual int FillOutputPortInformation(int, vtkInformation*);

    // The observer to modify this object when array selections are modified
    vtkCallbackCommand* SelectionObserver;

    char* FileName;

private:

    //- Disallow default bitwise copy construct
    vtkPV3FoamBlockMeshReader(const vtkPV3FoamBlockMeshReader&);

    //- Disallow default bitwise assignment
    void operator=(const vtkPV3FoamBlockMeshReader&);

    //- Add/remove point numbers to/from the view
    void updatePointNumbersView(const bool show);

    int ShowPointNumbers;

    //- Dummy variable/switch to invoke a reader update
    int UpdateGUI;

    vtkDataArraySelection* PartSelection;

    //BTX
    Foam::vtkPV3FoamBlockMesh* foamData_;
    //ETX
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
