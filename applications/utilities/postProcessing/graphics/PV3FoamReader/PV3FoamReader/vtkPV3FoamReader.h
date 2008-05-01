/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPV3FoamReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPV3FoamReader - reads a dataset in OpenFOAM format
// .SECTION Description
// vtkPV3FoamReader creates an multiblock dataset. It reads a controlDict
// file, mesh information, and time dependent data.  The controlDict file
// contains timestep information. The polyMesh folders contain mesh information
// The time folders contain transient data for the cells  Each folder can
// contain any number of data files.

#ifndef __vtkPV3FoamReader_h
#define __vtkPV3FoamReader_h

#include "vtkMultiBlockDataSetAlgorithm.h"

// Foam forward declarations
namespace Foam
{
    class vtkPV3Foam;
}

// VTK forward declarations
class vtkUnstructuredGrid;
class vtkPoints;
class vtkIntArray;
class vtkFloatArray;
class vtkDoubleArray;
class vtkDataArraySelection;
class vtkCallbackCommand;


class VTK_IO_EXPORT vtkPV3FoamReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:

    static vtkPV3FoamReader* New();

    vtkTypeRevisionMacro
    (
        vtkPV3FoamReader,
        vtkMultiBlockDataSetAlgorithm
    );

    void PrintSelf
    (
        ostream& os,
        vtkIndent indent
    );

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // GUI update control
    vtkSetMacro(UpdateGUI, int);
    vtkGetMacro(UpdateGUI, int);

    // FOAM mesh caching control
    vtkSetMacro(CacheMesh, int);
    vtkGetMacro(CacheMesh, int);

    // FOAM read sets control
    vtkSetMacro(IncludeSets, int);
    vtkGetMacro(IncludeSets, int);

    // FOAM patch names control
    vtkSetMacro(ShowPatchNames, int);
    vtkGetMacro(ShowPatchNames, int);

    // Time-step slider control
    vtkSetMacro(TimeStep, int);
    vtkGetMacro(TimeStep, int);
    vtkSetVector2Macro(TimeStepRange, int);
    vtkGetVector2Macro(TimeStepRange, int);

    // Time selection list control
    vtkDataArraySelection* GetTimeSelection();
    int GetNumberOfTimeArrays();
    const char* GetTimeArrayName(int index);
    int GetTimeArrayStatus(const char* name);
    void SetTimeArrayStatus(const char* name, int status);


    // Region selection list control
    vtkDataArraySelection* GetRegionSelection();
    int GetNumberOfRegionArrays();
    const char* GetRegionArrayName(int index);
    int GetRegionArrayStatus(const char* name);
    void SetRegionArrayStatus(const char* name, int status);

    // volField selection list control
    vtkDataArraySelection* GetVolFieldSelection();
    int GetNumberOfVolFieldArrays();
    const char* GetVolFieldArrayName(int index);
    int GetVolFieldArrayStatus(const char* name);
    void SetVolFieldArrayStatus(const char* name, int status);

    // pointField selection list control
    vtkDataArraySelection* GetPointFieldSelection();
    int GetNumberOfPointFieldArrays();
    const char* GetPointFieldArrayName(int index);
    int GetPointFieldArrayStatus(const char* name);
    void SetPointFieldArrayStatus(const char* name, int status);

    // lagrangianField selection list control
    vtkDataArraySelection* GetLagrangianFieldSelection();
    int GetNumberOfLagrangianFieldArrays();
    const char* GetLagrangianFieldArrayName(int index);
    int GetLagrangianFieldArrayStatus(const char* name);
    void SetLagrangianFieldArrayStatus(const char* name, int status);

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

    vtkPV3FoamReader();
    ~vtkPV3FoamReader();

    char* FileName;

    virtual int RequestData
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    virtual int RequestInformation
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    // The observer to modify this object when the array selections
    // are modified
    vtkCallbackCommand* SelectionObserver;


private:

    vtkPV3FoamReader(const vtkPV3FoamReader&);  // Not implemented.
    void operator=(const vtkPV3FoamReader&);  // Not implemented.

    //- Add patch names to the view
    void addPatchNamesToView();

    //- Remove patch names from the view
    void removePatchNamesFromView();

    int CacheMesh;
    int IncludeSets;
    int ShowPatchNames;

    int UpdateGUI;
    int UpdateGUIOld;
    int TimeStep;
    int TimeStepRange[2];

    vtkDataArraySelection* TimeSelection;
    vtkDataArraySelection* RegionSelection;
    vtkDataArraySelection* VolFieldSelection;
    vtkDataArraySelection* PointFieldSelection;
    vtkDataArraySelection* LagrangianFieldSelection;

    //BTX
    Foam::vtkPV3Foam* foamData_;
    //ETX
};

#endif
