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
// vtkPV3FoamReader creates an multiblock dataset.
// It uses the OpenFOAM infrastructure (fvMesh, etc) to
// handle mesh and field data.

#ifndef __vtkPV3FoamReader_h
#define __vtkPV3FoamReader_h

// Foam forward declarations
namespace Foam
{
    class vtkPV3Foam;
}

// VTK includes
#include "vtkMultiBlockDataSetAlgorithm.h"

// VTK forward declarations
class vtkDataArraySelection;
class vtkCallbackCommand;


class VTK_IO_EXPORT vtkPV3FoamReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:
    vtkTypeRevisionMacro(vtkPV3FoamReader,vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream&, vtkIndent);

    static vtkPV3FoamReader* New();

    // Description:
    // Get the current timestep and the timestep range.
    vtkGetVector2Macro(TimeStepRange, int);

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Description:
    // GUI update control
    vtkSetMacro(UpdateGUI, int);
    vtkGetMacro(UpdateGUI, int);

    // Description:
    // FOAM mesh caching control
    vtkSetMacro(CacheMesh, int);
    vtkGetMacro(CacheMesh, int);

    // Description:
    // FOAM extrapolate internal values onto the walls
    vtkSetMacro(ExtrapolateWalls, int);
    vtkGetMacro(ExtrapolateWalls, int);

    // FOAM read sets control
    vtkSetMacro(IncludeSets, int);
    vtkGetMacro(IncludeSets, int);

    // Description:
    // FOAM read zones control
    vtkSetMacro(IncludeZones, int);
    vtkGetMacro(IncludeZones, int);

    // Description:
    // FOAM display patch names control
    vtkSetMacro(ShowPatchNames, int);
    vtkGetMacro(ShowPatchNames, int);

    // Description:
    // Get the current timestep
    int  GetTimeStep();

    // Description:
    // Region selection list control
    vtkDataArraySelection* GetRegionSelection();
    int  GetNumberOfRegionArrays();
    int  GetRegionArrayStatus(const char* name);
    void SetRegionArrayStatus(const char* name, int status);
    const char* GetRegionArrayName(int index);

    // Description:
    // volField selection list control
    vtkDataArraySelection* GetVolFieldSelection();
    int  GetNumberOfVolFieldArrays();
    int  GetVolFieldArrayStatus(const char* name);
    void SetVolFieldArrayStatus(const char* name, int status);
    const char* GetVolFieldArrayName(int index);

    // Description:
    // pointField selection list control
    vtkDataArraySelection* GetPointFieldSelection();
    int  GetNumberOfPointFieldArrays();
    int  GetPointFieldArrayStatus(const char* name);
    void SetPointFieldArrayStatus(const char* name, int status);
    const char* GetPointFieldArrayName(int index);

    // Description:
    // lagrangianField selection list control
    vtkDataArraySelection* GetLagrangianFieldSelection();
    int  GetNumberOfLagrangianFieldArrays();
    int  GetLagrangianFieldArrayStatus(const char* name);
    void SetLagrangianFieldArrayStatus(const char* name, int status);
    const char* GetLagrangianFieldArrayName(int index);

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
    vtkPV3FoamReader();

    //- Destructor
    ~vtkPV3FoamReader();

    //- Return information about mesh, times, etc without loading anything
    virtual int RequestInformation
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    //- Get the mesh/fields for a particular time
    //- Destructor
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
    vtkPV3FoamReader(const vtkPV3FoamReader&);

    //- Disallow default bitwise assignment
    void operator=(const vtkPV3FoamReader&);

    //- Add patch names to the view
    void addPatchNamesToView();

    //- Remove patch names from the view
    void removePatchNamesFromView();

    int TimeStepRange[2];
    int CacheMesh;

    int ExtrapolateWalls;
    int IncludeSets;
    int IncludeZones;
    int ShowPatchNames;

    //- Dummy variable/switch for invoke a reader update
    int UpdateGUI;

    vtkDataArraySelection* RegionSelection;
    vtkDataArraySelection* VolFieldSelection;
    vtkDataArraySelection* PointFieldSelection;
    vtkDataArraySelection* LagrangianFieldSelection;

    //- Cached data for output port0 (experimental!)
    vtkMultiBlockDataSet* output0_;

    //BTX
    Foam::vtkPV3Foam* foamData_;
    //ETX
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
