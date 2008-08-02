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

Class
    Foam::vtkPV3Foam

Description
    Provides a reader interface for OpenFOAM to VTK interaction.

SourceFiles
    vtkPV3Foam.C
    vtkPV3Foam.H
    vtkPV3FoamI.H
    vtkPV3FoamFields.C
    vtkPV3FoamMesh.C
    vtkPV3FoamMeshLagrangian.C
    vtkPV3FoamMeshPatch.C
    vtkPV3FoamMeshSet.C
    vtkPV3FoamMeshVolume.C
    vtkPV3FoamMeshZone.C
    vtkPV3FoamFaceField.H
    vtkPV3FoamLagrangianFields.H
    vtkPV3FoamPatchField.H
    vtkPV3FoamPointFields.H
    vtkPV3FoamPoints.H
    vtkPV3FoamUpdateInfo.C
    vtkPV3FoamUpdateInfoFields.H
    vtkPV3FoamVolFields.H
    vtkPV3FoamAddToSelection.H

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3Foam_H
#define vtkPV3Foam_H

#include "className.H"
#include "fileName.H"
#include "volPointInterpolation.H"
#include "stringList.H"
#include "wordList.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkPoints;
class vtkPV3FoamReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkUnstructuredGrid;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Foam class forward declarations
class argList;
class Time;
class fvMesh;
class IOobjectList;
class polyPatch;
class faceSet;
class pointSet;

template<class Type>
class IOField;

template<class Type>
class List;

/*---------------------------------------------------------------------------*\
                        Class vtkPV3Foam Declaration
\*---------------------------------------------------------------------------*/

class vtkPV3Foam
{
    // Private data

        //- Access to the controlling vtkPV3FoamReader
        vtkPV3FoamReader* reader_;

        //- Foam time control
        autoPtr<Time> dbPtr_;

        //- Foam mesh
        fvMesh* meshPtr_;

        //- Number of meshes
        // TODO - for info only - only set up to process ONE mesh
        int nMesh_;

        //- The time index
        int timeIndex_;

        //- Track changes in mesh geometry
        bool meshChanged_;

        //- Track changes in fields
        bool fieldsChanged_;

        //- Cell-centre labels used as additional points for decomposed cells
        labelList addPointCellLabels_;

        //- Label of original cell for decomposed cells (internal mesh)
        labelList superCells_;

        //- Label of original cell for decomposed cells (cellZone meshes)
        List<labelList> zoneSuperCells_;

        //- Label of original cell for decomposed cells (cellSet meshes)
        List<labelList> csetSuperCells_;

        //- List of patch names
        List<vtkTextActor*> patchTextActorsPtrs_;


        //- Selected geometrical pieces
        //  [0] = internal mesh, patches, zones, sets, lagrangian
        boolList regionStatus_;

        //- Selected regions indices in each respective block
        labelList regionDataset_;


    // Private classes

        //- bookkeeping for the GUI checklists and the multi-block organization
        class selectionInfo
        {
            const char *name_;
            int block_;
            int start_;
            int size_;

        public:

            selectionInfo(const char *name, const int blockNo=0)
            :
                name_(name),
                block_(blockNo),
                start_(-1),
                size_(0)
            {}

            //- Return the block holding these datasets
            int block() const
            {
                return block_;
            }

            int block(int blockNo)
            {
                int prev = block_;
                block_ = blockNo;
                return prev;
            }

            const char* name() const
            {
                return name_;
            }

            int start() const
            {
                return start_;
            }

            int end() const
            {
                return start_ + size_;
            }

            int size() const
            {
                return size_;
            }

            void reset()
            {
                start_ = -1;
                size_ = 0;
            }

            //- Assign new start and reset the size
            void operator=(const int i)
            {
                start_ = i;
                size_ = 0;
            }

            //- Increment the size
            void operator+=(const int n)
            {
                size_ += n;
            }
        };

    // Private data

        //- First instance and size of various regions
        selectionInfo regionInfoVolume_;
        selectionInfo regionInfoPatches_;
        selectionInfo regionInfoLagrangian_;
        selectionInfo regionInfoCellZones_;
        selectionInfo regionInfoFaceZones_;
        selectionInfo regionInfoPointZones_;
        selectionInfo regionInfoCellSets_;
        selectionInfo regionInfoFaceSets_;
        selectionInfo regionInfoPointSets_;


    // Private Member Functions

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static void AddToBlock
        (
            vtkMultiBlockDataSet* output,
            const selectionInfo&,
            const label datasetNo,
            vtkDataSet* dataset,
            const string& blockName=string::null
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static vtkDataSet* GetDataSetFromBlock
        (
            vtkMultiBlockDataSet* output,
            const selectionInfo&,
            const label datasetNo
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static label GetNumberOfDataSets
        (
            vtkMultiBlockDataSet* output,
            const selectionInfo&
        );

        //- Reset data counters
        void resetCounters();

        // Update information helper functions

            //- Update the regions selected in the GUI
            void updateRegionStatus();

            //- Internal mesh info
            void updateInfoInternalMesh();

            //- Lagrangian info
            void updateInfoLagrangian();

            //- Patch info
            void updateInfoPatches();

            //- Set info
            void updateInfoSets();

            //- Zone info
            void updateInfoZones();

            //- Read zone names for zoneType from file
            wordList readZoneNames(const word& zoneType);

            //- Add to paraview array selection
            template<class Type>
            label addToSelection
            (
                vtkDataArraySelection*,
                const IOobjectList&,
                const string& suffix=string::null
            );

            //- Field info
            template<template<class> class patchType, class meshType>
            void updateInfoFields(vtkDataArraySelection*);

            //- Lagrangian field info
            void updateInfoLagrangianFields();


        // Update helper functions

            //- Foam mesh
            void updateFoamMesh();

            //- Volume fields
            void updateVolFields(vtkMultiBlockDataSet*);

            //- Point fields
            void updatePointFields(vtkMultiBlockDataSet*);

            //- Lagrangian fields
            void updateLagrangianFields(vtkMultiBlockDataSet*);


        // Mesh conversion functions

            //- Volume mesh
            void convertMeshVolume(vtkMultiBlockDataSet*, int& blockNo);

            //- Lagrangian mesh
            void convertMeshLagrangian(vtkMultiBlockDataSet*, int& blockNo);

            //- Patch meshes
            void convertMeshPatches(vtkMultiBlockDataSet*, int& blockNo);

            //- Cell zone meshes
            void convertMeshCellZones(vtkMultiBlockDataSet*, int& blockNo);

            //- Face zone meshes
            void convertMeshFaceZones(vtkMultiBlockDataSet*, int& blockNo);

            //- Point zone meshes
            void convertMeshPointZones(vtkMultiBlockDataSet*, int& blockNo);

            //- Cell set meshes
            void convertMeshCellSets(vtkMultiBlockDataSet*, int& blockNo);

            //- Face set meshes
            void convertMeshFaceSets(vtkMultiBlockDataSet*, int& blockNo);

            //- Point set meshes
            void convertMeshPointSets(vtkMultiBlockDataSet*, int& blockNo);


        // Add mesh functions

            //- Add internal mesh/cell set meshes
            vtkUnstructuredGrid* volumeVTKMesh
            (
                const fvMesh&,
                labelList& superCells
            );

            //- Add Lagrangian mesh
            vtkPolyData* lagrangianVTKMesh
            (
                const fvMesh&,
                const word& cloudName
            );

            //- Add patch mesh
            vtkPolyData* patchVTKMesh(const polyPatch&);

            //- Add face zone mesh
            vtkPolyData* faceZoneVTKMesh
            (
                const fvMesh&,
                const labelList& faceLabels
            );

            //- Add point zone
            vtkPolyData* pointZoneVTKMesh
            (
                const fvMesh&,
                const labelList& pointLabels
            );

            //- Add face set mesh
            vtkPolyData* faceSetVTKMesh
            (
                const fvMesh&,
                const faceSet&
            );

            //- Add point mesh
            vtkPolyData* pointSetVTKMesh
            (
                const fvMesh&,
                const pointSet&
            );

        // Field conversion functions

            //- Convert volume fields
            void convertVolFields(vtkMultiBlockDataSet*);

            //- Convert point fields
            void convertPointFields(vtkMultiBlockDataSet*);

            //- Convert Lagrangian fields
            void convertLagrangianFields(vtkMultiBlockDataSet*);


        //- Add the fields in the selected time directory to the selection
        //  lists
        template<class GeoField>
        label addObjectsToSelection
        (
            vtkDataArraySelection*,
            const IOobjectList&,
            const string& suffix=string::null
        );


        // Convert Foam fields

            //- Volume fields - all types
            template<class Type>
            void convertVolFields
            (
                const fvMesh&,
                const volPointInterpolation&,
                const PtrList<PrimitivePatchInterpolation<primitivePatch> >&,
                const IOobjectList&,
                vtkMultiBlockDataSet* output
            );

            //- Volume field
            template<class Type>
            void convertVolField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo,
                labelList& superCells
            );

            //- Patch field
            template<class Type>
            void convertPatchField
            (
                const word& name,
                const Field<Type>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo
            );

            //- face set/zone field
            template<class Type>
            void convertFaceField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo,
                const fvMesh&,
                const labelList& faceLabels
            );

            //- face set/zone field
            template<class Type>
            void convertFaceField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo,
                const fvMesh&,
                const faceSet&
            );

            //- Lagrangian fields - all types
            template<class Type>
            void convertLagrangianFields
            (
                const IOobjectList&,
                vtkMultiBlockDataSet* output,
                const label datasetNo
            );

            //- Lagrangian field
            template<class Type>
            void convertLagrangianField
            (
                const IOField<Type>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo
            );

            //- Point fields - all types
            template<class Type>
            void convertPointFields
            (
                const fvMesh&,
                const pointMesh&,
                const IOobjectList&,
                vtkMultiBlockDataSet* output
            );

            //- Point fields
            template<class Type>
            void convertPointField
            (
                const GeometricField<Type, pointPatchField, pointMesh>&,
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo
            );

            //- Patch point field
            template<class Type>
            void convertPatchPointField
            (
                const word& name,
                const Field<Type>&,
                vtkMultiBlockDataSet* output,
                const selectionInfo&,
                const label datasetNo
            );


       // GUI selection helper functions

            //- Extract up to the first non-word characters
            inline static word getFirstWord(const char*);

            //- Subset an IOobjectList based on a hashSet
            static void subsetObjectList
            (
                IOobjectList&,
                const wordHashSet&
            );

            //- Store the current selection(s)
            static wordHashSet getSelected(vtkDataArraySelection*);

            //- Store the current selection(s)
            static stringList getSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const bool firstWord=false
            );

            //- Store the current selection(s) for a sub-selection
            static stringList getSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const selectionInfo&,
                const bool firstWord=false
            );

            //- Set selection(s)
            static void setSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const stringList&
            );


        //- Disallow default bitwise copy construct
        vtkPV3Foam(const vtkPV3Foam&);

        //- Disallow default bitwise assignment
        void operator=(const vtkPV3Foam&);


public:

    //- Static data members

        ClassName("vtkPV3Foam");


    // Constructors

        //- Construct from components
        vtkPV3Foam
        (
            const char* const FileName,
            vtkPV3FoamReader* reader
        );


    //- Destructor

        ~vtkPV3Foam();


    // Member Functions

        //- Update
        void updateInfo();

        void Update
        (
            vtkMultiBlockDataSet* output,
            vtkMultiBlockDataSet* lagrangianOutput
        );

        //- Allocate and return a list of selected times
        //  returns the count via the parameter
        double* findTimes(int& nTimeSteps);

        //- Add patch names to the display
        void addPatchNames(vtkRenderer* renderer);

        //- Remove patch names from the display
        void removePatchNames(vtkRenderer* renderer);

        //- set the runTime to the requested time, returns the timeIndex
        //  sets to "constant" on error and returns -1
        int setTime(const double& requestedTime);

        //- The current time index
        int timeIndex() const
        {
           return timeIndex_;
        }

     // Access

        //- Debug information
        void PrintSelf(ostream&, vtkIndent) const;

        //- Simple memory used debugging information
        static void printMemory();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#    include "vtkPV3FoamI.H"

#endif

// ************************************************************************* //
