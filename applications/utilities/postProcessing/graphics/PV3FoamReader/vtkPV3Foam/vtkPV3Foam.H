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

Class
    vtkPV3Foam

Description

SourceFiles
    interpolatePointToCell.C
    interpolatePointToCell.H
    vtkPV3Foam.C
    vtkPV3Foam.H
    vtkPV3FoamAddFaceSetMesh.C
    vtkPV3FoamAddFields.H
    vtkPV3FoamAddLagrangianMesh.C
    vtkPV3FoamAddPatchMesh.C
    vtkPV3FoamAddPointSetMesh.C
    vtkPV3FoamAddVolumeMesh.C
    vtkPV3FoamConvertLagrangianFields.H
    vtkPV3FoamConvertMesh.C
    vtkPV3FoamConvertPatchFaceField.H
    vtkPV3FoamConvertPatchPointField.H
    vtkPV3FoamConvertPointFields.H
    vtkPV3FoamConvertVolFields.H
    vtkPV3FoamInsertNextPoint.H
    vtkPV3FoamUpdate.C
    vtkPV3FoamUpdateInformation.C
    vtkPV3FoamUpdateInformationFields.H

    // Needed by VTK?
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3Foam_H
#define vtkPV3Foam_H

#include "className.H"
#include "fileName.H"
#include "volPointInterpolation.H"


// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkMultiBlockDataSet;
class vtkPoints;
class vtkPV3FoamReader;
class vtkRenderer;
class vtkTextActor;
class vtkUnstructuredGrid;

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
class WordList;

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

        //BTX
        //- Indices for datasets in vtkMultiBlockDataSet
        enum
        {
            VOLUME      = 0, // internal mesh AND patches
            LAGRANGIAN  = 1,
            CELLSET     = 2,
            FACESET     = 3,
            POINTSET    = 4
        };
        //ETX

        //- Indices of (first instance of) various regions
        int idRegionVolume_;
        int idRegionLagrangian_;
        int idRegionPatches_;
        int idRegionCellSets_;
        int idRegionFaceSets_;
        int idRegionPointSets_;

        //- Access to the controlling vtkPV3FoamReader
        vtkPV3FoamReader* reader_;

        //- Selected regions
        //  [0] = internal mesh, then lagrangian, patches, and sets
        boolList selectedRegions_;

        //- Selected regions indices in each respective block
        labelList selectedRegionDatasetIds_;

        //- Foam argument handler
        autoPtr<argList> argsPtr_;

        //- Path to root and case
        fileName rootPath_;
        fileName casePath_;

        //- Foam time control
        autoPtr<Time> dbPtr_;

        //- Foam mesh
        fvMesh* meshPtr_;

        //- Labels of cell-centres used as additional points when decomposing
        // polyhedra
        labelList addPointCellLabels_;

        //- Label of original cell the decomposed cells are split from
        //  - internal mesh
        labelList superCells_;

        //- Label of original cell the decomposed cells are split from
        //  - cell set meshes
        List<labelList> superCellSetCells_;

        //- List of patch names
        List<vtkTextActor*> patchTextActorsPtrs_;


        // Dataset sizes

            //- Number of meshes
            // TODO - for info only - only set up to process ONE mesh
            int nMesh_;

            //- Lagrangian data size (i.e. 0 or 1)
            int lagrangianDataSize_;

            //- Cloud name to be processed
            // TODO - currently only set up to process ONE cloud
            word cloudName_;

            //- Patch data size, 0 to n
            int patchDataSize_;

            //- Cell set data size, 0 to n
            int cellSetDataSize_;

            //- Face set data size, 0 to n
            int faceSetDataSize_;

            //- Point set data size, 0 to n
            int pointSetDataSize_;


    // Private Member Functions

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        void AddToBlock
        (
            vtkMultiBlockDataSet* output,
            unsigned int blockNo,
            unsigned int datasetNo,
            vtkDataSet* dataset
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        vtkDataSet* GetDataSetFromBlock
        (
            vtkMultiBlockDataSet* output,
            unsigned int blockNo,
            unsigned int datasetNo
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        label GetNumberOfDataSets
        (
            vtkMultiBlockDataSet* output,
            unsigned int blockNo
        );

        //- Reset data counters
        void resetCounters();

        //- Find and set the selected time from all the methods of selection
        static void setSelectedTime
        (
            Time& runTime,
            vtkPV3FoamReader* reader
        );

        // Update information helper functions

            //- Update the regions selected in the GUI
            void updateSelectedRegions();

            //- Internal mesh info
            void updateInformationInternalMesh();

            //- Lagrangian info
            void updateInformationLagrangian();

            //- Patch info
            void updateInformationPatches();

            //- Set info
            void updateInformationSets();

            //- Field info
            template<template<class> class patchType, class meshType>
            void updateInformationFields
            (
                vtkDataArraySelection *fieldSelection
            );

            //- Lagrangian field info
            void updateInformationLagrangianFields();


        // Update helper functions

            //- Foam mesh
            void updateFoamMesh();

            //- Volume fields
            void updateVolFields(vtkMultiBlockDataSet* output);

            //- Point fields
            void updatePointFields(vtkMultiBlockDataSet* output);

            //- Lagrangian fields
            void updateLagrangianFields(vtkMultiBlockDataSet* output);


        // Mesh conversion functions

            //- Volume mesh
            void convertMeshVolume(vtkMultiBlockDataSet* output);

            //- Lagrangian mesh
            void convertMeshLagrangian(vtkMultiBlockDataSet* output);

            //- Patch meshes
            void convertMeshPatches(vtkMultiBlockDataSet* output);

            //- Cell set meshes
            void convertMeshCellSet(vtkMultiBlockDataSet* output);

            //- Face set meshes
            void convertMeshFaceSet(vtkMultiBlockDataSet* output);

            //- Point set meshes
            void convertMeshPointSet(vtkMultiBlockDataSet* output);


        // Add mesh functions

            //- Add internal mesh/cell set meshes
            void addVolumeMesh
            (
                const fvMesh& mesh,
                vtkUnstructuredGrid* vtkMesh,
                labelList& superCells
            );

            //- Add Lagrangian mesh
            void addLagrangianMesh
            (
                const fvMesh& mesh,
                vtkUnstructuredGrid* vtkMesh
            );

            //- Add patch mesh
            void addPatchMesh
            (
                 const polyPatch& p,
                 vtkUnstructuredGrid* vtkMesh
            );

            //- Add cell set mesh
            void addCellSetMesh
            (
                const fvMesh& mesh,
                vtkUnstructuredGrid* vtkMesh
            );

            //- Add face set mesh
            void addFaceSetMesh
            (
                const fvMesh& mesh,
                const faceSet& fSet,
                vtkUnstructuredGrid* vtkMesh
            );

            //- Add point mesh
            void addPointSetMesh
            (
                const fvMesh& mesh,
                const pointSet& pSet,
                vtkUnstructuredGrid* vtkMesh
            );


        //- Add the fields in the selected time directory to the selection
        //  lists
        template<class GeoField>
        void addFields
        (
            vtkDataArraySelection* fieldSelection,
            const IOobjectList& objects
        );


        // Convert Foam fields

            //- Volume fields - all types
            template<class Type>
            void convertVolFields
            (
                const fvMesh& mesh,
                const volPointInterpolation& pInterp,
                const IOobjectList& objects,
                vtkDataArraySelection* fieldSelection,
                vtkMultiBlockDataSet* output
            );

            //- Volume field - all types except scalar
            template<class Type>
            void convertVolField
            (
                const GeometricField<Type, fvPatchField, volMesh>& tf,
                vtkMultiBlockDataSet* output
            );

            //- Patch field
            template<class Type>
            void convertPatchFaceField
            (
                const word& name,
                const Field<Type>& tf,
                const label regioni,
                vtkMultiBlockDataSet* output
            );

            //- Lagrangian fields - all types
            template<class Type>
            void convertLagrangianFields
            (
                const fvMesh& mesh,
                const IOobjectList& objects,
                vtkDataArraySelection *fieldSelection,
                vtkMultiBlockDataSet* output
            );

            //- Lagrangian field - all types except scalar
            template<class Type>
            void convertLagrangianField
            (
                const IOField<Type>& tf,
                vtkMultiBlockDataSet* output
            );

            //- Point fields - all types
            template<class Type>
            void convertPointFields
            (
                const fvMesh& mesh,
                const IOobjectList& objects,
                vtkDataArraySelection *fieldSelection,
                vtkMultiBlockDataSet* output
            );

            //- Point fields - all types except scalar
            template<class Type>
            void convertPointField
            (
                const GeometricField<Type, pointPatchField, pointMesh>& ptf,
                const GeometricField<Type, fvPatchField, volMesh>& tf,
                vtkMultiBlockDataSet* output
            );

            //- Patch point field - all types except scalar
            template<class Type>
            void convertPatchPointField
            (
                const word& name,
                const Field<Type>& tf,
                const label regioni,
                vtkMultiBlockDataSet* output
            );


       // GUI selection helper functions

            //- Store selection(s)
            wordList getSelectedArrayEntries
            (
                vtkDataArraySelection* arraySelection
            );

            //- Set selection(s)
            void setSelectedArrayEntries
            (
                vtkDataArraySelection* arraySelection,
                const wordList& selectedEntries
            );


        //- Set the name VTK objects
        void SetName(vtkUnstructuredGrid* vtkMesh, const char* name);

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
            vtkPV3FoamReader* reader,
            vtkMultiBlockDataSet* output
        );


    //- Destructor

        ~vtkPV3Foam();


    // Member Functions

        void UpdateInformation();

        void Update(vtkMultiBlockDataSet* output);

        double* timeSteps();

        //- Add patch names to the display
        void addPatchNames(vtkRenderer* renderer);

        //- Remove patch names from the display
        void removePatchNames(vtkRenderer* renderer);


     // Access

        int numberOfTimeSteps();

        int numberOfCells();

        int numberOfPoints();

        int numberOfMeshes() const;

};


// * * * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * //

template<>
void vtkPV3Foam::convertVolField
(
    const GeometricField<scalar, fvPatchField, volMesh>& sf,
    vtkMultiBlockDataSet* output
);


template<>
void vtkPV3Foam::convertPatchFaceField
(
    const word& name,
    const Field<scalar>& sf,
    const label regioni,
    vtkMultiBlockDataSet* output
);


template<>
void vtkPV3Foam::convertPointField
(
    const GeometricField<scalar, pointPatchField, pointMesh>& psf,
    const GeometricField<scalar, fvPatchField, volMesh>& sf,
    vtkMultiBlockDataSet* output
);


template<>
void vtkPV3Foam::convertPatchPointField
(
    const word& name,
    const Field<scalar>& sf,
    const label regioni,
    vtkMultiBlockDataSet* output
);


template<>
void vtkPV3Foam::convertLagrangianField
(
    const IOField<scalar>& sf,
    vtkMultiBlockDataSet* output
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
