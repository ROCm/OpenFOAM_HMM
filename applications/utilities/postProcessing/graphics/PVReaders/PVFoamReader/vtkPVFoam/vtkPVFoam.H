/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::vtkPVFoam

Description
    Provides a reader interface for OpenFOAM to VTK interaction.

SourceFiles
    vtkPVFoam.C
    vtkPVFoamFields.C
    vtkPVFoamMesh.C
    vtkPVFoamMeshLagrangian.C
    vtkPVFoamMeshVolume.C
    vtkPVFoamTemplates.C
    vtkPVFoamUpdateInfo.C
    vtkPVFoamUtils.C
    vtkPVFoamFieldTemplates.C
    vtkPVFoamUpdateTemplates.C

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoam_H
#define vtkPVFoam_H

#include "className.H"
#include "fileName.H"
#include "stringList.H"
#include "wordList.H"
#include "primitivePatch.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "foamPvCore.H"

#undef VTKPVFOAM_DUALPORT

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkFloatArray;
class vtkPoints;
class vtkPVFoamReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkUnstructuredGrid;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// OpenFOAM class forward declarations
class argList;
class Time;
class fvMesh;
class IOobjectList;
class polyPatch;

template<class Type> class IOField;
template<class Type> class Field;
template<class Type> class List;

/*---------------------------------------------------------------------------*\
                        Class vtkPVFoam Declaration
\*---------------------------------------------------------------------------*/

class vtkPVFoam
:
    private foamPvCore
{
    // Convenience typedefs
    typedef PrimitivePatchInterpolation<primitivePatch> patchInterpolator;


    // Private classes

        //- Bookkeeping for polyhedral cell decomposition
        //  hide in extra pointMap (cellSet/cellZone) for now
        class polyDecomp
        {
            labelList superCells_;
            labelList addPointCellLabels_;
            labelList pointMap_;

        public:

            polyDecomp()
            {}

            //- Label of original cell for decomposed cells
            labelList& superCells()
            {
                return superCells_;
            }

            //- Label of original cell for decomposed cells
            const labelList& superCells() const
            {
                return superCells_;
            }

            //- Cell-centre labels for additional points of decomposed cells
            labelList& addPointCellLabels()
            {
                return addPointCellLabels_;
            }

            //- Cell-centre labels for additional points of decomposed cells
            const labelList& addPointCellLabels() const
            {
                return addPointCellLabels_;
            }

            //- Point labels for subsetted meshes
            labelList& pointMap()
            {
                return pointMap_;
            }

            //- Point labels for subsetted meshes
            const labelList& pointMap() const
            {
                return pointMap_;
            }


            //- Clear
            void clear()
            {
                superCells_.clear();
                addPointCellLabels_.clear();
                pointMap_.clear();
            }
        };


    // Private Data

        //- Access to the controlling vtkPVFoamReader
        vtkPVFoamReader* reader_;

        //- OpenFOAM time control
        autoPtr<Time> dbPtr_;

        //- OpenFOAM mesh
        fvMesh* meshPtr_;

        //- The mesh region
        word meshRegion_;

        //- The mesh directory for the region
        fileName meshDir_;

        //- The time index
        int timeIndex_;

        //- Track changes in mesh geometry
        bool meshChanged_;

        //- Track changes in fields
        bool fieldsChanged_;

        //- Selected geometrical parts (internalMesh, patches, ...)
        boolList partStatus_;

        //- Datasets corresponding to selected geometrical pieces
        //  a negative number indicates that no vtkmesh exists for this piece
        labelList partDataset_;

        //- First instance and size of various mesh parts
        //  used to index into partStatus_ and partDataset_
        arrayRange arrayRangeVolume_;
        arrayRange arrayRangePatches_;
        arrayRange arrayRangeLagrangian_;
        arrayRange arrayRangeCellZones_;
        arrayRange arrayRangeFaceZones_;
        arrayRange arrayRangePointZones_;
        arrayRange arrayRangeCellSets_;
        arrayRange arrayRangeFaceSets_;
        arrayRange arrayRangePointSets_;

        //- Decomposed cells information (mesh regions)
        //  TODO: regions
        List<polyDecomp> regionPolyDecomp_;

        //- Decomposed cells information (cellZone meshes)
        List<polyDecomp> zonePolyDecomp_;

        //- Decomposed cells information (cellSet meshes)
        List<polyDecomp> csetPolyDecomp_;

        //- List of patch names for rendering to window
        List<vtkTextActor*> patchTextActorsPtrs_;


    // Private Member Functions

        //- Reset data counters
        void resetCounters();

      // Update information helper functions

        //- Update the mesh parts selected in the GUI
        void updateMeshPartsStatus();

        //- Internal mesh info
        void updateInfoInternalMesh(vtkDataArraySelection*);

        //- Lagrangian info
        void updateInfoLagrangian(vtkDataArraySelection*);

        //- Patch info
        void updateInfoPatches(vtkDataArraySelection*, stringList&);

        //- Set info
        void updateInfoSets(vtkDataArraySelection*);

        //- Zone info
        void updateInfoZones(vtkDataArraySelection*);

        //- Get non-empty zone names for zoneType from file
        wordList getZoneNames(const word& zoneType) const;

        //- Get names of on non-empty zones from the mesh info
        template<class ZoneType>
        static wordList getZoneNames
        (
            const ZoneMesh<ZoneType, polyMesh>& zmesh
        );

        //- Field (volume, point, lagrangian) info
        void updateInfoFields();

        //- Field info
        template<template<class> class patchType, class meshType>
        void updateInfoFields
        (
            vtkDataArraySelection* select
        );

        //- Lagrangian field info
        void updateInfoLagrangianFields(vtkDataArraySelection* select);


      // Update helper functions

        //- OpenFOAM mesh
        void updateFoamMesh();

        //- Reduce memory footprint after conversion
        void reduceMemory();


      // Mesh conversion functions

        //- Convert volume mesh
        void convertMeshVolume(vtkMultiBlockDataSet* output, int& blockNo);

        //- Convert Lagrangian points
        void convertMeshLagrangian(vtkMultiBlockDataSet* output, int& blockNo);

        //- Convert mesh patches
        void convertMeshPatches(vtkMultiBlockDataSet* output, int& blockNo);

        //- Convert cell zones
        void convertMeshCellZones(vtkMultiBlockDataSet* output, int& blockNo);

        //- Convert face zones
        void convertMeshFaceZones(vtkMultiBlockDataSet* output, int& blockNo);

        //- Convert point zones
        void convertMeshPointZones(vtkMultiBlockDataSet* output, int& blockNo);

        //- Convert cell sets
        void convertMeshCellSets(vtkMultiBlockDataSet* output, int& blockNo);

        //- Convert face sets
        void convertMeshFaceSets(vtkMultiBlockDataSet* output, int& blockNo);

        //- Convert point sets
        void convertMeshPointSets(vtkMultiBlockDataSet* output, int& blockNo);


      // Add mesh functions

        //- Volume meshes as vtkUnstructuredGrid
        vtkUnstructuredGrid* volumeVTKMesh
        (
            const fvMesh& mesh,
            polyDecomp& decompInfo
        );

        //- Lagrangian positions as vtkPolyData
        vtkPolyData* lagrangianVTKMesh
        (
            const polyMesh& mesh,
            const word& cloudName
        );

        //- Patches (mesh or primitive) as vtkPolyData
        template<class PatchType>
        vtkPolyData* patchVTKMesh
        (
            const word& name,
            const PatchType& p
        );


      // Field conversion functions

        //- Convert Field to VTK field
        template<class Type>
        vtkFloatArray* convertFieldToVTK
        (
            const word& name,
            const Field<Type>& fld
        );

        //- Face set/zone field
        template<class Type>
        vtkFloatArray* convertFaceFieldToVTK
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const labelUList& faceLabels
        );

        //- Volume field
        template<class Type>
        vtkFloatArray* convertVolFieldToVTK
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const polyDecomp& decompInfo
        );


        //- Convert volume fields
        void convertVolFields(vtkMultiBlockDataSet* output);

        //- Convert point fields
        void convertPointFields(vtkMultiBlockDataSet* output);

        //- Convert Lagrangian fields
        void convertLagrangianFields(vtkMultiBlockDataSet* output);


      // Convert OpenFOAM fields

        //- Volume field - all types
        template<class Type>
        void convertVolField
        (
            const PtrList<patchInterpolator>& patchInterpList,
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            vtkMultiBlockDataSet* output
        );

        //- Volume fields - all types
        template<class Type>
        void convertVolFields
        (
            const fvMesh& mesh,
            const PtrList<patchInterpolator>& patchInterpList,
            const IOobjectList& objects,
            vtkMultiBlockDataSet* output
        );

        //- Volume internal fields (DimensionedField)- all types
        template<class Type>
        void convertDimFields
        (
            const fvMesh& mesh,
            const PtrList<patchInterpolator>& patchInterpList,
            const IOobjectList& objects,
            vtkMultiBlockDataSet* output
        );

        //- Volume field - all selected parts
        template<class Type>
        void convertVolFieldBlock
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr,
            vtkMultiBlockDataSet* output,
            const arrayRange& range,
            const List<polyDecomp>& decompLst
        );

        //- Lagrangian fields - all types
        template<class Type>
        void convertLagrangianFields
        (
            const IOobjectList& objects,
            vtkMultiBlockDataSet* output,
            const label datasetNo
        );

        //- Point fields - all types
        template<class Type>
        void convertPointFields
        (
            const pointMesh& pMesh,
            const IOobjectList& objects,
            vtkMultiBlockDataSet* output
        );

        //- Point field - all selected parts
        template<class Type>
        void convertPointFieldBlock
        (
            const GeometricField<Type, pointPatchField, pointMesh>& pfld,
            vtkMultiBlockDataSet* output,
            const arrayRange& range,
            const List<polyDecomp>& decompLst
        );

        //- Point field
        template<class Type>
        void convertPointField
        (
            vtkUnstructuredGrid* vtkmesh,
            const GeometricField<Type, pointPatchField, pointMesh>& pfld,
            const GeometricField<Type, fvPatchField, volMesh>& vfld,
            const polyDecomp& decomp
        );


      // GUI selection helper functions

        //- Only retain specified fields
        static void pruneObjectList
        (
            IOobjectList& objects,
            const hashedWordList& retain
        );

        //- Get the first word from the reader 'parts' selection
        word getPartName(const int partId);


    // Constructors

        //- Disallow default bitwise copy construct
        vtkPVFoam(const vtkPVFoam&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const vtkPVFoam&) = delete;


public:

    //- Static data members

        ClassName("vtkPVFoam");


    // Constructors

        //- Construct from components
        vtkPVFoam
        (
            const char* const FileName,
            vtkPVFoamReader* reader
        );


    //- Destructor
    ~vtkPVFoam();


    // Member Functions

        //- Update
        void updateInfo();

        void Update
        (
            vtkMultiBlockDataSet* output,
            vtkMultiBlockDataSet* lagrangianOutput
        );

        //- Clean any storage
        void CleanUp();

        //- Return a list of selected times.
        //  Use STL container since these values are used by the plugin
        std::vector<double> findTimes(const bool skipZero = false) const;

        //- Add/remove patch names to/from the view
        void renderPatchNames(vtkRenderer*, const bool show);

        //- Set the runTime to the first plausible request time,
        //  returns the timeIndex
        //  sets to "constant" on error
        int setTime(int count, const double requestTimes[]);


        //- The current time index
        int timeIndex() const
        {
           return timeIndex_;
        }


     // Access

        //- Debug information
        void PrintSelf(ostream&, vtkIndent) const;

        void printInfo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vtkPVFoamTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
