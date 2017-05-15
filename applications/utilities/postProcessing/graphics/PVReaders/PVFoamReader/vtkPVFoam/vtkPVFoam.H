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
#include "foamVtkMeshMaps.H"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkFloatArray;
class vtkPoints;
class vtkPVFoamReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
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
class fvMeshSubset;

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

        //- Bookkeeping for vtkPolyData
        class foamVtpData
        {
        public:
            vtkSmartPointer<vtkPolyData> vtkmesh;
        };


        //- Bookkeeping for vtkUnstructuredGrid
        class foamVtuData : public foamVtkMeshMaps
        {
        public:
            vtkSmartPointer<vtkUnstructuredGrid> vtkmesh;
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
        enum polyMesh::readUpdateState meshState_;

        //- Track changes in fields
        bool fieldsChanged_;

        //- The index of selected parts mapped to their names
        Map<string> selectedPartIds_;

        //- Any information for 2D (VTP) geometries
        HashTable<foamVtpData, string> cachedVtp_;

        //- Cell maps and other information for 3D (VTU) geometries
        HashTable<foamVtuData, string> cachedVtu_;

        //- First instance and size of various mesh parts
        //  used to index into selectedPartIds and thus indirectly into
        //  cachedVtp, cachedVtu
        arrayRange rangeVolume_;
        arrayRange rangePatches_;
        arrayRange rangeLagrangian_;
        arrayRange rangeCellZones_;
        arrayRange rangeFaceZones_;
        arrayRange rangePointZones_;
        arrayRange rangeCellSets_;
        arrayRange rangeFaceSets_;
        arrayRange rangePointSets_;

        //- List of patch names for rendering to window
        List<vtkSmartPointer<vtkTextActor>> patchTextActors_;


    // Private Member Functions

        template<class Container>
        bool addOutputBlock
        (
            vtkMultiBlockDataSet* output,
            const HashTable<Container, string>& cache,
            const arrayRange& selector,
            const bool singleDataset = false
        ) const;


        //- Reset data counters
        void resetCounters();

      // Update information helper functions

        //- Internal mesh info
        void updateInfoInternalMesh(vtkDataArraySelection* select);

        //- Lagrangian info
        void updateInfoLagrangian(vtkDataArraySelection* select);

        //- Patch info, modifies enabledEntries
        void updateInfoPatches
        (
            vtkDataArraySelection* select,
            HashSet<string>& enabledEntries
        );

        //- Set info
        void updateInfoSets(vtkDataArraySelection* select);

        //- Zone info
        void updateInfoZones(vtkDataArraySelection* select);


        //- Get non-empty zone names for zoneType from file
        wordList getZoneNames(const word& zoneType) const;

        //- Get names of on non-empty zones from the mesh info
        template<class ZoneType>
        static wordList getZoneNames
        (
            const ZoneMesh<ZoneType, polyMesh>& zmesh
        );

        //- Field info
        template<template<class> class patchType, class meshType>
        void updateInfoFields
        (
            vtkDataArraySelection* select
        );

        //- Lagrangian field info
        void updateInfoLagrangianFields(vtkDataArraySelection* select);


      // Mesh conversion functions

        //- Convert volume mesh
        void convertMeshVolume();

        //- Convert Lagrangian points
        void convertMeshLagrangian();

        //- Convert mesh patches
        void convertMeshPatches();

        //- Convert subsetted mesh
        void convertMeshSubset
        (
            const fvMeshSubset& subsetter,
            const string& longName
        );

        //- Convert cell zones
        void convertMeshCellZones();

        //- Convert face zones
        void convertMeshFaceZones();

        //- Convert point zones
        void convertMeshPointZones();

        //- Convert cell sets
        void convertMeshCellSets();

        //- Convert face sets
        void convertMeshFaceSets();

        //- Convert point sets
        void convertMeshPointSets();


      // Add mesh functions

        //- Volume meshes as vtkUnstructuredGrid
        vtkSmartPointer<vtkUnstructuredGrid> volumeVTKMesh
        (
            const fvMesh& mesh,
            foamVtuData& vtuData
        );

        //- Lagrangian positions as vtkPolyData
        vtkSmartPointer<vtkPolyData> lagrangianVTKMesh
        (
            const polyMesh& mesh,
            const word& cloudName
        );

        //- Patches (mesh or primitive) as vtkPolyData
        template<class PatchType>
        vtkSmartPointer<vtkPolyData> patchVTKMesh
        (
            const string& name,
            const PatchType& p
        );


      // Field conversion functions

        //- Convert Field to VTK field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertFieldToVTK
        (
            const word& name,
            const Field<Type>& fld
        );

        //- Face set/zone field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertFaceFieldToVTK
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const labelUList& faceLabels
        );

        //- Volume field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertVolFieldToVTK
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const foamVtuData& vtuData
        );


        //- Convert volume fields
        void convertVolFields();

        //- Convert point fields
        void convertPointFields();

        //- Convert Lagrangian fields
        void convertLagrangianFields();


      // Convert OpenFOAM fields

        //- Volume field - all types
        template<class Type>
        void convertVolField
        (
            const PtrList<patchInterpolator>& patchInterpList,
            const GeometricField<Type, fvPatchField, volMesh>& fld
        );

        //- Volume fields - all types
        template<class Type>
        void convertVolFields
        (
            const fvMesh& mesh,
            const PtrList<patchInterpolator>& patchInterpList,
            const IOobjectList& objects
        );

        //- Volume internal fields (DimensionedField)- all types
        template<class Type>
        void convertDimFields
        (
            const fvMesh& mesh,
            const PtrList<patchInterpolator>& patchInterpList,
            const IOobjectList& objects
        );

        //- Volume field - all selected parts
        template<class Type>
        void convertVolFieldBlock
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr,
            const arrayRange& range
        );

        //- Lagrangian fields - all types
        template<class Type>
        void convertLagrangianFields
        (
            const IOobjectList& objects,
            vtkPolyData* vtkmesh
        );

        //- Point fields - all types
        template<class Type>
        void convertPointFields
        (
            const pointMesh& pMesh,
            const IOobjectList& objectst
        );

        //- Point field - all selected parts
        template<class Type>
        void convertPointFieldBlock
        (
            const GeometricField<Type, pointPatchField, pointMesh>& pfld,
            const arrayRange& range
        );

        //- Point field
        template<class Type>
        void convertPointField
        (
            vtkUnstructuredGrid* vtkmesh,
            const GeometricField<Type, pointPatchField, pointMesh>& pfld,
            const GeometricField<Type, fvPatchField, volMesh>& vfld,
            const foamVtuData& vtuData
        );


      // GUI selection helper functions

        //- Only retain specified fields
        static void pruneObjectList
        (
            IOobjectList& objects,
            const hashedWordList& retain
        );


        //- Get the first word from selectedPartIds_
        word getPartName(const int partId) const;

        //- Get the first word from the reader 'parts' selection
        word getReaderPartName(const int partId) const;


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
            vtkMultiBlockDataSet* outputLagrangian
        );

        //- Final part of Update(), after any last minute rendering.
        void UpdateFinalize();

        //- Add/remove patch names to/from the view
        void renderPatchNames(vtkRenderer* renderer, const bool show);

        //- Return a list of selected times.
        //  Use STL container since these values are used by the plugin
        std::vector<double> findTimes(const bool skipZero = false) const;

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
