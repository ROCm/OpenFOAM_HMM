/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::faMesh

Description
    Finite area mesh.  Used for 2-D non-Euclidian finite area method.

SourceFiles
    faMesh.C
    faMeshDemandDrivenData.C
    faMeshPatches.C
    faMeshTopology.C
    faMeshUpdate.C

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faMesh_H
#define faMesh_H

#include "MeshObject.H"
#include "polyMesh.H"
#include "lduMesh.H"
#include "faBoundaryMesh.H"
#include "edgeList.H"
#include "faceList.H"
#include "primitiveFieldsFwd.H"
#include "DimensionedField.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "uindirectPrimitivePatch.H"
#include "edgeInterpolation.H"
#include "labelIOList.H"
#include "FieldFields.H"
#include "faGlobalMeshData.H"
#include "faSchemes.H"
#include "faSolution.H"
#include "data.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMeshBoundaryHalo;
class faMeshLduAddressing;
class faMeshMapper;
class faPatchData;

/*---------------------------------------------------------------------------*\
                           Class faMesh Declaration
\*---------------------------------------------------------------------------*/

class faMesh
:
    public MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>,
    public lduMesh,
    public edgeInterpolation,
    public faSchemes,
    public faSolution,
    public data
{
    // Private (internal) classes/structures

        //- A (proc, patchi, patchEdgei) tuple used internally for managing
        //- patch/patch bookkeeping during construction.
        //  Finite-area patches are stored with negated indices, which makes
        //  them readily identifiable and always sort before normal patches.
        //  Note
        struct patchTuple
        :
            public FixedList<label, 4>
        {

        // Constructors

            // Inherit constructors
            using FixedList<label, 4>::FixedList;

            //- Default construct as 'invalid'
            patchTuple()
            {
                clear();
            }


        // Static Member Functions

            // Globally consistent ordering
            //
            // 1. sort left/right as lower/higher processor connection
            // 2. sort by proc/patch/patch index
            static void sort(UList<Pair<patchTuple>>& list)
            {
                for (auto& tuples : list)
                {
                    tuples.sort();
                }
                Foam::stableSort(list);
            }


        // Member Functions

            //- Reset to 'invalid'
            void clear()
            {
                procNo(-1);
                patchi(labelMax);
                patchEdgei(-1);
                meshFacei(-1);
            }

            //- Valid if proc and edge are non-negative
            bool valid() const noexcept
            {
                return (procNo() >= 0 && patchEdgei() >= 0);
            }

            // Processor is the first sort index
            label procNo() const { return (*this)[0]; }
            void procNo(label val) { (*this)[0] = val; }

            // PatchId (-ve for finiteArea patches) is the second sort index
            label patchi() const { return (*this)[1]; }
            void patchi(label val) { (*this)[1] = val; }

            // The patch edge index (on the finiteArea patch)
            // is the third sort index
            label patchEdgei() const { return (*this)[2]; }
            void patchEdgei(label val) { (*this)[2] = val; }

            // The processor-local mesh face is the fourth sort index
            label meshFacei() const { return (*this)[3]; }
            void meshFacei(label val) { (*this)[3] = val; }

            //- Return the real patchId
            label realPatchi() const
            {
                const label id = patchi();
                return (id < 0 ? -(id + 1) : id);
            }

            //- Set patchId as finiteArea
            void faPatchi(label val)
            {
                patchi(-(val + 1));
            }

            //- Considered to be finiteArea if patchi < 0
            bool is_finiteArea() const noexcept
            {
                return (patchi() < 0);
            }

            //- Considered to be processor local
            bool is_localProc() const noexcept
            {
                return (procNo() == Pstream::myProcNo());
            }
        };


    // Private Data

        //- Face labels
        labelIOList faceLabels_;

        //- Boundary mesh
        faBoundaryMesh boundary_;


        // Primitive mesh data

            //- Edges, addressing into local point list
            edgeList edges_;

            //- Edge owner
            labelList edgeOwner_;

            //- Edge neighbour
            labelList edgeNeighbour_;


        // Primitive size data

            //- Number of points
            mutable label nPoints_;

            //- Number of edges
            mutable label nEdges_;

            //- Number of internal edges
            mutable label nInternalEdges_;

            //- Number of faces
            mutable label nFaces_;


        // Communication support

            //- Communicator used for parallel communication
            label comm_;


    // Demand-driven data

        //- Primitive patch
        mutable std::unique_ptr<uindirectPrimitivePatch> patchPtr_;

        //- List of proc/mesh-face for boundary edge neighbours
        mutable std::unique_ptr<List<labelPair>> bndConnectPtr_;

        //- Ldu addressing data
        mutable faMeshLduAddressing* lduPtr_;

        //- Current time index for motion
        //  Note.  The whole mechanism will be replaced once the
        //  dimensionedField is created and the dimensionedField
        //  will take care of the old-time levels.
        mutable label curTimeIndex_;

        //- Face areas
        mutable DimensionedField<scalar, areaMesh>* SPtr_;

        //- Face areas old time level
        mutable DimensionedField<scalar, areaMesh>* S0Ptr_;

        //- Face areas old-old time level
        mutable DimensionedField<scalar, areaMesh>* S00Ptr_;

        //- Patch starts in the edge list
        mutable labelList* patchStartsPtr_;

        //- Edge length vectors
        mutable edgeVectorField* LePtr_;

        //- Mag edge length vectors
        mutable edgeScalarField* magLePtr_;

        //- Face centres
        mutable areaVectorField* centresPtr_;

        //- Edge centres
        mutable edgeVectorField* edgeCentresPtr_;

        //- Face area normals
        mutable areaVectorField* faceAreaNormalsPtr_;

        //- Edge area normals
        mutable edgeVectorField* edgeAreaNormalsPtr_;

        //- Point area normals
        mutable std::unique_ptr<vectorField> pointAreaNormalsPtr_;

        //- Face curvatures
        mutable areaScalarField* faceCurvaturesPtr_;

        //- Edge transformation tensors
        mutable FieldField<Field, tensor>* edgeTransformTensorsPtr_;

        //- Whether point normals must be corrected for a patch
        mutable boolList* correctPatchPointNormalsPtr_;


    // Other mesh-related data

        //- Parallel info
        mutable autoPtr<faGlobalMeshData> globalMeshDataPtr_;

        //- Mapping/swapping for boundary edge halo neighbours
        mutable std::unique_ptr<faMeshBoundaryHalo> haloMapPtr_;

        //- Face centres for boundary edge halo neighbours
        mutable std::unique_ptr<pointField> haloFaceCentresPtr_;

        //- Face normals for boundary edge halo neighbours
        mutable std::unique_ptr<vectorField> haloFaceNormalsPtr_;


    // Static Private Data

        //- Use the original method for point normals
        static int origPointAreaMethod_;

        //- Use quadrics fit
        static const int quadricsFit_;


    // Private Member Functions

        //- No copy construct
        faMesh(const faMesh&) = delete;

        //- No copy assignment
        void operator=(const faMesh&) = delete;

        //- Set indirect patch, removing any old one
        void initPatch() const;

        //- Set primitive mesh data
        void setPrimitiveMeshData();

        //- Get list of (proc/patchi/patchEdgei/meshFacei) tuple pairs in an
        //- globally consistent ordering
        List<Pair<patchTuple>> getBoundaryEdgeConnections() const;

        //- Determine the boundary edge neighbour connections
        void calcBoundaryConnections() const;

        //- Define boundary edge neighbours (proc/face) based on
        //- gathered topology information
        void setBoundaryConnections
        (
            const List<Pair<patchTuple>>& bndEdgeConnections
        ) const;


        // Private member functions to calculate demand driven data

            //- Calculate ldu addressing
            void calcLduAddressing() const;

            //- Calculate patch starts in the edge list
            void calcPatchStarts() const;

            //- Calculate edge lengths
            void calcLe() const;

            //- Calculate mag edge lengths
            void calcMagLe() const;

            //- Calculate face centres
            void calcAreaCentres() const;

            //- Calculate edge centres
            void calcEdgeCentres() const;

            //- Calculate face areas
            void calcS() const;

            //- Calculate face area normals
            void calcFaceAreaNormals() const;

            //- Calculate edge area normals
            void calcEdgeAreaNormals() const;

            //- Calculate point area normals
            void calcPointAreaNormals_orig(vectorField& result) const;

            //- Calculate point area normals
            void calcPointAreaNormals(vectorField& result) const;

            //- Calculate point area normals by quadrics fit
            void calcPointAreaNormalsByQuadricsFit(vectorField& result) const;

            //- Calculate face curvatures
            void calcFaceCurvatures() const;

            //- Calculate edge transformation tensors
            void calcEdgeTransformTensors() const;

            //- Clear geometry but not the face areas
            void clearGeomNotAreas() const;

            //- Clear boundary halo information
            void clearHalo() const;

            //- Clear geometry
            void clearGeom() const;

            //- Clear addressing
            void clearAddressing() const;

            //- Clear demand-driven data
            void clearOut() const;


    // Halo handling

        //- Calculate halo centres/normals
        void calcHaloFaceGeometry() const;


    // Helpers

        //- Create a single patch
        PtrList<faPatch> createOnePatch
        (
            const word& patchName,
            const word& patchType = ""
        ) const;

        //- Create list of patches from boundary definition
        PtrList<faPatch> createPatchList
        (
            const dictionary& bndDict,
            const word& emptyPatchName = "",
            const dictionary* defaultPatchDefinition = nullptr
        ) const;


        //- Fatal error if edge labels are out of range
        void checkBoundaryEdgeLabelRange(const labelUList& edgeLabels) const;

        //- Extract list from contiguous (unordered) boundary data
        //- to the locally sorted order.
        template<class T>
        List<T> boundarySubset
        (
            const UList<T>& bndField,
            const labelUList& edgeLabels
        ) const
        {
            #ifdef FULLDEBUG
            checkBoundaryEdgeLabelRange(edgeLabels);
            #endif
            // Like UIndirectList but with an offset
            List<T> result(edgeLabels.size());
            forAll(edgeLabels, i)
            {
                result[i] = bndField[edgeLabels[i] - nInternalEdges_];
            }
            return result;
        }


public:

    // Public Typedefs

        typedef faMesh Mesh;
        typedef faBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("faMesh");

        //- The prefix to local: %finite-area
        static const word prefix;

        //- The mesh sub-directory name (usually "faMesh")
        static word meshSubDir;


    // Constructors

        //- Construct zero-sized from polyMesh
        //  Boundary is added using addFaPatches() member function
        faMesh(const polyMesh& pMesh, const Foam::zero);

        //- Construct from polyMesh
        explicit faMesh(const polyMesh& pMesh);

        //- Construct for specified face labels without boundary.
        //  Boundary is added using addFaPatches() member function
        faMesh
        (
            const polyMesh& pMesh,
            const UList<label>& faceLabels
        );

        //- Construct from single polyPatch
        explicit faMesh(const polyPatch& pp);

        //- Construct from definition
        faMesh
        (
            const polyMesh& pMesh,
            const dictionary& faMeshDefinition
        );


    //- Destructor
    virtual ~faMesh();


    // Member Functions

    // Helpers

        //- Add boundary patches. Constructor helper
        void addFaPatches
        (
            PtrList<faPatch>& plist,
            const bool validBoundary = true
        );

        //- Add boundary patches. Constructor helper
        void addFaPatches
        (
            const List<faPatch*>& p,
            const bool validBoundary = true
        );


        // Database

            //- Return access to polyMesh
            inline const polyMesh& mesh() const;

            //- Interface to referenced polyMesh (similar to GeoMesh)
            const polyMesh& operator()() const { return mesh(); }

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshDir() const;

            //- Return reference to time
            const Time& time() const;

            //- Return the current instance directory for points
            //  Used in the construction of geometric mesh data dependent
            //  on points
            const fileName& pointsInstance() const;

            //- Return the current instance directory for faces
            const fileName& facesInstance() const;


        // Communication support

            //- Return communicator used for parallel communication
            inline label comm() const noexcept;

            //- Return communicator used for parallel communication
            inline label& comm() noexcept;


            // Mesh size parameters

                //- Number of local mesh points
                inline label nPoints() const noexcept;

                //- Number of local mesh edges
                inline label nEdges() const noexcept;

                //- Number of internal faces
                inline label nInternalEdges() const noexcept;

                //- Number of boundary edges (== nEdges - nInternalEdges)
                inline label nBoundaryEdges() const noexcept;

                //- Number of patch faces
                inline label nFaces() const noexcept;


            // Primitive mesh data

                //- Return local patch points
                inline const pointField& points() const;

                //- Return local patch edges with reordered boundary
                inline const edgeList& edges() const noexcept;

                //- Return local patch faces
                inline const faceList& faces() const;

                //- Edge owner addressing
                inline const labelList& edgeOwner() const noexcept;

                //- Edge neighbour addressing
                inline const labelList& edgeNeighbour() const noexcept;



        // Access

            //- Return true if thisDb() is a valid DB
            virtual bool hasDb() const;

            //- Return reference to the mesh database
            virtual const objectRegistry& thisDb() const;

            //- Name function is needed to disambiguate those inherited
            //  from base classes
            const word& name() const
            {
                return thisDb().name();
            }

            //- Return constant reference to boundary mesh
            inline const faBoundaryMesh& boundary() const noexcept;

            //- Return the underlying polyMesh face labels
            inline const labelList& faceLabels() const noexcept;


            //- Return parallel info
            const faGlobalMeshData& globalData() const;

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return boundary().interfaces();
            }

            //- Internal face owner
            const labelUList& owner() const
            {
                return lduAddr().lowerAddr();
            }

            //- Internal face neighbour
            const labelUList& neighbour() const
            {
                return lduAddr().upperAddr();
            }

            //- True if given edge label is internal to the mesh
            bool isInternalEdge(const label edgeIndex) const
            {
                return (edgeIndex < nInternalEdges_);
            }

            //- List of proc/face for the boundary edge neighbours
            //- using primitive patch edge numbering.
            inline const List<labelPair>& boundaryConnections() const;

            //- Boundary edge neighbour processors
            //- (does not include own proc)
            labelList boundaryProcs() const;

            //- List of proc/size for the boundary edge neighbour processors
            //- (does not include own proc)
            List<labelPair> boundaryProcSizes() const;

            //- Mapping/swapping for boundary halo neighbours
            const faMeshBoundaryHalo& boundaryHaloMap() const;

            //- Face centres of boundary halo neighbours
            const pointField& haloFaceCentres() const;

            //- Face normals of boundary halo neighbours
            const vectorField& haloFaceNormals() const;

            //- Face centres of boundary halo neighbours for specified patch
            tmp<pointField> haloFaceCentres(const label patchi) const;

            //- Face normals of boundary halo neighbours for specified patch
            tmp<vectorField> haloFaceNormals(const label patchi) const;


        // Mesh motion and morphing

            //- Is mesh moving
            bool moving() const
            {
                return mesh().moving();
            }

            //- Update after mesh motion
            virtual bool movePoints();

            //- Update after topo change
            virtual void updateMesh(const mapPolyMesh&);


        // Mapping

            //- Map all fields in time using given map.
            virtual void mapFields(const faMeshMapper& mapper) const;

            //- Map face areas in time using given map.
            virtual void mapOldAreas(const faMeshMapper& mapper) const;


        // Demand-driven data

            //- Return constant reference to primitive patch
            inline const uindirectPrimitivePatch& patch() const;

            //- Return reference to primitive patch
            inline uindirectPrimitivePatch& patch();

            //- Return patch starts
            const labelList& patchStarts() const;

            //- Return edge length vectors
            const edgeVectorField& Le() const;

            //- Return edge length magnitudes
            const edgeScalarField& magLe() const;

            //- Return face centres as areaVectorField
            const areaVectorField& areaCentres() const;

            //- Return edge centres as edgeVectorField
            const edgeVectorField& edgeCentres() const;

            //- Return face areas
            const DimensionedField<scalar, areaMesh>& S() const;

            //- Return old-time face areas
            const DimensionedField<scalar, areaMesh>& S0() const;

            //- Return old-old-time face areas
            const DimensionedField<scalar, areaMesh>& S00() const;

            //- Return face area normals
            const areaVectorField& faceAreaNormals() const;

            //- Return edge area normals
            const edgeVectorField& edgeAreaNormals() const;

            //- Return point area normals
            const vectorField& pointAreaNormals() const;

            //- Return face curvatures
            const areaScalarField& faceCurvatures() const;

            //- Return edge transformation tensors
            const FieldField<Field, tensor>& edgeTransformTensors() const;

            //- Return internal point labels
            labelList internalPoints() const;

            //- Return boundary point labels
            labelList boundaryPoints() const;

            //- Return edge length correction
            tmp<edgeScalarField> edgeLengthCorrection() const;

            //- Whether point normals should be corrected for a patch
            bool correctPatchPointNormals(const label patchID) const;

            //- Set whether point normals should be corrected for a patch
            boolList& correctPatchPointNormals() const;


        //- Write mesh
        virtual bool write(const bool valid = true) const;


    // Member Operators

        bool operator!=(const faMesh& m) const;

        bool operator==(const faMesh& m) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faMeshI.H"

#ifdef NoRepository
    #include "faPatchFaMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
