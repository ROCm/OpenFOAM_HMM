/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

    As a special exception, you have permission to link this program with the
    CGAL library and distribute executables, as long as you follow the
    requirements of the GNU GPL in regard to all of the software in the
    executable aside from CGAL.

Class
    Foam::conformalVoronoiMesh

Description

SourceFiles
    conformalVoronoiMeshI.H
    conformalVoronoiMesh.C
    conformalVoronoiMeshIO.C
    conformalVoronoiMeshConformToSurface.C
    conformalVoronoiMeshFeaturePoints.C
    conformalVoronoiMeshFeaturePointSpecialisations.C
    conformalVoronoiMeshCalcDualMesh.C

\*---------------------------------------------------------------------------*/

#ifndef conformalVoronoiMesh_H
#define conformalVoronoiMesh_H

#define CGAL_INEXACT

#include "CGALTriangulation3Ddefs.H"
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include "uint.H"
#include "ulong.H"
#include "searchableSurfaces.H"
#include "conformationSurfaces.H"
#include "cellSizeControlSurfaces.H"
#include "cvControls.H"
#include "DynamicList.H"
#include "PackedBoolList.H"
#include "Time.H"
#include "polyMesh.H"
#include "plane.H"
#include "SortableList.H"
#include "meshTools.H"
#include "dynamicIndexedOctree.H"
#include "dynamicTreeDataPoint.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "unitConversion.H"
#include "transform.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Histogram.H"
#include "labelPair.H"
#include "HashSet.H"
#include "memInfo.H"
#include "point.H"
#include "cellSet.H"
#include "wallPolyPatch.H"
#include "processorPolyPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "globalIndex.H"
#include "pointFeatureEdgesTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class initialPointsMethod;

class relaxationModel;

class faceAreaWeightModel;

class backgroundMeshDecomposition;

/*---------------------------------------------------------------------------*\
                    Class conformalVoronoiMesh Declaration
\*---------------------------------------------------------------------------*/

class conformalVoronoiMesh
:
    public Delaunay
{
public:

    enum reconformationMode
    {
        rmNone,   // Do not rebuild the surface conformation
        rmCoarse, // Rebuild the conformation with coarse tolerances (faster)
        rmFine    // Rebuild the conformation with fine tolerances (slower)
    };

    enum faceCollapseMode
    {
        fcmNone,             // Do not collapse face
        fcmEdge,             // Collapse face to a single edge
        fcmPoint,            // Collapse face to a point
        fcmDeferredMultiEdge // Collapse face to several edges
    };


private:

    // Static data

        static const scalar tolParallel;

        static const scalar searchConeAngle;

        static const scalar searchAngleOppositeSurface;


    // Private data

        //- The time registry of the application
        const Time& runTime_;

        //- Random number generator
        mutable Random rndGen_;

        //- Controls for the conformal Voronoi meshing process
        cvControls cvMeshControls_;

        //- All geometry of the meshing process, including surfaces to be
        //  conformed to and those to be used for refinement
        searchableSurfaces allGeometry_;

        //- The surfaces to conform to
        conformationSurfaces geometryToConformTo_;

        //- The cell size control object
        cellSizeControlSurfaces cellSizeControl_;

        //- Keep track of the start index of the internal points.  Needs
        //  careful bookkeeping within several functions.
        label startOfInternalPoints_;

        //- Limiting bound box before infinity begins
        treeBoundBox limitBounds_;

        //- Store the feature constraining points to be reinserted after a
        //  triangulation clear. Maintained with relative types and indices.
        List<Vb> featureVertices_;

        //- Storing the locations of all of the features to be conformed to.
        //  Single pointField required by the featurePointTree.
        pointField featurePointLocations_;

        //- Search tree for feature point locations
        mutable autoPtr<indexedOctree<treeDataPoint> > featurePointTreePtr_;

        //- Search tree for edge point locations
        mutable autoPtr<dynamicIndexedOctree<dynamicTreeDataPoint> >
        edgeLocationTreePtr_;

        //- Search tree for surface point locations
        mutable autoPtr<dynamicIndexedOctree<dynamicTreeDataPoint> >
        surfacePtLocationTreePtr_;

        //- Store locations where the cell size and alignments will be
        //  pre-calculated and looked up
        pointField sizeAndAlignmentLocations_;

        //- Stored cell size at sizeAndAlignmentLocations_
        scalarField storedSizes_;

        //- Stored alignments at sizeAndAlignmentLocations_
        tensorField storedAlignments_;

        //- Search tree for size and alignment lookup points
        mutable autoPtr<indexedOctree<treeDataPoint> > sizeAndAlignmentTreePtr_;

        //- Store the surface and feature edge conformation locations to be
        //  reinserted
        List<Vb> surfaceConformationVertices_;

        //- Method for inserting initial points.  Runtime selectable.
        autoPtr<initialPointsMethod> initialPointsMethod_;

        //- Relaxation coefficient model.  Runtime selectable.
        autoPtr<relaxationModel> relaxationModel_;

        //-  Face area weight function.  Runtime selectable.
        autoPtr<faceAreaWeightModel> faceAreaWeightModel_;

        //- Background mesh decomposition, only available in parallel.
        autoPtr<backgroundMeshDecomposition> decomposition_;


    // Private Member Functions

        //- Return the local target cell size at the given location.  Takes
        //  boolean argument to allow speed-up of queries if the point is going
        //  to be on a surface.
        inline scalar targetCellSize(const Foam::point& pt) const;

        //- Return the target cell size from that stored on a pair of
        //  Delaunay vertices, using a mean function.
        inline scalar averageCellSize
        (
            const Vertex_handle& vA,
            const Vertex_handle& vB
        ) const;

        //- Return the target cell size from that stored on a pair of
        //  Delaunay vertices, including the possibility that one of
        //  them is not an internalOrBoundaryPoint, and so will not
        //  have valid data.
        inline scalar averageAnyCellSize
        (
            const Vertex_handle& vA,
            const Vertex_handle& vB
        ) const;

        //- The average target cell size of a Delaunay facet, i.e., of
        //  a dual edge
        inline scalar averageAnyCellSize
        (
            const Delaunay::Finite_facets_iterator& fit
        ) const;

        //- Return the local point pair separation at the given location
        inline scalar pointPairDistance(const Foam::point& pt) const;

        //- Return the local mixed feature point placement distance
        inline scalar mixedFeaturePointDistance(const Foam::point& pt) const;

        //- Return the square of the local feature point exclusion distance
        inline scalar featurePointExclusionDistanceSqr
        (
            const Foam::point& pt
        ) const;

        //- Return the square of the local feature edge exclusion distance
        inline scalar featureEdgeExclusionDistanceSqr
        (
            const Foam::point& pt
        ) const;

        //- Return the square of the local surface point exclusion distance
        inline scalar surfacePtExclusionDistanceSqr
        (
            const Foam::point& pt
        ) const;

        //- Return the square of the local surface search distance
        inline scalar surfaceSearchDistanceSqr(const Foam::point& pt) const;

        //- Return the local maximum surface protrusion distance
        inline scalar maxSurfaceProtrusion(const Foam::point& pt) const;

        //- Return the required alignment directions at the given location
        tensor requiredAlignment(const Foam::point& pt) const;

        //- Insert Foam::point and return its auto-generated index
        inline label insertPoint
        (
            const Foam::point& p,
            const label type
        );

        //- Insert Point and return its auto-generated index
        inline label insertPoint
        (
            const Point& P,
            const label type
        );

        //- Insert Foam::point with specified index and type
        inline void insertPoint
        (
            const Foam::point& p,
            const label index,
            const label type
        );

        //- Insert Point with specified index and type
        inline void insertPoint
        (
            const Point& P,
            const label index,
            const label type
        );

        //- Insert Delaunay vertices using the CGAL range insertion method,
        //  optionally check processor occupancy and distribute to other
        //  processors
        void insertPoints
        (
            List<Point>& points,
            bool distribute = true
        );

        //- Insert indexed and typed Delaunay vertices, optionally check
        //  processor occupancy and distribute to other processors
        void insertPoints
        (
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types,
            bool distribute = true
        );

        //- Insert a point-pair at a ppDist distance either side of
        //  surface point surfPt, in the direction n
        inline void insertPointPair
        (
            const scalar ppDist,
            const Foam::point& surfPt,
            const vector& n
        );

        //- Create a point-pair at a ppDist distance either side of
        //  surface point surfPt, in the direction n
        inline void createPointPair
        (
            const scalar ppDist,
            const Foam::point& surfPt,
            const vector& n,
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Insert pairs of points on the surface with the given normals, at the
        //  specified spacing
        void insertSurfacePointPairs
        (
            const List<pointIndexHit>& surfaceHits,
            const List<label>& hitSurfaces,
            const fileName fName = fileName::null
        );

        //- Insert groups of points to conform to an edge given a list of
        //  pointIndexHits specifying the location and edge index of the point
        //  to be conformed to on the corresponding entry in featureHit
        void insertEdgePointGroups
        (
            const List<pointIndexHit>& edgeHits,
            const labelList& featuresHit,
            const fileName fName = fileName::null
        );

        //- Call the appropriate function to conform to an edge
        void createEdgePointGroup
        (
            const extendedFeatureEdgeMesh& feMesh,
            const pointIndexHit& edHit,
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Create points to conform to an external edge
        void createExternalEdgePointGroup
        (
            const extendedFeatureEdgeMesh& feMesh,
            const pointIndexHit& edHit,
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Create points to conform to an internal edge
        void createInternalEdgePointGroup
        (
            const extendedFeatureEdgeMesh& feMesh,
            const pointIndexHit& edHit,
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Create points to conform to a flat edge
        void createFlatEdgePointGroup
        (
            const extendedFeatureEdgeMesh& feMesh,
            const pointIndexHit& edHit,
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Create points to conform to an open edge
        void createOpenEdgePointGroup
        (
            const extendedFeatureEdgeMesh& feMesh,
            const pointIndexHit& edHit,
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Create points to conform to multiply connected edge
        void createMultipleEdgePointGroup
        (
            const extendedFeatureEdgeMesh& feMesh,
            const pointIndexHit& edHit,
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Determine and insert point groups at the feature points
        void insertFeaturePoints();

        bool edgesShareNormal(const label e1, const label e2) const;

        //- Create point groups at convex feature points
        void createConvexFeaturePoints
        (
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Create point groups at concave feature points
        void createConcaveFeaturePoints
        (
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Create point groups at mixed feature points
        void createMixedFeaturePoints
        (
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Fill the pointFeatureEdgesType struct with the types of feature
        //  edges that are attached to the point.
        List<extendedFeatureEdgeMesh::edgeStatus> calcPointFeatureEdgesTypes
        (
            const extendedFeatureEdgeMesh& feMesh,
            const labelList& pEds,
            pointFeatureEdgesTypes& pFEdgeTypes
        );

        //- Create feature point groups if a specialisation exists for the
        //  structure
        bool createSpecialisedFeaturePoint
        (
            const extendedFeatureEdgeMesh& feMesh,
            const labelList& pEds,
            const pointFeatureEdgesTypes& pFEdgeTypes,
            const List<extendedFeatureEdgeMesh::edgeStatus>& allEdStat,
            const label ptI,
            DynamicList<Foam::point>& pts,
            DynamicList<label>& indices,
            DynamicList<label>& types
        );

        //- Store the locations of all of the features to be conformed to
        void constructFeaturePointLocations();

        List<pointIndexHit> findSurfacePtLocationsNearFeaturePoint
        (
            const Foam::point& featurePoint
        ) const;

        //- Reinsert stored feature point defining points
        void reinsertFeaturePoints(bool distribute = false);

        //- Demand driven construction of octree for feature points
        const indexedOctree<treeDataPoint>& featurePointTree() const;

        //- Check if a location is in exclusion range around a feature point
        bool nearFeaturePt(const Foam::point& pt) const;

        //- Clear the entire tesselation
        //  Reinsert bounding points, feature points and recalculate
        //  startOfInternalPoints_
        void reset(const bool distribute = false);

        //- Insert far points in a large bounding box to avoid dual edges
        //  spanning huge distances
        void insertBoundingPoints();

        //- Reinsert the bounding points
        void reinsertBoundingPoints();

        //- Insert the initial points into the triangulation, based on the
        //  initialPointsMethod
        void insertInitialPoints();

        //- In parallel redistribute the backgroundMeshDecomposition and
        //  vertices to balance the number of vertices on each processor.
        //  Returns true if the background mesh changes as this removes all
        //  referred vertices, so the parallel interface may need rebuilt.
        bool distributeBackground();

        //- Store data for sizeAndAlignmentLocations_, storedSizes_ and
        //  storedAlignments_ and initialise the sizeAndAlignmentTreePtr_,
        //  determining the appropriate sizeAndAlignmentLocations_
        //  automatically
        void storeSizesAndAlignments();

        //- Store data for sizeAndAlignmentLocations_, storedSizes_ and
        //  storedAlignments_ and initialise the sizeAndAlignmentTreePtr_
        void storeSizesAndAlignments(const List<Point>& storePts);

        //- Restore the sizes and alignments if required
        void updateSizesAndAlignments(const List<Point>& storePts);

        //- Demand driven construction of octree for and alignment points
        const indexedOctree<treeDataPoint>& sizeAndAlignmentTree() const;

        //- Set the size and alignment data for each vertex
        void setVertexSizeAndAlignment();

        //- Builds a dual face by circulating around the supplied edge.
        face buildDualFace
        (
            const Delaunay::Finite_edges_iterator& eit
        ) const;

        //- Finds the maximum filterCount of the dual vertices
        //  (Delaunay cells) that form the dual face produced by the
        //  supplied edge
        label maxFilterCount
        (
            const Delaunay::Finite_edges_iterator& eit
        ) const;

        //- Determines the owner and neighbour labels for dual cells
        //  corresponding to the dual face formed by the supplied
        //  Delaunay vertices.  If the dual face is a boundary face
        //  then neighbour = -1.  Returns true if the dual face
        //  created by vA -> vB needs to be reversed to be correctly
        //  orientated.
        bool ownerAndNeighbour
        (
            Vertex_handle vA,
            Vertex_handle vB,
            label& owner,
            label& neighbour
        ) const;

        //- Insert the necessary point pairs to conform to the surface, either
        //  from stored results, or trigger a re-conformation
        void conformToSurface();

        //- Decision making function for when to rebuild the surface
        //  conformation
        reconformationMode reconformationControl() const;

        //- Create and insert the necessary point pairs to conform to the
        //  surface, then store the result
        void buildSurfaceConformation(reconformationMode reconfMode);

        //- Check to see if dual cell specified by given vertex iterator
        //  intersects the boundary and hence reqires a point-pair
        bool dualCellSurfaceAnyIntersection
        (
            const Delaunay::Finite_vertices_iterator& vit
        ) const;

        //- Return all intersections
        bool dualCellSurfaceAllIntersections
        (
            const Delaunay::Finite_vertices_iterator& vit,
            DynamicList<pointIndexHit>& info,
            DynamicList<label>& hitSurface
        ) const;

        //- Return false if the line is entirely outside the current processor
        //  domain, true is either point is inside, or the processor domain
        //  bounadry is intersected (i.e. the points are box outside but the
        //  line cuts.  The points will be moved onto the box where they
        //  intersect.
        bool clipLineToProc
        (
            const Foam::point& pt,
            Foam::point& a,
            Foam::point& b
        ) const;

        //- Build the parallelInterfaces of the mesh
        void buildParallelInterface
        (
            const word& outputName
        );

        //- Build the parallelInterfaces of the mesh, supply hash sets
        //  externally to allow updates
        void buildParallelInterface
        (
            List<labelHashSet>& referralVertices,
            List<labelHashSet>& receivedVertices,
            bool initialEdgeReferral,
            const word& outputName
        );

        //- Refer all vertices to all processors
        void buildParallelInterfaceAll
        (
            List<labelHashSet>& referralVertices,
            List<labelHashSet>& receivedVertices,
            const word& outputName
        );

        //- Refer vertices where there dual edge pierces other processors
        void buildParallelInterfaceIntersection
        (
            List<labelHashSet>& referralVertices,
            List<labelHashSet>& receivedVertices,
            const word& outputName
        );

        //- Refer vertices that are attached to Delaunay tets whose
        //  circumspheres touch other processors
        void buildParallelInterfaceInfluence
        (
            List<labelHashSet>& referralVertices,
            List<labelHashSet>& receivedVertices,
            const word& outputName
        );

        //- Refer vertices to their required processors
        void referVertices
        (
            const DynamicList<label>& targetProcessor,
            DynamicList<Foam::point>& parallelPoints,
            DynamicList<label>& parallelIndices,
            List<labelHashSet>& receivedVertices,
            const word& stageName,
            const word& outputName
        );

        //- Find the "worst" protrusion of a dual cell through the surface,
        //  subject to the maxSurfaceProtrusion tolerance
        void dualCellLargestSurfaceProtrusion
        (
            const Delaunay::Finite_vertices_iterator& vit,
            pointIndexHit& surfHit,
            label& hitSurface
        ) const;

        //- Find the "worst" incursion of the dual cell of a non-internal or
        //  boundary point through the surface, subject to the
        //  maxSurfaceProtrusion tolerance
        void dualCellLargestSurfaceIncursion
        (
            const Delaunay::Finite_vertices_iterator& vit,
            pointIndexHit& surfHitLargest,
            label& hitSurfaceLargest
        ) const;

        //- Write out vertex-processor occupancy information for debugging
        void reportProcessorOccupancy();

        //- Write out debugging information about the surface conformation
        //  quality
        void reportSurfaceConformationQuality();

        //- Limit the displacement of a point so that it doesn't penetrate the
        //  surface to be meshed or come too close to it
        void limitDisplacement
        (
            const Delaunay::Finite_vertices_iterator& vit,
            vector& displacement,
            label callCount = 0
        ) const;

        //- Find angle between the normals of two close surface points.
        scalar angleBetweenSurfacePoints(Foam::point pA, Foam::point pB) const;

        //- Check if a surface point is near another.
        bool nearSurfacePoint
        (
            pointIndexHit& pHit,
            label& surfaceHit,
            DynamicList<Foam::point>& existingSurfacePtLocations
        ) const;

        //- Append a point to the surface point tree and the existing list
        bool appendToSurfacePtTree
        (
            const Foam::point& pt,
            DynamicList<Foam::point>& existingSurfacePtLocations
        ) const;

        //- Append a point to the edge location tree and the existing list
        bool appendToEdgeLocationTree
        (
            const Foam::point& pt,
            DynamicList<Foam::point>& existingEdgeLocations
        ) const;

        //- Check if a point is near any feature edge points.
        bool pointIsNearFeatureEdgeLocation(const Foam::point& pt) const;

        bool pointIsNearFeatureEdgeLocation
        (
            const Foam::point& pt,
            pointIndexHit& info
        ) const;

        //- Check if a point is near any surface conformation points.
        bool pointIsNearSurfaceLocation(const Foam::point& pt) const;

        bool pointIsNearSurfaceLocation
        (
            const Foam::point& pt,
            pointIndexHit& info
        ) const;

        //- Check if a location is in the exclusion range of an existing feature
        //- edge conformation location
        bool nearFeatureEdgeLocation
        (
            pointIndexHit& pHit,
            DynamicList<Foam::point>& existingEdgeLocations
        ) const;

        //- Build or rebuild the edge location tree
        void buildEdgeLocationTree
        (
            const DynamicList<Foam::point>& existingEdgeLocations
        ) const;

        //- Build or rebuild the surface point location tree
        void buildSurfacePtLocationTree
        (
            const DynamicList<Foam::point>& existingSurfacePtLocations
        ) const;

        //- Build or rebuild the sizeAndAlignmentTree
        void buildSizeAndAlignmentTree() const;

        //- Process the surface conformation locations to decide which surface
        //  and edge conformation locations to add
        void addSurfaceAndEdgeHits
        (
            const Delaunay::Finite_vertices_iterator& vit,
            const Foam::point& vert,
            const DynamicList<pointIndexHit>& surfHit,
            const DynamicList<label>& hitSurface,
            scalar surfacePtReplaceDistCoeffSqr,
            scalar edgeSearchDistCoeffSqr,
            DynamicList<pointIndexHit>& surfaceHits,
            DynamicList<label>& hitSurfaces,
            DynamicList<pointIndexHit>& featureEdgeHits,
            DynamicList<label>& featureEdgeFeaturesHit,
            DynamicList<Foam::point>& newEdgeLocations,
            DynamicList<Foam::point>& existingEdgeLocations,
            DynamicList<Foam::point>& existingSurfacePtLocations
        ) const;

        //- Store the surface conformation with the indices offset to be
        //  relative to zero
        void storeSurfaceConformation();

        //- Reinsert the surface conformation re-offsetting indices to be
        //  relative to new number of internal vertices
        void reinsertSurfaceConformation();

        //- Dual calculation
        void calcDualMesh
        (
            pointField& points,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            wordList& patchTypes,
            wordList& patchNames,
            labelList& patchSizes,
            labelList& patchStarts,
            labelList& procNeighbours,
            pointField& cellCentres,
            labelList& cellToDelaunayVertex,
            labelListList& patchToDelaunayVertex,
            bool filterFaces
        );

        //- Tet mesh calculation
        void calcTetMesh
        (
            pointField& points,
            labelList& pointToDelaunayVertex,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            wordList& patchTypes,
            wordList& patchNames,
            labelList& patchSizes,
            labelList& patchStarts
        );

        //- Determines if the dual face constructed by the Delaunay
        //  edge is a boundary face
        inline bool isBoundaryDualFace
        (
            const Delaunay::Finite_edges_iterator& eit
        ) const;

        //- Which processors are attached to the dual edge represented by this
        //  Delaunay facet
        inline List<label> processorsAttached
        (
            const Delaunay::Finite_facets_iterator& fit
        ) const;

        //- Determines if the edge constructed from the face is on
        //  a processor patch
        inline bool isParallelDualEdge
        (
            const Delaunay::Finite_facets_iterator& fit
        ) const;

        //- Merge vertices that are very close together
        void mergeCloseDualVertices
        (
            const pointField& pts,
            const PackedBoolList& boundaryPts
        );

        label mergeCloseDualVertices
        (
            const pointField& pts,
            const PackedBoolList& boundaryPts,
            Map<label>& dualPtIndexMap
        ) const;

        //- Smooth the dual vertices of the dual faces on the boundary
        //  so that they conform to the surface and remove any
        //  small, normal oriented faces
        void smoothSurface
        (
            pointField& pts,
            const PackedBoolList& boundaryPts
        );

        label smoothSurfaceDualFaces
        (
            pointField& pts,
            const PackedBoolList& boundaryPts,
            Map<label>& dualPtIndexMap
        ) const;

        //- Collapse a face to a point or an edge, modifying and
        //  mapping the points, returns the true if the face was
        //  collapsed in this operation
        void collapseFaces
        (
            pointField& pts,
            const PackedBoolList& boundaryPts,
            HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
        );

        label collapseFaces
        (
            pointField& pts,
            const PackedBoolList& boundaryPts,
            Map<label>& dualPtIndexMap,
            HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
        ) const;

        //- Collapse a face to an edge, updating the point and point
        //  map. Returns the collapse mode that was applied.
        faceCollapseMode collapseFace
        (
            const face& f,
            pointField& pts,
            const PackedBoolList& boundaryPts,
            Map<label>& dualPtIndexMap,
            scalar targetFaceSize,
            scalar collapseSizeLimitCoeff,
            label maxFC
        ) const;

        //- Identify the index of the longest edge on the face
        label longestEdge(const face& f, const pointField& pts) const;

        //- Identify the face labels of the deferred collapse faces
        void deferredCollapseFaceSet
        (
            labelList& owner,
            labelList& neighbour,
            const HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
        ) const;

        //- Create a polyMesh and check its quality, reports which
        //  elements damage the mesh quality, allowing backtracking.
        labelHashSet checkPolyMeshQuality(const pointField& pts) const;

        //- Index all of the the Delaunay cells and calculate their
        //- dual points
        void indexDualVertices
        (
            pointField& pts,
            PackedBoolList& boundaryPts
        );

        //- Re-index all of the the Delaunay cells
        void reindexDualVertices(const Map<label>& dualPtIndexMap);

        //- Create all of the internal and boundary faces
        void createFacesOwnerNeighbourAndPatches
        (
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            wordList& patchTypes,
            wordList& patchNames,
            labelList& patchSizes,
            labelList& patchStarts,
            labelList& procNeighbours,
            labelListList& patchPointPairSlaves,
            bool includeEmptyPatches = false
        ) const;

        //- Create the cell centres to use for the mesh
        void createCellCentres(pointField& cellCentres) const;

        //- Extract all points in vertex-index order
        tmp<pointField> allPoints() const;

        //- Sort the faces, owner and neighbour lists into
        //  upper-triangular order.  For internal faces only, use
        //  before adding patch faces
        void sortFaces
        (
            faceList& faces,
            labelList& owner,
            labelList& neighbour
        ) const;

        //- Sort the processor patches so that the faces are in the same order
        //  on both processors
        void sortProcPatches
        (
            List<DynamicList<face> >& patchFaces,
            List<DynamicList<label> >& patchOwners,
            List<DynamicList<label> >& patchPointPairSlaves,
            List<Pair<DynamicList<label> > >& patchSortingIndices
        ) const;

        //- Add the faces and owner information for the patches
        void addPatches
        (
            const label nInternalFaces,
            faceList& faces,
            labelList& owner,
            labelList& patchSizes,
            labelList& patchStarts,
            List<DynamicList<face> >& patchFaces,
            List<DynamicList<label> >& patchOwners
        ) const;

        //- Remove points that are no longer used by any faces
        void removeUnusedPoints
        (
            faceList& faces,
            pointField& pts
        ) const;

        //- Remove dual cells that are not used by any faces. Return compaction
        //  map.
        labelList removeUnusedCells
        (
            labelList& owner,
            labelList& neighbour
        ) const;

        //- Disallow default bitwise copy construct
        conformalVoronoiMesh(const conformalVoronoiMesh&);

        //- Disallow default bitwise assignment
        void operator=(const conformalVoronoiMesh&);


public:

    //- Runtime type information
    ClassName("conformalVoronoiMesh");


    // Constructors

        //- Construct from Time and cvMeshDict
        conformalVoronoiMesh
        (
            const Time& runTime,
            const dictionary& cvMeshDict
        );


    //- Destructor
    ~conformalVoronoiMesh();


    // Member Functions

        //- Move the vertices according to the controller, re-conforming to the
        //  surface as required
        void move();

        //- Check if the point is in the domain handled by this processor
        bool positionOnThisProc(const Foam::point& pt) const;

        //- Check if the point is in the domain handled by this processor
        boolList positionOnThisProc(const Foam::List<Foam::point>& pts) const;

        //- Which processor's domain handles this point
        labelList positionProc(const Foam::List<Foam::point>& pts) const;

        //- Which other processors does each line segment intersect
        List<List<pointIndexHit> > intersectsProc
        (
            const List<Foam::point>& starts,
            const List<Foam::point>& ends
        ) const;

        //- Which other processors does each sphere overlap
        labelListList overlapsProc
        (
            const List<Foam::point>& centres,
            const List<scalar>& radiusSqrs
        ) const;

        //- Conversion functions between point (FOAM::) and Point (CGAL)

#       ifdef CGAL_INEXACT
            typedef const Foam::point& pointFromPoint;
            typedef const Point& PointFrompoint;
#       else
            typedef Foam::point pointFromPoint;
            typedef Point PointFrompoint;
#       endif

        inline pointFromPoint topoint(const Point&) const;
        inline PointFrompoint toPoint(const Foam::point&) const;

        typedef K::Vector_3 CGALVector;

        inline CGALVector toCGALVector(const Foam::vector& v) const;

        // Access

            //- Return the Time object
            inline const Time& time() const;

            //- Return the random number generator
            inline Random& rndGen() const;

            //- Return the allGeometry object
            inline const searchableSurfaces& allGeometry() const;

            //- Return the conformationSurfaces object
            inline const conformationSurfaces& geometryToConformTo() const;

            //- Return the backgroundMeshDecomposition
            inline const backgroundMeshDecomposition& decomposition() const;

            //- Return the cellSizeControlSurfaces object
            inline const cellSizeControlSurfaces& cellSizeControl() const;

            //- Return the cvMeshControls object
            inline const cvControls& cvMeshControls() const;


        // Write

            //- Write the elapsedCpuTime and memory usage, with an optional
            //  description
            void timeCheck
            (
                const string& description = string::null
            ) const;

            //- Write the Delaunay cell
            void drawDelaunayCell
            (
                Ostream& os,
                const Cell_handle& c,
                label offset = 0
            ) const;

            //- Write Delaunay points to .obj file
            void writePoints(const fileName& fName, bool internalOnly) const;

            //- Write the boundary Delaunay points to .obj file
            void writeBoundaryPoints(const fileName& fName) const;

            //- Write list of points to file
            void writePoints
            (
                const fileName& fName,
                const List<Foam::point>& points
            ) const;

            //- Write the internal Delaunay vertices of the tessellation as a
            //  pointField that may be used to restart the meshing process
            void writeInternalDelaunayVertices(const fileName& instance) const;

            //- Prepare data and call writeMesh for polyMesh and
            //  tetDualMesh
            void writeMesh
            (
                const fileName& instance,
                bool filterFaces = true
            );

            //- Write mesh to disk
            void writeMesh
            (
                const word& meshName,
                const fileName& instance,
                pointField& points,
                faceList& faces,
                labelList& owner,
                labelList& neighbour,
                const wordList& patchTypes,
                const wordList& patchNames,
                const labelList& patchSizes,
                const labelList& patchStarts,
                const labelList& procNeighbours,
                const pointField& cellCentres
            ) const;

            //- Write points and faces as .obj file
            void writeObjMesh
            (
                const pointField& points,
                const faceList& faces,
                const fileName& fName
            ) const;

            //- Calculate and write a field of the target cell size,
            //  target cell volume, actual cell volume and equivalent
            //  actual cell size (cbrt(actual cell volume)).
            void writeCellSizes(const fvMesh& mesh) const;

            //- Calculate and write the cell centres.
            void writeCellCentres(const fvMesh& mesh) const;

            //- Find the cellSet of the boundary cells which have points that
            //  protrude out of the surface beyond a tolerance.
            void findRemainingProtrusionSet(const fvMesh& mesh) const;


            //- Function inserting points into a triangulation and setting the
            //  index and type data of the point in the correct order. This is
            //  faster than inserting points individually.
            //
            //  Adapted from a post on the CGAL lists: 2010-01/msg00004.html by
            //  Sebastien Loriot (Geometry Factory).
            //
            //  @todo problems putting declaration in the .C file. Function
            //        prototype is not recognised.
            template<class Triangulation, class Point_iterator>
            void rangeInsertWithInfo
            (
                Point_iterator begin,
                Point_iterator end,
                Triangulation& T,
                DynamicList<label>& indices,
                DynamicList<label>& types
            )
            {
                typedef std::vector
                <
                    std::pair<const typename Triangulation::Point*, label>
                > vectorPairPointIndex;

                vectorPairPointIndex points;
                label index = 0;

                for (Point_iterator it = begin; it != end; ++it)
                {
                    points.push_back
                    (
                        std::make_pair(&(toPoint(*it)), index++)
                    );
                }

                std::random_shuffle(points.begin(), points.end());

                spatial_sort
                (
                    points.begin(),
                    points.end(),
                    Traits_for_spatial_sort<Triangulation>()
                );

                typename Triangulation::Cell_handle hint;

                for
                (
                    typename vectorPairPointIndex::const_iterator
                    p = points.begin();
                    p != points.end();
                    ++p
                )
                {
                    typename Triangulation::Locate_type lt;
                    typename Triangulation::Cell_handle c;
                    label li, lj;

                    c = T.locate(*(p->first), lt, li, lj, hint);

                    typename Triangulation::Vertex_handle v
                        = T.insert(*(p->first), lt, c, li, lj);

                    label oldIndex = p->second;

                    label type = types[oldIndex];

                    if (type > Vb::vtFar)
                    {
                        // This is a member of a point pair, don't use the type
                        // directly (note that this routine never gets called
                        // for referredPoints so type will never be -procI)
                        type += T.number_of_vertices() - 1;
                    }

                    v->index() = indices[oldIndex] + T.number_of_vertices() - 1;
                    v->type() = type;
                }
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "conformalVoronoiMeshI.H"

//#ifdef NoRepository
//#   include "conformalVoronoiMeshTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
