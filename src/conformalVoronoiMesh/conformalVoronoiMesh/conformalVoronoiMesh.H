/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::conformalVoronoiMesh

Description

SourceFiles
    conformalVoronoiMeshI.H
    conformalVoronoiMesh.C
    conformalVoronoiMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef conformalVoronoiMesh_H
#define conformalVoronoiMesh_H

//#define CGAL_INEXACT
#define CGAL_HIERARCHY

#include "CGALTriangulation3Ddefs.H"
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
#include "triSurfaceTools.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "mathematicalConstants.H"
#include "transform.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Histogram.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class initialPointsMethod;

class relaxationModel;

class faceAreaWeightModel;


/*---------------------------------------------------------------------------*\
                    Class conformalVoronoiMesh Declaration
\*---------------------------------------------------------------------------*/

class conformalVoronoiMesh
:
    public HTriangulation
{
public:

    enum reconformationMode
    {
        NONE,       // Do not rebuild the surface conformation
        COARSE,     // Rebuild the conformation with coarse tolerances (faster)
        FINE        // Rebuild the conformation with fine tolerances (slower)
    };

private:

    // Private data

        //- The time registry of the application
        const Time& runTime_;

        //- All geometry of the meshing process, including surfaces to be
        //  conformed to and those to be used for refinement
        searchableSurfaces allGeometry_;

        //- The surfaces to conform to
        conformationSurfaces geometryToConformTo_;

        //- The cell size control object
        cellSizeControlSurfaces cellSizeControl_;

        //- Controls for the conformal Voronoi meshing process
        cvControls cvMeshControls_;

        //- Keep track of the start of the internal points
        label startOfInternalPoints_;

        //- Keep track of the start of the surface point-pairs
        label startOfSurfacePointPairs_;

        //- Store the feature constraining points to be reinserted after a
        //  triangulation clear
        List<Vb> featureVertices_;

        //- Storing the locations of all of the features to be conformed to.
        //  Single pointField required by the featurePointTree_.
        pointField featurePointLocations_;

        //- Search tree for feature point locations
        mutable autoPtr<indexedOctree<treeDataPoint> > featurePointTree_;

        //- Store locations where the cell size and alignments will be
        //  pre-calculated and looked up
        pointField sizeAndAlignmentLocations_;

        //- Stored cell size at sizeAndAlignmentLocations_
        scalarField storedSizes_;

        //- Stored alignments at sizeAndAlignmentLocations_
        tensorField storedAlignments_;

        //- Search tree for size and alignment lookup points
        mutable autoPtr<indexedOctree<treeDataPoint> > sizeAndAlignmentTree_;

        //- Store the surface and feature edge conformation locations to be
        //  reinserted
        List<Vb> surfaceConformationVertices_;

        //- Method for inserting initial points.  Runtime selectable.
        autoPtr<initialPointsMethod> initialPointsMethod_;

        //- Relaxation coefficient model.  Runtime selectable
        autoPtr<relaxationModel> relaxationModel_;

        //-  Face area weight function.  Runtime selectable.
        autoPtr<faceAreaWeightModel> faceAreaWeightModel_;


    // Private Member Functions

        //- Write the elapsedCpuTime
        inline void timeCheck() const;

        //- Return the local target cell size at the given location.  Takes
        //  boolean argument to allow speed-up of queries if the point is going
        //  to be on a surface.
        inline scalar targetCellSize
        (
            const point& pt,
            bool isSurfacePoint = false
        ) const;

        //- Return the local point pair separation at the given location
        inline scalar pointPairDistance(const point& pt) const;

        //- Return the local mixed feature point placement distance
        inline scalar mixedFeaturePointDistance(const point& pt) const;

        //- Return the square of the local feature point exclusion distance
        inline scalar featurePointExclusionDistanceSqr(const point& pt) const;

        //- Return the square of the local feature edge exclusion distance
        inline scalar featureEdgeExclusionDistanceSqr(const point& pt) const;

        //- Return the square of the local surface search distance
        inline scalar surfaceSearchDistanceSqr(const point& pt) const;

        //- Return the local maximum surface protrusion distance
        inline scalar maxSurfaceProtrusion(const point& pt) const;

        //- Return the local minimum allowed dual edge length
        inline scalar minimumEdgeLength(const point& pt) const;

        //- Return the required alignment directions at the given location
        tensor requiredAlignment(const point& pt) const;

        //- Insert point and return it's index
        inline label insertPoint
        (
            const point& pt,
            const label type
        );

        //- Insert a std::vector of CGAL Points using the CGAL range insertion
        //  method
        inline void insertPoints(const std::vector<Point>& points);

        //- Insert a point-pair at a ppDist distance either side of
        //  surface point surfPt, in the direction n
        inline void insertPointPair
        (
            const scalar ppDist,
            const point& surfPt,
            const vector& n
        );

        //- Insert a Vb (a typedef of CGAL::indexedVertex<K>) with the
        //  possibility of an offset for the index and the type
        inline void insertVb(const Vb& v, label offset = 0);

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
        void insertEdgePointGroup
        (
            const featureEdgeMesh& feMesh,
            const pointIndexHit& edHit
        );

        //- Insert points to conform to an external edge
        void insertExternalEdgePointGroup
        (
            const featureEdgeMesh& feMesh,
            const pointIndexHit& edHit
        );

        //- Insert points to conform to an internal edge
        void insertInternalEdgePointGroup
        (
            const featureEdgeMesh& feMesh,
            const pointIndexHit& edHit
        );

        //- Insert points to conform to a flat edge
        void insertFlatEdgePointGroup
        (
            const featureEdgeMesh& feMesh,
            const pointIndexHit& edHit
        );

        //- Insert points to conform to an open edge
        void insertOpenEdgePointGroup
        (
            const featureEdgeMesh& feMesh,
            const pointIndexHit& edHit
        );

        //- Insert points to conform to multiply connected edge
        void insertMultipleEdgePointGroup
        (
            const featureEdgeMesh& feMesh,
            const pointIndexHit& edHit
        );

        //- Determine and insert point groups at the feature points.
        void createFeaturePoints();

        //- Insert point groups at convex feature points
        void insertConvexFeaturesPoints();

        //- Insert point groups at concave feature points
        void insertConcaveFeaturePoints();

        //- Insert point groups at mixed feature points
        void insertMixedFeaturePoints();

        //- Store the locations of all of the features to be conformed to.
        void constructFeaturePointLocations();

        //- Reinsert stored feature point defining points.
        void reinsertFeaturePoints();

        //- Demand driven construction of octree for feature points
        const indexedOctree<treeDataPoint>& featurePointTree() const;

        //- Check if a location is in exclusion range around a feature point
        bool nearFeaturePt(const point& pt) const;

        //- Insert the initial points into the triangulation, based on the
        //  initialPointsMethod
        void insertInitialPoints();

        //- Store data for sizeAndAlignmentLocations_, storedSizes_ and
        //  storedAlignments_ and initialise the sizeAndAlignmentTree_
        void storeSizesAndAlignments(const std::vector<Point>& initPts);

        //- Demand driven construction of octree for and alignment points
        const indexedOctree<treeDataPoint>& sizeAndAlignmentTree() const;

        //- Set the size and alignment data for each vertex
        void setVertexSizeAndAlignment();

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
        //  intersects the boundary and hence reqires a point-pair.
        bool dualCellSurfaceAnyIntersection
        (
            const Triangulation::Finite_vertices_iterator& vit
        ) const;

        //- Find the "worst" protrusion of a dual cell through the surface,
        //  subject to the tolerance
        void dualCellLargestSurfaceProtrusion
        (
            const Triangulation::Finite_vertices_iterator& vit,
            pointIndexHit& surfHit,
            label& hitSurface
        ) const;

        //- Limit the displacement of a point so that it doesn't penetrate the
        //  surface to be meshed or come too close to it.
        void limitDisplacement
        (
            const Triangulation::Finite_vertices_iterator& vit,
            vector& displacement
        ) const;

        //- Check if a location is in the exclusion range of an existing feature
        //- edge conformation location
        bool nearFeatureEdgeLocation
        (
            const point& pt,
            DynamicList<point>& newEdgeLocations,
            pointField& existingEdgeLocations,
            autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree
        ) const;

        //- Build or rebuild the edgeLocationTree
        void buildEdgeLocationTree
        (
            autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree,
            const pointField& existingEdgeLocations
        ) const;

        //- Build or rebuild the sizeAndAlignmentTree
        void buildSizeAndAlignmentTree() const;

        //- Process the surface conformation locations to decide which surface
        //  and edge conformation locations to add.
        void addSurfaceAndEdgeHits
        (
            const Triangulation::Finite_vertices_iterator& vit,
            const point& vert,
            const pointIndexHit& surfHit,
            label hitSurface,
            scalar surfacePtReplaceDistCoeffSqr,
            scalar edgeSearchDistCoeffSqr,
            DynamicList<pointIndexHit>& surfaceHits,
            DynamicList<label>& hitSurfaces,
            DynamicList<pointIndexHit>& featureEdgeHits,
            DynamicList<label>& featureEdgeFeaturesHit,
            DynamicList<point>& newEdgeLocations,
            pointField& existingEdgeLocations,
            autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree
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
            wordList& patchNames,
            labelList& patchSizes,
            labelList& patchStarts
        );

        //- Disallow default bitwise copy construct
        conformalVoronoiMesh(const conformalVoronoiMesh&);

        //- Disallow default bitwise assignment
        void operator=(const conformalVoronoiMesh&);


public:

    // Constructors

        //- Construct null
        conformalVoronoiMesh();

        //- Construct from Time and cvMeshDict
        conformalVoronoiMesh
        (
            const Time& runTime,
            const IOdictionary& cvMeshDict
        );


    //- Destructor
    ~conformalVoronoiMesh();


    // Member Functions

        //- Move the vertices according to the controller, re-conforming to the
        //  surface as required
        void move();

        //- Conversion functions between point (OpenFOAM) and Point (CGAL)

#       ifdef CGAL_INEXACT
            typedef const point& pointFromPoint;
            typedef const Point& PointFrompoint;
#       else
            typedef point pointFromPoint;
            typedef Point PointFrompoint;
#       endif

        inline pointFromPoint topoint(const Point&) const;
        inline PointFrompoint toPoint(const point&) const;

        typedef K::Vector_3 CGALVector;

        inline CGALVector toCGALVector(const point& pt) const;

        // Access

            //- Return the Time object
            inline const Time& time() const;

            //- Return the conformationSurfaces object
            inline const conformationSurfaces& geometryToConformTo() const;

            //- Return the cellSizeControlSurfaces object
            inline const cellSizeControlSurfaces& cellSizeControl() const;

            //- Return the cvMeshControls object
            inline const cvControls& cvMeshControls() const;


        // Write

            //- Write Delaunay points to .obj file
            void writePoints(const fileName& fName, bool internalOnly) const;

            //- Write list of points to file
            void writePoints
            (
                const fileName& fName,
                const List<point>& points
            ) const;

            //- Write the internal Delaunay vertices of the tessellation as a
            //  pointField that may be used to restart the meshing process
            void writeInternalDelaunayVertices(bool writeToConstant) const;

            //- Write polyMesh
            void writeMesh(bool writeToConstant = true);

            //- Write dual points and faces as .obj file
            void writeDual
            (
                const pointField& points,
                const faceList& faces,
                const fileName& fName
            ) const;

            //- Read mesh from file as an fvMesh then calculate and write a
            //  field of the target cell size
            void writeTargetCellSize() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "conformalVoronoiMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
