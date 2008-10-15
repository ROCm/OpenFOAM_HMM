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
    CV3D

Description

SourceFiles
    CV3DI.H
    CV3D.C
    CV3DIO.C
    controls.C
    tolerances.C
    insertFeaturePoints.C
    insertSurfaceNearestPointPairs.C
    insertSurfaceNearPointPairs.C
    insertBoundaryConformPointPairs.C

\*---------------------------------------------------------------------------*/

#ifndef CV3D_H
#define CV3D_H

//#define CGAL_INEXACT
#define CGAL_HIERARCHY

#include "CGALTriangulation3Ddefs.H"

#include "querySurface.H"
#include "IOdictionary.H"
#include "pointIOField.H"
#include "argList.H"
#include "DynamicList.H"
#include "Switch.H"
#include "Time.H"
#include "polyMesh.H"
#include "SortableList.H"
#include "plane.H"
#include "triSurfaceTools.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class CV3D Declaration
\*---------------------------------------------------------------------------*/

class CV3D
:
    public HTriangulation
{
public:

    class controls
    {
    public:

        //- Relaxation factor at the start of the iteration
        scalar relaxationFactorStart;

        //- Relaxation factor at the end of the iteration
        scalar relaxationFactorEnd;

        //- Minimum cell size below which protusions through the surface are
        //  not split
        scalar minCellSize;

        //- Square of minCellSize
        scalar minCellSize2;

        //- The feature angle used to select corners to be
        //  explicitly represented in the mesh.
        //  0 = no features, 180 = all features
        scalar includedAngle;

        //- Maximum quadrant angle allowed at a concave corner before
        //  additional "mitering" lines are added
        scalar maxQuadAngle;

        //- Should the mesh be square-dominated or of unbiased hexagons
        Switch squares;

        //- Near-wall region where cells are aligned with the wall
        scalar nearWallAlignedDist;

        //- Square of nearWallAlignedDist
        scalar nearWallAlignedDist2;

        //- Chose if the cell orientation should relax during the iterations
        //  or remain fixed to the x-y directions
        Switch relaxOrientation;

        //- Insert near-boundary point mirror or point-pairs
        Switch insertSurfaceNearestPointPairs;

        //- Mirror near-boundary points rather than insert point-pairs
        Switch mirrorPoints;

        //- Insert point-pairs for dual-cell vertices very near the surface
        Switch insertSurfaceNearPointPairs;

        Switch writeInitialTriangulation;
        Switch writeFeatureTriangulation;
        Switch writeNearestTriangulation;
        Switch writeInsertedPointPairs;
        Switch writeFinalTriangulation;

        Switch randomiseInitialGrid;
        scalar randomPerturbation;

        controls(const IOdictionary& controlDict);
    };

    class tolerances
    {
    public:

        //- Maximum cartesian span of the geometry
        scalar span;

        //- Square of span
        scalar span2;

        //- Minumum edge-length of the cell size below which protusions
        //  through the surface are not split
        scalar minEdgeLen;

        //- Square of minEdgeLen
        scalar minEdgeLen2;

        //- Distance between boundary conforming point-pairs
        scalar ppDist;

        //- Square of ppDist
        scalar ppDist2;

        //- Maximum distance between adjacent control points on a feature edge
        //- before they are considered a group
        scalar edgeGroupSpacing;

        //- Guard distance from a feature point within which an edge control
        //- point may not be placed
        scalar featurePointGuard;

        //- Guard distance from a feature edge within which a surface control
        //- point may not be placed
        scalar featureEdgeGuard;

        //- minimum and maximum distances along a feature edge allowed between
        //- pairs of points.  Eventually these should be able to adapt to local
        //- grading requirements.
        scalar minEdgeSpacing;
        scalar maxEdgeSpacing;

        tolerances
        (
            const IOdictionary& controlDict,
            scalar minCellSize,
            const boundBox&
        );
    };

private:

// Private data

        //- The surface to mesh
        const querySurface& qSurf_;

        //- The time registry of the application
        const Time& runTime_;

        //- Meshing controls
        controls controls_;

        //- Meshing tolerances
        tolerances tols_;

        //- Keep track of the start of the internal points
        label startOfInternalPoints_;

        //- Keep track of the start of the surface point-pairs
        label startOfSurfacePointPairs_;

        //- Keep track of the boundary conform point-pairs
        //  stored after the insertion of the surface point-pairs in case
        //  the boundary conform function is called more than once without
        //  removing and insertin the surface point-pairs
        label startOfBoundaryConformPointPairs_;

        //- Store the feature constraining points to be reinserted after a
        //- triangulation clear
        List<Vb> featureConstrainingVertices_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CV3D(const CV3D&);

        //- Disallow default bitwise assignment
        void operator=(const CV3D&);

        //- Insert point and return it's index
        inline label insertPoint
        (
            const point& pt,
            const label type
        );

        inline bool insertMirrorPoint
        (
            const point& nearSurfPt,
            const point& surfPt
        );

        //- Insert a point-pair at a distance ppDist either side of
        //  surface point point surfPt in the direction n
        inline void insertPointPair
        (
            const scalar mirrorDist,
            const point& surfPt,
            const vector& n
        );

        inline void movePoint(const Vertex_handle& vh, const point& p);

        //- Create the initial mesh from the bounding-box
        void insertBoundingBox();

        //- Insert point groups at the feature points.
        void insertFeaturePoints();

        //- Reinsert stored feature point defining points.
        void reinsertFeaturePoints();

        //- Reinsert points that have been saved to clear the mesh.
        void reinsertPoints (const pointField& points);

        //- Insert point-pairs at the given set of points using the surface
        //  normals corresponding to the given set of surface triangles
        //  and write the inserted point locations to the given file.
        void insertPointPairs
        (
            const DynamicList<point>& nearSurfacePoints,
            const DynamicList<point>& surfacePoints,
            const DynamicList<label>& surfaceTris,
            const fileName fName
        );

        //- Insert edge control groups at the given edge points using the
        //  normals of the faces attached to the edge.
        void insertEdgePointGroups
        (
            const DynamicList<point>& edgePoints,
            const DynamicList<label>& edgeLabels,
            const fileName fName
        );

        //- Check to see if dual cell specified by given vertex iterator
        //  intersects the boundary and hence reqires a point-pair.
        bool dualCellSurfaceIntersection
        (
            const Triangulation::Finite_vertices_iterator& vit
        ) const;

        //- Smooths the points positioned along a feature edge by collecting
        //- bunches together and deleting excess points.
        void smoothEdgePositions
        (
            DynamicList<point>& edgePoints,
            DynamicList<label>& edgeLabels
        ) const;

        void smoothEdge
        (
            List<point>& edgePoints,
            List<scalar>& edgeDistances,
            const label edgeI
        ) const;

        //- Insert point-pairs at the nearest points on the surface to the
        //  control vertex of dual-cells which intersect the boundary in order
        //  to provide a boundary-layer mesh.
        //  NB: This is not guaranteed to close the boundary
        void insertSurfaceNearestPointPairs();

        //- Insert point-pairs at small dual-cell edges on the surface in order
        //  to improve the boundary-layer mesh generated by
        //  insertSurfaceNearestPointPairs.
        void insertSurfaceNearPointPairs();

        void markNearBoundaryPoints();

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


public:

    // Constructors

        //- Construct for given surface
        CV3D
        (
            const Time& runTime,
            const querySurface& qSurf
        );

    // Destructor

        ~CV3D();


    // Member Functions

        // Access

            const controls& meshingControls() const
            {
                return controls_;
            }

            // Conversion functions between point and Point

#           ifdef CGAL_INEXACT
                typedef const point& pointFromPoint;
                typedef const Point& PointFrompoint;
#           else
                typedef point pointFromPoint;
                typedef Point PointFrompoint;
#           endif

            inline pointFromPoint topoint(const Point&) const;
            inline PointFrompoint toPoint(const point&) const;

        // Point insertion

            //- Create the initial mesh from the given internal points.
            //  Points must be inside the boundary by at least nearness
            //  otherwise they are ignored.
            void insertPoints
            (
                const pointField& points,
                const scalar nearness
            );

            //- Create the initial mesh from the internal points in the given
            //  file.  Points outside the geometry are ignored.
            void insertPoints(const fileName& pointFileName);

            //- Create the initial mesh as a BCC lattice of points.
            //  Points outside the geometry are ignored.
            void insertGrid();


        // Point motion

            //- Calculate the displacements to relax the points
            void relaxPoints(const scalar relaxation);

            //- Insert all surface point-pairs from
            //  insertSurfaceNearestPointPairs and
            //  findIntersectionForOutsideCentroid
            void insertSurfacePointPairs();

            //- Insert point-pairs where there are protrusions into
            //  or out of the surface
            void boundaryConform();

        // Point removal

            //- Remove the point-pairs introduced by insertSurfacePointPairs
            //  and boundaryConform
            void removeSurfacePointPairs();

        // Check

        // Edit

        // Write

            //- Write Delaunay points to .obj file
            void writePoints(const fileName& fName, bool internalOnly) const;

            //- Write Delaunay triangles as .obj file
            void writeTriangles(const fileName& fName, bool internalOnly) const;

            //- Write dual points and faces as .obj file
            void writeDual(const fileName& fName) const;

            //- Write polyMesh
            void writeMesh(bool writeToTimestep = false);

            void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CV3DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
