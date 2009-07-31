/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2007-2009 OpenCFD Ltd.
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
    CV2D

Description
    Conformal-Voronoi 2D automatic mesher with grid or read initial points
    and point position relaxation with optional "squarification".

    There are a substantial number of options to this mesher read from
    CV2DMesherDict file e.g.:

    // Min cell size used in tolerances when inserting points for
    // boundary conforming.
    // Also used to as the grid spacing usind in insertGrid.
    minCellSize  0.05;

    // Feature angle used to inser feature points
    // 0 = all features, 180 = no features
    featureAngle 45;

    // Maximum quadrant angle allowed at a concave corner before
    // additional "mitering" lines are added
    maxQuadAngle 110;

    // Should the mesh be square-dominated or of unbiased hexagons
    squares      yes;

    // Near-wall region where cells are aligned with the wall specified as a
    // number of cell layers
    nearWallAlignedDist 3;

    // Chose if the cell orientation should relax during the iterations
    //  or remain fixed to the x-y directions
    relaxOrientation    no;

    // Insert near-boundary point mirror or point-pairs
    insertSurfaceNearestPointPairs yes;

    // Mirror near-boundary points rather than insert point-pairs
    mirrorPoints   no;

    // Insert point-pairs vor dual-cell vertices very near the surface
    insertSurfaceNearPointPairs yes;

    // Choose if to randomise the initial grid created by insertGrid.
    randomiseInitialGrid yes;

    // Perturbation fraction, 1 = cell-size.
    randomPurturbation   0.1;

    // Number of relaxation iterations.
    nIterations      5;

    // Relaxation factor at the start of the iteration sequence.
    // 0.5 is a sensible maximum and < 0.2 converges better.
    relaxationFactorStart 0.8;

    // Relaxation factor at the end of the iteration sequence.
    // Should be <= relaxationFactorStart
    relaxationFactorEnd   0;

    writeInitialTriangulation no;
    writeFeatureTriangulation no;
    writeNearestTriangulation no;
    writeInsertedPointPairs   no;
    writeFinalTriangulation   yes;

    // Maximum number of iterations used in boundaryConform.
    maxBoundaryConformingIter 5;

    minEdgeLenCoeff           0.5;
    maxNotchLenCoeff          0.3;
    minNearPointDistCoeff     0.25;
    ppDistCoeff               0.05;

SourceFiles
    CGALTriangulation2Ddefs.H
    indexedVertex.H
    indexedFace.H
    CV2DI.H
    CV2D.C
    CV2DIO.C
    tolerances.C
    controls.C
    insertFeaturePoints.C
    insertSurfaceNearestPointPairs.C
    insertSurfaceNearPointPairs.C
    insertBoundaryConformPointPairs.C

\*---------------------------------------------------------------------------*/

#ifndef CV2D_H
#define CV2D_H

#define CGAL_INEXACT
#define CGAL_HIERARCHY

#include "CGALTriangulation2Ddefs.H"

#include "querySurface.H"
#include "point2DFieldFwd.H"
#include "dictionary.H"
#include "Switch.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CV2D Declaration
\*---------------------------------------------------------------------------*/

class CV2D
:
    public HTriangulation
{
public:

    class controls
    {
    public:

        //- Minimum cell size below which protusions through the surface are
        //  not split
        scalar minCellSize;

        //- Square of minCellSize
        scalar minCellSize2;

        //- The feature angle used to select corners to be
        //  explicitly represented in the mesh.
        //  0 = all features, 180 = no features
        scalar featAngle;

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

        //- Insert point-pairs vor dual-cell vertices very near the surface
        Switch insertSurfaceNearPointPairs;

        Switch writeInitialTriangulation;
        Switch writeFeatureTriangulation;
        Switch writeNearestTriangulation;
        Switch writeInsertedPointPairs;
        Switch writeFinalTriangulation;

        Switch randomiseInitialGrid;
        scalar randomPurturbation;

        label maxBoundaryConformingIter;

        //- Relaxation factor at the start of the iteration
        scalar relaxationFactorStart;

        //- Relaxation factor at the end of the iteration
        scalar relaxationFactorEnd;

        controls(const dictionary& controlDict);
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

        //- Maximum notch size below which protusions into the surface are
        //  not filled
        scalar maxNotchLen;

        //- Square of maxNotchLen
        scalar maxNotchLen2;

        //- The minimum distance alowed between a dual-cell vertex
        // and the surface before a point-pair is introduced
        scalar minNearPointDist;

        //- Square of minNearPoint
        scalar minNearPointDist2;

        //- Distance between boundary conforming point-pairs
        scalar ppDist;

        //- Square of ppDist
        scalar ppDist2;

        tolerances
        (
            const dictionary& controlDict,
            scalar minCellSize,
            const boundBox&
        );
    };


private:

    // Private data

        //- The surface to mesh
        const querySurface& qSurf_;

        //- Meshing controls
        controls controls_;

        //- Meshing tolerances
        tolerances tols_;

        //- z-level
        scalar z_;

        //- Keep track of the start of the internal points
        label startOfInternalPoints_;

        //- Keep track of the start of the surface point-pairs
        label startOfSurfacePointPairs_;

        //- Keep track of the boundary conform point-pairs
        //  stored after the insertion of the surface point-pairs in case
        //  the boundary conform function is called more than once without
        //  removing and insertin the surface point-pairs
        label startOfBoundaryConformPointPairs_;

        //- Temporary storage for a dual-cell
        static const label maxNvert = 20;
        mutable point2D vertices[maxNvert+1];
        mutable vector2D edges[maxNvert+1];


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CV2D(const CV2D&);

        //- Disallow default bitwise assignment
        void operator=(const CV2D&);


        //- Insert point and return it's index
        inline label insertPoint
        (
            const point2D& pt,
            const label type
        );

        inline bool insertMirrorPoint
        (
            const point2D& nearSurfPt,
            const point2D& surfPt
        );

        //- Insert a point-pair at a distance ppDist either side of
        //  surface point point surfPt in the direction n
        inline void insertPointPair
        (
            const scalar mirrorDist,
            const point2D& surfPt,
            const vector2D& n
        );

        //- Create the initial mesh from the bounding-box
        void insertBoundingBox();

        //- Insert point groups at the feature points.
        void insertFeaturePoints();

        //- Insert point-pairs at the given set of points using the surface
        //  normals corresponding to the given set of surface triangles
        //  and write the inserted point locations to the given file.
        void insertPointPairs
        (
            const DynamicList<point2D>& nearSurfacePoints,
            const DynamicList<point2D>& surfacePoints,
            const DynamicList<label>& surfaceTris,
            const fileName fName
        );

        //- Check to see if dual cell specified by given vertex iterator
        //  intersects the boundary and hence reqires a point-pair.
        bool dualCellSurfaceIntersection
        (
            const Triangulation::Finite_vertices_iterator& vit
        ) const;

        //- Insert point-pairs at the nearest points on the surface to the
        //  control vertex of dual-cells which intersect the boundary in order
        //  to provide a boundary-layer mesh.
        //  NB: This is not guaranteed to close the boundary
        void insertSurfaceNearestPointPairs();

        //- Insert point-pairs at small duak-cell edges on the surface in order
        //  to improve the boundary-layer mesh generated by
        //  insertSurfaceNearestPointPairs.
        void insertSurfaceNearPointPairs();

        //- Insert point-pair and correcting the Finite_vertices_iterator
        //  to account for the additional vertices
        void insertPointPair
        (
            Triangulation::Finite_vertices_iterator& vit,
            const point2D& p,
            const label trii
        );

        //- Insert point-pair at the best intersection point between the lines
        //  from the dual-cell real centroid and it's vertices and the surface.
        bool insertPointPairAtIntersection
        (
            Triangulation::Finite_vertices_iterator& vit,
            const point2D& defVert,
            const point2D vertices[],
            const scalar maxProtSize
        );

        //- Insert point-pairs corresponding to dual-cells which intersect
        //  the boundary surface
        label insertBoundaryConformPointPairs(const fileName& fName);

        void markNearBoundaryPoints();

        //- Restore the Delaunay contraint
        void fast_restore_Delaunay(Vertex_handle vh);

        // Flip operations used by fast_restore_Delaunay
        void external_flip(Face_handle& f, int i);
        bool internal_flip(Face_handle& f, int i);


public:

    // Constructors

        //- Construct for given surface
        CV2D(const dictionary& controlDict, const querySurface& qSurf);


    // Destructor

        ~CV2D();


    // Member Functions

        // Access

            const controls& meshingControls() const
            {
                return controls_;
            }


        // Conversion functions between point2D, point and Point

            inline const point2D& toPoint2D(const point&) const;
            inline point toPoint3D(const point2D&) const;

#           ifdef CGAL_INEXACT
                typedef const point2D& point2DFromPoint;
                typedef const Point& PointFromPoint2D;
#           else
                typedef point2D point2DFromPoint;
                typedef Point PointFromPoint2D;
#           endif

            inline point2DFromPoint toPoint2D(const Point&) const;
            inline PointFromPoint2D toPoint(const point2D&) const;
            inline point toPoint3D(const Point&) const;


        // Point insertion

            //- Create the initial mesh from the given internal points.
            //  Points must be inside the boundary by at least nearness
            //  otherwise they are ignored.
            void insertPoints
            (
                const point2DField& points,
                const scalar nearness
            );

            //- Create the initial mesh from the internal points in the given
            //  file.  Points outside the geometry are ignored.
            void insertPoints(const fileName& pointFileName);

            //- Create the initial mesh as a regular grid of points.
            //  Points outside the geometry are ignored.
            void insertGrid();

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


        // Point motion

            inline void movePoint(const Vertex_handle& vh, const Point& P);

            //- Move the internal points to the given new locations and update
            //  the triangulation to ensure it is Delaunay
            // void moveInternalPoints(const point2DField& newPoints);

            //- Calculate the displacements to create the new points
            void newPoints(const scalar relaxation);


        // Write

            //- Write internal points to .obj file
            void writePoints(const fileName& fName, bool internalOnly) const;

            //- Write triangles as .obj file
            void writeTriangles(const fileName& fName, bool internalOnly) const;

            //- Write dual faces as .obj file
            void writeFaces(const fileName& fName, bool internalOnly) const;

            //- Calculates dual points (circumcentres of tets) and faces
            //  (point-cell walk of tets). Returns
            // - dualPoints (in triangle ordering)
            // - dualFaces (compacted)
            void calcDual(point2DField& dualPoints, faceList& dualFaces) const;

            //- Write patch
            void writePatch(const fileName& fName) const;

            void write() const;
};


inline bool boundaryTriangle(const CV2D::Face_handle fc);
inline bool outsideTriangle(const CV2D::Face_handle fc);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CV2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
