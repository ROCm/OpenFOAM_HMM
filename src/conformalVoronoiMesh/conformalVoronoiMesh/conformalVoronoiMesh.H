/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "searchableSurfaces.H"
#include "conformationSurfaces.H"
#include "cvControls.H"
#include "DynamicList.H"
#include "Time.H"
#include "polyMesh.H"
#include "SortableList.H"
#include "meshTools.H"
#include "triSurfaceTools.H"
#include "mathematicalConstants.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class initialPointsMethod;


/*---------------------------------------------------------------------------*\
                    Class conformalVoronoiMesh Declaration
\*---------------------------------------------------------------------------*/

class conformalVoronoiMesh
:
    public HTriangulation
{
    // Private data

        //- The time registry of the application
        const Time& runTime_;

        //- All geometry of the meshing process, including surfaces to be
        //  conformed to and those to be used for refinement
        searchableSurfaces allGeometry_;

        //- The surfaces to conform to
        conformationSurfaces geometryToConformTo_;

        //- Controls for the conformal Voronoi meshing process
        cvControls cvMeshControls_;

        //- Keep track of the start of the internal points
        label startOfInternalPoints_;

        //- Keep track of the start of the surface point-pairs
        label startOfSurfacePointPairs_;

        //- Store the feature constraining points to be reinserted after a
        //- triangulation clear
        List<Vb> featureConstrainingVertices_;

        //- Method for inserting initial points.  Runtime selectable.
        autoPtr<initialPointsMethod> initialPointsMethod_;

    // Private Member Functions

        //- Write the elapsedCpuTime
        void timeCheck() const;

        //- Insert point and return it's index
        inline label insertPoint
        (
            const point& pt,
            const label type
        );

        //- Insert a point-pair at a ppDist distance either side of
        //  surface point surfPt, in the direction n
        inline void insertPointPair
        (
            const scalar ppDist,
            const point& surfPt,
            const vector& n
        );

        //- Insert a Vb (a typedef of CGAL::indexedVertex<K>)
        inline void insertVb(const Vb& v);

        //- Insert pairs of points on the surface with the given normals, at the
        //  specified spacing
        void insertSurfacePointPairs
        (
            const List<scalar>& surfacePpDist,
            const List<point>& surfacePoints,
            const List<vector>& surfaceNormals,
            const fileName fName = fileName::null
        );

        //- Insert point groups at the feature points.
        void insertFeaturePoints();

        //- Reinsert stored feature point defining points.
        void reinsertFeaturePoints();

        //- Insert the initial points into the triangulation, based on the
        //  initialPointsMethod
        void insertInitialPoints();

        //- Check to see if dual cell specified by given vertex iterator
        //  intersects the boundary and hence reqires a point-pair.
        bool dualCellSurfaceIntersection
        (
            const Triangulation::Finite_vertices_iterator& vit
        ) const;

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

        //- Insert the necessary point pairs to conform to the surface
        void conformToSurface();

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

        // Access

            //- Return the Time object
            inline const Time& time() const;

            //- Return the conformationSurfaces object
            inline const conformationSurfaces& geometryToConformTo() const;

            //- Return the cvMeshControls object
            inline const cvControls& cvMeshControls() const;

        // Check

        // Edit

        // Write

            //- Write Delaunay points to .obj file
            void writePoints(const fileName& fName, bool internalOnly) const;

            //- Write list of points to file
            void writePoints
            (
                const fileName& fName,
                const List<point>& points
            ) const;

            //- Write polyMesh
            void writeMesh();

            //- Write dual points and faces as .obj file
            void writeDual
            (
                const pointField& points,
                const faceList& faces,
                const fileName& fName
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "conformalVoronoiMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
