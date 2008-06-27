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
    Foam::triSurfaceMeshes

Description
    Container for triSurfaces read from files.

SourceFiles
    triSurfaceMeshes.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMeshes_H
#define triSurfaceMeshes_H

#include "triSurfaceMesh.H"
#include "fileNameList.H"
#include "treeDataTriSurface.H"
#include "indexedOctree.H"
#include "surfaceLocation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class plane;

/*---------------------------------------------------------------------------*\
                           Class triSurfaceMeshes Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMeshes
:
    public PtrList<triSurfaceMesh>
{
    // Private data

        //- Indices of all surfaces. Precalculated and stored.
        const labelList allSurfaces_;


    // Private Member Functions

        //- Calculate sum of distances to nearest point on surfaces. Is used
        //  in minimisation to find intersection. Returns sum of (square of)
        // distances to the surfaces.
        scalar sumDistSqr
        (
            const labelList& surfacesToTest,
            const scalar initialDistSqr,    // search box
            const point& pt
        ) const;

        //- Takes the tet (points p) and reflects the point with the
        //  highest value around the centre (pSum). Checks if it gets closer
        //  and updates p, y if so.
        scalar tryMorphTet
        (
            const labelList& surfacesToTest,
            const scalar initialDistSqr,
            List<vector>& p,
            List<scalar>& y,
            vector& pSum,
            const label ihi,
            const scalar fac
        ) const;

        //- Downhill simplex method: find the point with min cumulative
        //  distance to all surfaces. Does so by morphing a tet (points p).
        //  Returns the point on the 0th surface or hit if not reached within
        //  maxIters iterations.
        bool morphTet
        (
            const labelList& surfacesToTest,
            const scalar initialDistSqr,
            const scalar convergenceDistSqr,
            const label maxIter,
            List<vector>& p,
            List<scalar>& y
        ) const;

        //- Disallow default bitwise copy construct
        triSurfaceMeshes(const triSurfaceMeshes&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceMeshes&);


public:

    ClassName("triSurfaceMeshes");

    // Constructors

        //- Construct from directory (io.instance()) and list of local filenames
        triSurfaceMeshes(const IOobject& io, const fileNameList&);


    // Member Functions

        //- Get all surfaces in directory
        static fileNameList allNames(const IOobject&);

        //- Get names of surfaces
        fileNameList names() const;


        //- Find any intersection. Return hit point information and surface
        //  number
        label findAnyIntersection
        (
            const labelList& surfacesToTest,
            const point& start,
            const point& end,
            pointIndexHit&
        ) const;

        label findAnyIntersection
        (
            const point& start,
            const point& end,
            pointIndexHit&
        ) const;


        //- Find all intersections in order from start to end. Returns for
        //  every hit the surface and the hit info.
        void findAllIntersections
        (
            const labelList& surfacesToTest,
            const point& start,
            const point& end,
            labelList& surfaces,
            List<pointIndexHit>& surfaceHits
        ) const;

        void findAllIntersections
        (
            const point& start,
            const point& end,
            labelList& surfaces,
            List<pointIndexHit>& surfaceHits
        ) const;


        //- Find intersections of edge nearest to both endpoints.
        void findNearestIntersection
        (
            const labelList& surfacesToTest,
            const point& start,
            const point& end,

            label& surface1,        // surface index
            pointIndexHit& hit1,    // hit point information
            label& surface2,
            pointIndexHit& hit2
        ) const;

        void findNearestIntersection
        (
            const point& start,
            const point& end,

            label& surface1,        // surface index
            pointIndexHit& hit1,    // hit point information
            label& surface2,
            pointIndexHit& hit2
        ) const;


        //- Find nearest. Return -1 (and a miss()) or surface and nearest point.
        label findNearest
        (
            const labelList& surfacesToTest,
            const point&,
            const scalar nearestDistSqr,
            pointIndexHit&
        ) const;

        label findNearest
        (
            const point&,
            const scalar nearestDistSqr,
            pointIndexHit&
        ) const;


        //- Find nearest point (like findNearest) but also return whether
        //  point is on edge or on point. Returns -1 if not found
        //  (within nearestDistSqr) on any of the surfaces.
        //  Returns surfaceLocation
        //  - hit : nearest point is inside triangle
        //  - elementType : none, edge  or point
        //  - index       : triI, edgeI or localPointI
        label findNearestAndClassify
        (
            const labelList& surfacesToTest,
            const point& pt,
            const scalar nearestDistSqr,
            surfaceLocation& nearest
        ) const;

        //- Calculate point which is on a set of surfaces. Takes
        //  - initialDistSqr : search dimensions to find point on surface
        //  - convergenceDistSqr : when point is converged
        // Will bomb out if no stable point reached after certain number
        // of iterations.
        surfaceLocation facesIntersectionTrack
        (
            const labelList& surfacesToTest,
            const scalar initialDistSqr,
            const scalar convergenceDistSqr,
            const point& start
        ) const;

        //- Calculate point which is on a set of surfaces.
        pointIndexHit facesIntersection
        (
            const labelList& surfacesToTest,
            const scalar initialDistSqr,
            const scalar convergenceDistSqr,
            const point& start
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
