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
    Foam::conformationSurfaces

Description

SourceFiles
    conformationSurfacesI.H
    conformationSurfaces.C
    conformationSurfacesIO.C

\*---------------------------------------------------------------------------*/

#ifndef conformationSurfaces_H
#define conformationSurfaces_H

#include "searchableSurfaces.H"
#include "searchableSurfacesQueries.H"
#include "featureEdgeMesh.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class conformalVoronoiMesh;


/*---------------------------------------------------------------------------*\
                       Class conformationSurfaces Declaration
\*---------------------------------------------------------------------------*/

class conformationSurfaces
{
    // Private data

        //- Reference to the conformalVoronoiMesh holding this cvControls object
        const conformalVoronoiMesh& cvMesh_;

        //- Reference to the searchableSurfaces object holding all geometry data
        const searchableSurfaces& allGeometry_;

        //- A list of featureEdge meshes
        PtrList<featureEdgeMesh> features_;

        //- The location in the mesh that specifies which portion of surfaces is
        //  to be meshed.
        point locationInMesh_;

        //- Indices of surfaces in allGeometry that are to be conformed to
        labelList surfaces_;

        //- A boolean value for each surface to be conformed to specifying if it
        //  is to be treated as a baffle
        boolList baffleSurfaces_;

        //- The overall boundBox of all of the surfaces to be conformed to
        boundBox bounds_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        conformationSurfaces(const conformationSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const conformationSurfaces&);


public:

    // Constructors

        //- Construct from references to conformalVoronoiMesh and
        //- searchableSurfaces
        conformationSurfaces
        (
            const conformalVoronoiMesh& cvMesh,
            const searchableSurfaces& allGeometry,
            const dictionary& surfaceConformationDict
        );

    //- Destructor
    ~conformationSurfaces();


    // Member Functions

        // Access

            //- Return reference to the searchableSurfaces object containing all
            //  of the geometry
            inline const searchableSurfaces& geometry() const;

            //- Return the object holding the feature points and edges
            inline const PtrList<featureEdgeMesh>& features() const;

            //- Return the surface indices
            inline const labelList& surfaces() const;

            //- Return the boundBox
            inline const boundBox& bounds() const;


        // Query

            //- Check if point is inside surfaces to conform to
            Field<bool> inside(const pointField& samplePts) const;

            //- Check if point is outside surfaces to conform to
            Field<bool> outside(const pointField& samplePts) const;

            //- Check if point is inside surfaces to conform to by at least
            //  dist2
            Field<bool> wellInside
            (
                const pointField& samplePts,
                const scalar dist2
            ) const;

            //- Check if point is outside surfaces to conform to by at least
            // dist2
            Field<bool> wellOutside
            (
                const pointField& samplePts,
                const scalar dist2
            ) const;

            // Finding if the line joining start and end intersects the surface
            bool findSurfaceAnyIntersection(point start, point end) const;

            //- Find the nearest point to the sample and return it to the
            //  pointIndexHit, and the normal at the hit location, if found.
            void findSurfaceNearestAndNormal
            (
                const point& sample,
                scalar nearestDistSqr,
                pointIndexHit& pHit,
                vector& normal
            ) const;

            void findSurfaceNearestAndNormal
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>& hitInfo,
                vectorField& normals
            ) const;


            //- Find the nearest point on any feature edge
            void findEdgeNearest
            (
                const point& sample,
                scalar nearestDistSqr,
                pointIndexHit& edgeHit,
                label& featureHit
            ) const;

            void findEdgeNearest
            (
                const pointField& samples,
                const scalarField& nearestDistsSqr,
                List<pointIndexHit>& edgeHits,
                labelList& featuresHit
            ) const;

        // Write

            //- Write all components of all the featureEdgeMeshes as an obj file
            void writeFeatureObj(const fileName& prefix) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "conformationSurfacesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
