/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

Class
    Foam::conformationSurfaces

Description

SourceFiles
    conformationSurfacesI.H
    conformationSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef conformationSurfaces_H
#define conformationSurfaces_H

#include "searchableSurfaces.H"
#include "searchableSurfacesQueries.H"
#include "extendedFeatureEdgeMesh.H"
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

        //- Reference to the conformalVoronoiMesh holding this object
        const conformalVoronoiMesh& cvMesh_;

        //- Reference to the searchableSurfaces object holding all geometry data
        const searchableSurfaces& allGeometry_;

        //- A list of extendedFeatureEdgeMesh
        PtrList<extendedFeatureEdgeMesh> features_;

        //- The location in the mesh that specifies which portion of surfaces is
        //  to be meshed.
        point locationInMesh_;

        //- Indices of surfaces in allGeometry that are to be conformed to
        labelList surfaces_;

        //- Reverse mapping, which entry in surfaces corresponds to the surface
        //  in allGeometry specified by the list index.  -1 for a surface that
        //  isn't used.
        labelList allGeometryToSurfaces_;

        //- A boolean value for each surface to be conformed to specifying if it
        //  is to be treated as a baffle
        boolList baffleSurfaces_;

        //- A flat list of all of the names of the patches from all of the
        //  surfaces to be reproduced in the meshed geometry
        List<word> patchNames_;

        //- The offset between the patch indices of the individual surface and
        //  the entry in the overall patch list
        labelList patchOffsets_;

        //- The boundBox of volume to be conformed to on this processor
        treeBoundBox bounds_;

        //- The overall boundBox of all of the surfaces to be conformed to
        treeBoundBox globalBounds_;

        //- Magnitude of the span of the domain
        scalar spanMag_;

        //- Square of span_
        scalar spanMagSqr_;

        //- All of the treeBoundBoxes representing volumes held by all
        //  processors.  A non-overlapping set of volumes.
        List<treeBoundBoxList> processorDomains_;

        //- The bounds of the mesh on each processor - can overlap as the mesh
        //  is generated, used to determine which vertices to refer.
        treeBoundBoxList processorMeshBounds_;

        //- The pattern/signature of volumeTypes representing a point in the
        // domain to be meshed
        List<searchableSurface::volumeType> referenceVolumeTypes_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        conformationSurfaces(const conformationSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const conformationSurfaces&);


public:

    // Constructors

        //- Construct from dictionary and references to conformalVoronoiMesh and
        //  searchableSurfaces
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
            inline const PtrList<extendedFeatureEdgeMesh>& features() const;

            //- Return the surface indices
            inline const labelList& surfaces() const;

            //- Return the patch names
            inline const List<word>& patchNames() const;

            //- Return the bounds
            inline const treeBoundBox& bounds() const;

            //- Return the global bounds
            inline const treeBoundBox& globalBounds() const;

            //- Return the spanMag
            inline scalar spanMag() const;

            //- Return spanSqr
            inline scalar spanMagSqr() const;

            //- Return the processorDomains
            inline const List<treeBoundBoxList>& processorDomains() const;

            //- Return the processorMeshBounds
            inline const treeBoundBoxList& processorMeshBounds() const;

            //- Return non-const access to the processorMeshBounds
            inline treeBoundBoxList& processorMeshBounds();


        // Query

            //- Check if the supplied bound box overlaps any part of any of
            //  the surfaces
            bool overlaps(const treeBoundBox& bb) const;

            //- Check if the point is in the domain handled by this processor
            bool positionOnThisProc(const point& pt) const;

            //- Which processor's domain handles this point
            label positionProc(const point& pt) const;

            //- Check if points are inside surfaces to conform to
            Field<bool> inside(const pointField& samplePts) const;

            //- Check if point is inside surfaces to conform to
            bool inside(const point& samplePt) const;

            //- Check if points are outside surfaces to conform to
            Field<bool> outside(const pointField& samplePts) const;

            //- Check if point is outside surfaces to conform to
            bool outside(const point& samplePt) const;

            //- Check if point is closer to the surfaces to conform to than
            //  testDistSqr, in which case return false, otherwise assess in or
            //  outside and erturn a result depending on the testForInside flag
            Field<bool> wellInOutSide
            (
                const pointField& samplePts,
                const scalarField& testDistSqr,
                bool testForInside
            ) const;

            //- Check if point is inside surfaces to conform to by at least
            //  testDistSqr
            Field<bool> wellInside
            (
                const pointField& samplePts,
                const scalarField& testDistSqr
            ) const;

            bool wellInside
            (
                const point& samplePt,
                scalar testDistSqr
            ) const;

            //- Check if point is outside surfaces to conform to by at least
            //  testDistSqr
            Field<bool> wellOutside
            (
                const pointField& samplePts,
                const scalarField& testDistSqr
            ) const;

            bool wellOutside
            (
                const point& samplePt,
                scalar testDistSqr
            ) const;

            // Finding if the line joining start and end intersects the surface
            bool findSurfaceAnyIntersection
            (
                const point& start,
                const point& end
            ) const;

            //- Finding if the line joining start and end intersects the surface
            //  and returning the hit and surface information
            void findSurfaceAnyIntersection
            (
                const point& start,
                const point& end,
                pointIndexHit& surfHit,
                label& hitSurface
            ) const;

            //- Finding the nearestIntersection of the surface to start
            void findSurfaceNearestIntersection
            (
                const point& start,
                const point& end,
                pointIndexHit& surfHit,
                label& hitSurface
            ) const;

                //- Find the nearest point to the sample and return it to the
            //  pointIndexHit
            void findSurfaceNearest
            (
                const point& sample,
                scalar nearestDistSqr,
                pointIndexHit& surfHit,
                label& hitSurface
            ) const;

            void findSurfaceNearest
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>& hitInfo,
                labelList& hitSurfaces
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

            //- Find the nearest point on each type of feature edge
            void findEdgeNearestByType
            (
                const point& sample,
                scalar nearestDistSqr,
                List<pointIndexHit>& edgeHit,
                List<label>& featureHit
            ) const;

            //- Find which patch is intersected by the line from one point to
            //  another
            label findPatch(const point& ptA, const point& ptB) const;

            //- Find which patch is closest to the point
            label findPatch(const point& pt) const;


        // Write

            //- Write all components of all the extendedFeatureEdgeMeshes as
            //  an obj file
            void writeFeatureObj(const fileName& prefix) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "conformationSurfacesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
