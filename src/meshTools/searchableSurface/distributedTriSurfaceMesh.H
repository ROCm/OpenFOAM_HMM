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
    Foam::distributedTriSurfaceMesh

Description
    IOoject and searching on triSurface

SourceFiles
    distributedTriSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef distributedTriSurfaceMesh_H
#define distributedTriSurfaceMesh_H

#include "triSurfaceMesh.H"
#include "IOdictionary.H"
#include "Pair.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistribute;

// Typedefs
typedef Pair<point> segment;
template<>
inline bool contiguous<segment>() {return contiguous<point>();}


/*---------------------------------------------------------------------------*\
                           Class distributedTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class distributedTriSurfaceMesh
:
    public triSurfaceMesh
{
private:

    // Static data

        //- Merging distance
        static scalar mergeDist_;


    // Private member data

        //- Bounding box settings
        IOdictionary dict_;

        //- bounding boxes of all processors
        List<List<treeBoundBox> > procBb_;

        //- Global triangle numbering
        mutable autoPtr<globalIndex> globalTris_;


    // Private Member Functions

        // Read

            //- Read my additional data
            bool read();


        // Line intersection

            static bool isLocal
            (
                const List<treeBoundBox>& myBbs,
                const point& start,
                const point& end
            );

            //- Split segment into subsegments overlapping the processor
            //  bounding boxes. Sort per processor.
            void splitSegment
            (
                const label,
                const point& start,
                const point& end,

                DynamicList<segment>&,
                DynamicList<label>&,
                List<DynamicList<label> >&
            ) const;

            //- Divide edges into local and remote segments. Construct map to
            //  distribute and collect data.
            autoPtr<mapDistribute> constructSegments
            (
                const pointField& start,
                const pointField& end,

                List<segment>& allSegments,
                List<label>& allSegmentMap
            ) const;

            //- Split edges, distribute, test and collect.
            void findLine
            (
                const bool nearestIntersection,
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>& info
            ) const;


        // Triangle index

            //- Obtains global indices from pointIndexHit and swaps them back
            //  to their original processor. Used to calculate local region
            //  and normal.
            autoPtr<mapDistribute> calcLocalQueries
            (
                const List<pointIndexHit>&,
                labelList& triangleIndex
            ) const;


        // Nearest

            label calcOverlappingProcs
            (
                const point& centre,
                const scalar radiusSqr,
                boolList& overlaps
            ) const;

            autoPtr<mapDistribute> calcLocalQueries
            (
                const pointField& centres,
                const scalarField& radiusSqr,

                pointField& allCentres,
                scalarField& allRadiusSqr,
                labelList& allSegmentMap
            ) const;


        // Surface redistribution

            //- Calculate surface bounding box
            void calcBounds(boundBox& bb, label& nPoints) const;

            //- Does any part of triangle overlap bb.
            static bool overlaps
            (
                const List<treeBoundBox>& bb,
                const point& p0,
                const point& p1,
                const point& p2
            );

            //- Find points used in subset
            static void subsetMeshMap
            (
                const triSurface& s,
                const boolList& include,
                const label nIncluded,
                labelList& newToOldPoints,
                labelList& oldToNewPoints,
                labelList& newToOldFaces
            );

            //- Construct subsetted surface
            static triSurface subsetMesh
            (
                const triSurface& s,
                const labelList& newToOldPoints,
                const labelList& oldToNewPoints,
                const labelList& newToOldFaces
            );

            //- Subset given marked faces
            static triSurface subsetMesh
            (
                const triSurface& s,
                const boolList& include,
                labelList& newToOldPoints,
                labelList& newToOldFaces
            );

            //- Subset given marked faces
            static triSurface subsetMesh
            (
                const triSurface& s,
                const labelList& newToOldFaces,
                labelList& newToOldPoints
            );

            //- Find triangle otherF in allFaces.
            static label findTriangle
            (
                const List<labelledTri>& allFaces,
                const labelListList& allPointFaces,
                const labelledTri& otherF
            );

            //- Merge triSurface (subTris, subPoints) into allTris, allPoints.
            static void merge
            (
                const scalar mergeDist,
                const List<labelledTri>& subTris,
                const pointField& subPoints,

                List<labelledTri>& allTris,
                pointField& allPoints,

                labelList& faceConstructMap,
                labelList& pointConstructMap
            );

            //- Distribute stored fields
            template<class Type>
            void distributeFields(const mapDistribute& map);


        //- Disallow default bitwise copy construct
        distributedTriSurfaceMesh(const distributedTriSurfaceMesh&);

        //- Disallow default bitwise assignment
        void operator=(const distributedTriSurfaceMesh&);


public:

    //- Runtime type information
    TypeName("distributedTriSurfaceMesh");


    // Constructors

        //- Construct from triSurface
        distributedTriSurfaceMesh
        (
            const IOobject&,
            const triSurface&,
            const dictionary& dict
        );

        //- Construct read
        distributedTriSurfaceMesh(const IOobject& io);

        //- Construct from dictionary (used by searchableSurface)
        distributedTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


    // Destructor

        virtual ~distributedTriSurfaceMesh();

        //- Clear storage
        void clearOut();


    // Member Functions

        //- Return merge distance
        static scalar mergeDist()
        {
            return mergeDist_;
        }

        //- Set the merge distance, returning the previous value
        static scalar setMergeDist(const scalar t)
        {
            if (t < -VSMALL)
            {
                FatalErrorIn
                (
                    "scalar distributedTriSurfaceMesh::setMergeDist"
                    "(const scalar t)"
                )   << "Negative merge distance tolerance. This is not allowed."
                    << abort(FatalError);
            }

            scalar oldTol = mergeDist_;
            mergeDist_ = t;

            return oldTol;
        }


        //- Triangle indexing (demand driven)
        const globalIndex& globalTris() const;


        // searchableSurface implementation

            //- Whether supports volume type below. I.e. whether is closed.
            //  Not supported.
            virtual bool hasVolumeType() const
            {
                return false;
            }

            //- Range of global indices that can be returned.
            virtual label globalSize() const
            {
                const labelList& offsets = globalTris().offsets();
                return offsets[offsets.size()-1];
            }

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit> >&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;

            //- Set bounds of surface. Bounds currently set as list of
            //  bounding boxes. Will do redistribution of surface to locally
            //  have all triangles overlapping bounds.
            //  Larger bounds: more triangles (memory), more fully local tests
            //  (quick).
            //  keepNonLocal = true : keep triangles that do not overlap
            //  any processor bounds.
            //  Should really be split into a routine to determine decomposition
            //  and one that does actual distribution but determining
            //  decomposition with duplicate triangle merging requires
            //  same amoun as work as actual distribution.
            virtual void distribute
            (
                const List<treeBoundBox>&,
                const bool keepNonLocal,
                autoPtr<mapDistribute>& faceMap,
                autoPtr<mapDistribute>& pointMap
            );


        // Other

            //- Specific to triSurfaceMesh: from a set of hits (points and
            //  indices) get the specified field. Misses do not get set.
            virtual void getField
            (
                const word& fieldName,
                const List<pointIndexHit>&,
                labelList& values
            ) const;

            //- Subset the part of surface that is overlapping bounds.
            static triSurface overlappingSurface
            (
                const triSurface&,
                const List<treeBoundBox>&,
                labelList& subPointMap,
                labelList& subFaceMap
            );

            //- Print some stats. Parallel aware version of
            //  triSurface::writeStats.
            void writeStats(Ostream& os) const;


        // regIOobject implementation

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "distributedTriSurfaceMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
