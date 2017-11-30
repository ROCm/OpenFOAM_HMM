/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::triSurface

Description
    Triangulated surface description with patch information.

SourceFiles
    triSurface.C

\*---------------------------------------------------------------------------*/

#ifndef triSurface_H
#define triSurface_H

#include "PrimitivePatch.H"
#include "pointField.H"
#include "labelledTri.H"
#include "boolList.H"
#include "geometricSurfacePatchList.H"
#include "surfacePatchList.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class Time;
class IFstream;
class surfZone;
class triSurface;

template<class Face> class MeshedSurface;

Istream& operator>>(Istream&, triSurface&);
Ostream& operator<<(Ostream&, const triSurface&);


/*---------------------------------------------------------------------------*\
                         Class triSurface Declaration
\*---------------------------------------------------------------------------*/

class triSurface
:
    public PrimitivePatch<labelledTri, ::Foam::List, pointField, point>
{
    // Private typedefs

    //- Typedefs for convenience
        typedef labelledTri Face;
        typedef PrimitivePatch
        <
            labelledTri,
            ::Foam::List,
            pointField,
            point
        >
        ParentType;


    // Private data

        //- Patch information
        //  (face ordering nFaces/startFace only used during reading, writing)
        geometricSurfacePatchList patches_;

        static wordHashSet readTypes_;

        static wordHashSet writeTypes_;


    // Demand driven private data.

        //- Edge-face addressing (sorted)
        mutable labelListList* sortedEdgeFacesPtr_;

        //- Label of face that 'owns' edge (i.e. e.vec() is righthanded walk
        //- along face)
        mutable labelList* edgeOwnerPtr_;


    // Private Member Functions

        //- Calculate sorted edgeFaces
        void calcSortedEdgeFaces() const;

        //- Calculate owner
        void calcEdgeOwner() const;

        //- Sort faces according to region.
        //  Returns patch list and sets faceMap to index of labelledTri
        //  inside *this.
        surfacePatchList calcPatches(labelList& faceMap) const;

        //- Sets default values for patches
        void setDefaultPatches();

        //- Function to stitch the triangles by removing duplicate points.
        //  Returns true if any points merged
        bool stitchTriangles
        (
            const scalar tol = SMALL,
            const bool verbose = false
        );

        //- Read in OpenFOAM format
        bool read(Istream& is);

        //- Read in STL format
        bool readSTL(const fileName& filename, bool forceBinary=false);

        //- Generic read routine. Chooses reader based on extension.
        bool read
        (
            const fileName& filename,
            const word& ext,
            const bool check = true
        );

        //- Write STL ASCII format.
        //  Each region becomes a 'solid' 'endsolid' block.
        void writeSTLASCII(const fileName& filename, const bool sort) const;

        //- Write STL BINARY format
        void writeSTLBINARY(const fileName& filename) const;

        //- Write GTS (Gnu Tri Surface library) format.
        void writeGTS(const fileName& filename, const bool sort) const;

        //- Generic write routine. Chooses writer based on extension.
        //  The sort option may not have an effect.
        void write
        (
            const fileName& filename,
            const word& ext,
            const bool sort
        ) const;


    // Static private functions

        //- Convert faces to labelledTri. All get same region.
        static List<labelledTri> convertToTri
        (
            const faceList& faces,
            const label defaultRegion = 0
        );

        //- Convert triFaces to labelledTri. All get same region.
        static List<labelledTri> convertToTri
        (
            const triFaceList& faces,
            const label defaultRegion = 0
        );

        //- Helper function to print triangle info
        static void printTriangle
        (
            Ostream& os,
            const string& pre,
            const labelledTri& f,
            const pointField& points
        );


protected:

    // Protected Member Functions

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return const_cast<pointField&>(ParentType::points());
        }

        //- Non-const access to the faces
        List<Face>& storedFaces()
        {
            return static_cast<List<Face>&>(*this);
        }


public:

    // Public typedefs

        //- Placeholder only, but do not remove - it is needed for GeoMesh
        typedef bool BoundaryMesh;

        //- Runtime type information
        ClassName("triSurface");


    // Static

        //- Name of triSurface directory to use.
        static fileName triSurfInstance(const Time&);

        //- Can we read this file format?
        static bool canRead(const fileName& name, const bool verbose=false);

        //- Can we read this file format?
        static bool canReadType(const word& ext, const bool verbose=false);

        //- Can we write this file format?
        static bool canWriteType(const word& ext, const bool verbose=false);

        //- Known readable file-types
        static const wordHashSet& readTypes();

        //- Known writable file-types
        static const wordHashSet& writeTypes();


    // Constructors

        //- Construct null
        triSurface();

        //- Copy construct
        triSurface(const triSurface& surf);

        //- Move construct
        triSurface(triSurface&& surf);

        //- Construct from triangles, patches, points.
        triSurface
        (
            const List<labelledTri>& triangles,
            const geometricSurfacePatchList& patches,
            const pointField& points
        );

        //- Construct from triangles, patches, points. Reuse storage.
        triSurface
        (
            List<labelledTri>& triangles,
            const geometricSurfacePatchList& patches,
            pointField& points,
            const bool reuse
        );

        //- Construct by transferring (triangles, points) components.
        triSurface
        (
            const Xfer<List<labelledTri>>& triangles,
            const geometricSurfacePatchList& patches,
            const Xfer<List<point>>& points
        );

        //- Construct from triangles, points. Set patch names to default.
        triSurface
        (
            const List<labelledTri>& triangles,
            const pointField& points
        );

        //- Construct from triangles, points Set region to 0 and default
        //- patchName.
        triSurface
        (
            const triFaceList& triangles,
            const pointField& points
        );

        //- Construct from file name (uses extension to determine type).
        //  Optional (positive, non-zero) point scaling is possible.
        triSurface(const fileName& name, const scalar scaleFactor = -1);

        //- Construct from Istream
        triSurface(Istream& is);

        //- Construct from objectRegistry
        triSurface(const Time& d);


    //- Destructor
    virtual ~triSurface();


    // Member Functions

        void clearOut();

        void clearTopology();

        void clearPatchMeshAddr();

        void swap(triSurface& surf);


      // Access

        const geometricSurfacePatchList& patches() const
        {
            return patches_;
        }

        geometricSurfacePatchList& patches()
        {
            return patches_;
        }

        //- Return const access to the faces
        inline const List<labelledTri>& surfFaces() const
        {
            return static_cast<const List<labelledTri>&>(*this);
        }

        //- Return edge-face addressing sorted (for edges with more than
        //  2 faces) according to the angle around the edge.
        //  Orientation is anticlockwise looking from
        //  edge.vec(localPoints())
        const labelListList& sortedEdgeFaces() const;

        //- If 2 face neighbours: label of face where ordering of edge
        //  is consistent with righthand walk.
        //  If 1 neighbour: label of only face.
        //  If >2 neighbours: undetermined.
        const labelList& edgeOwner() const;


        //- Face area vectors (normals)
        inline const vectorField& Sf() const
        {
            return ParentType::faceAreas();
        }

        //- Face area magnitudes
        inline const scalarField& magSf() const
        {
            return ParentType::magFaceAreas();
        }

        //- Face centres
        inline const vectorField& Cf() const
        {
            return ParentType::faceCentres();
        }


      // Interoperability with other surface mesh classes

        //- Sort faces according to zoneIds
        //  Returns a surfZoneList and sets faceMap to index within faces()
        //  (i.e. map from original,unsorted to sorted)
        List<surfZone> sortedZones(labelList& faceMap) const;

        //- Create a list of faces from the triFaces
        void triFaceFaces(List<face>& plainFaceList) const;


      // Edit

        //- Move points
        virtual void movePoints(const pointField& newPoints);

        //- Scale points. A non-positive factor is ignored.
        virtual void scalePoints(const scalar scaleFactor);

        //- Check/remove duplicate/degenerate triangles
        void checkTriangles(const bool verbose);

        //- Check triply (or more) connected edges.
        void checkEdges(const bool verbose);

        //- Remove non-valid triangles
        void cleanup(const bool verbose);

        //- Fill faceZone with currentZone for every face reachable
        //  from facei without crossing edge marked in borderEdge.
        //  Note: faceZone has to be sized nFaces before calling this fun.
        void markZone
        (
            const boolList& borderEdge,
            const label facei,
            const label currentZone,
            labelList& faceZone
        ) const;

        //- (size and) fills faceZone with zone of face. Zone is area
        //  reachable by edge crossing without crossing borderEdge
        //  (bool for every edge in surface). Returns number of zones.
        label markZones
        (
            const boolList& borderEdge,
            labelList& faceZone
        ) const;

        //- 'Create' sub mesh, including only faces for which boolList
        //- entry is true
        //  Sets: pointMap: from new to old localPoints
        //        faceMap: new to old faces
        void subsetMeshMap
        (
            const boolList& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return new surface. Returns pointMap, faceMap from
        //  subsetMeshMap
        triSurface subsetMesh
        (
            const boolList& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;


        //- Transfer stored faces to an Xfer container
        Xfer<List<labelledTri>> xferFaces();

        //- Transfer stored points to an Xfer container
        Xfer<List<point>> xferPoints();

        //- Alter contents by transferring (triangles, points) components.
        //  Patch information is small and therefore just copied.
        void reset
        (
            const Xfer<List<labelledTri>>& triangles,
            const geometricSurfacePatchList& patches,
            const Xfer<List<point>>& points
        );

        //- Alter contents by transferring (triangles, points) components.
        //  Patch information is small and therefore just copied.
        void reset(MeshedSurface<labelledTri>& input);


      // Write

        //- Write to Ostream in simple FOAM format
        void write(Ostream& os) const;

        //- Generic write routine. Chooses writer based on extension.
        void write(const fileName&, const bool sortByRegion = false) const;

        //- Write to database
        void write(const Time& d) const;

        //- Write some statistics
        void writeStats(Ostream& os) const;


    // Member operators

        //- Copy assignment
        void operator=(const triSurface& surf);

        //- Move assignment
        void operator=(triSurface&& surf);


    // IOstream Operators

        friend Istream& operator>>(Istream&, triSurface&);
        friend Ostream& operator<<(Ostream&, const triSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
