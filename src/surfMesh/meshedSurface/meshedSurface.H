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
    Foam::meshedSurface

Description
    A surface geometry mesh with patch information, not to be confused
    with a similarily named surfaceMesh, which actually refers to
    the cell faces of a volume mesh!

   The meshedSurface is intended to surfaces from a variety of sources.
   - A set of points and faces without any patch information.
   - A set of points and faces with randomly sorted patch information.
     This could arise, for example, from reading external file formats
     such as STL, etc.

SourceFiles
    meshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef meshedSurface_H
#define meshedSurface_H

#include "pointField.H"
#include "PrimitivePatchExtra.H"
#include "boolList.H"
#include "surfGroupList.H"
#include "face.H"
#include "Keyed.H"
#include "xfer.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Time;
class keyedSurface;
class meshedSurface;
class polyBoundaryMesh;
class surfMesh;

Istream& operator>>(Istream&, meshedSurface&);
Ostream& operator<<(Ostream&, const meshedSurface&);

/*---------------------------------------------------------------------------*\
                      Class meshedSurface Declaration
\*---------------------------------------------------------------------------*/

class meshedSurface
:
    public PrimitivePatchExtra<face, ::Foam::List, pointField, point>
{
    friend class keyedSurface;

protected:

    // Protected Member Data

        //- Typedef if this type has not already been defined
        typedef Keyed<face> keyedFace;

        //- Typedef for similar code in keyedSurface and meshedSurface
        typedef face FaceType;

private:

    // Private Member Data

    // Private typedefs

        typedef PrimitivePatchExtra
        <
            FaceType,
            ::Foam::List,
            pointField,
            point
        >
        MeshStorage;

    // Private Member Data

        //- Patch information (face ordering nFaces/startFace only used
        //  during reading and writing)
        List<surfGroup> patches_;

    // Private member functions

        //- Sort faces by regionIds and set patches
        void sortFacesByRegion(const List<label>&, const Map<word>&);

        //- set a single patch
        void onePatch();

        //- basic sanity check on patches
        void checkPatches();

        //- Read OpenFOAM Surface format
        bool read(Istream&, const bool triangulate=false);


protected:

    // Protected Member functions

        //- Return non-const access to global points
        pointField& points()
        {
            return const_cast<pointField&>(MeshStorage::points());
        }

        //- Return non-const access to the faces
        List<FaceType>& faces()
        {
            return static_cast<List<FaceType> &>(*this);
        }

public:

        //- Runtime type information
        ClassName("meshedSurface");

    // Static

        //- Name of meshedSurface directory to use.
        static fileName triSurfInstance(const Time&);

        //- Name of meshedSurface directory to use.
        static fileName triSurfName(const Time&);


    // Constructors

        //- Construct null
        meshedSurface();

        //- Construct by transferring components
        //  (points, faces and patches).
        meshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<face> >&,
            const xfer<surfGroupList>&
        );

        //- Construct from points, faces, and patch information
        meshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<face> >&,
            const List<label>& patchSizes,
            const List<word>& patchNames,
            const List<word>& patchTypes
        );

        //- Construct by transferring points, faces.
        //  Sort faces according to regionIds with the names of the regions,
        //  or jsut set a single default patch.
        meshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<face> >&,
            const List<label>& regionIds = List<label>::null(),
            const Map<word>& regionNames = Map<word>::null()
        );

        //- Construct by transferring points, faces.
        //  Sort faces according to regionIds with patch-names from hash
        meshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<face> >&,
            const List<label>& regionIds,
            const HashTable<label>& nameToRegionMapping
        );

        //- Construct from a boundary mesh with local points/faces
        meshedSurface
        (
            const polyBoundaryMesh&,
            const bool globalPoints=false
        );

        //- Construct from a keyedSurface
        meshedSurface(const keyedSurface&);

        //- Construct from a surfMesh
        meshedSurface(const surfMesh&);

        //- Construct by transferring the contents from a keyedSurface
        meshedSurface(const xfer<keyedSurface>&);

        //- Construct by transferring the contents from a meshedSurface
        meshedSurface(const xfer<meshedSurface>&);


        //- Construct from file name (uses extension to determine type)
        meshedSurface
        (
            const fileName&,
            const bool triangulate=false
        );

        //- Construct from file name (uses extension to determine type)
        meshedSurface
        (
            const fileName&,
            const word& ext,
            const bool triangulate=false
        );

        //- Construct from Istream
        meshedSurface(Istream&, const bool triangulate=false);

        //- Construct from objectRegistry
        meshedSurface(const Time&);

        //- Construct as copy
        meshedSurface(const meshedSurface&);

    // Selectors

        //- Select constructed from filename (implicit extension)
        static autoPtr<meshedSurface> New
        (
            const fileName&,
            const bool triangulate=false
        );

        //- Select constructed from filename (explicit extension)
        static autoPtr<meshedSurface> New
        (
            const fileName&,
            const word&,
            const bool triangulate=false
        );


    // Destructor

        virtual ~meshedSurface();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            meshedSurface,
            write,
            fileExtension,
            (
                const fileName& fName,
                const meshedSurface& surf
            ),
            (fName, surf)
        );

        //- Write to file
        static void write(const fileName&, const meshedSurface&);

        //- Can we read this file format?
        static bool canRead(const word& ext, const bool verbose=false);

        //- Can we write this file format?
        static bool canWrite(const word& ext, const bool verbose=false);


    // Member Functions

    // Access

        //- Return the number of points
        label nPoints() const
        {
            return points().size();
        }

        //- Return the number of faces
        label nFaces() const
        {
            return MeshStorage::size();
        }

        //- The surface size is the number of faces
        label size() const
        {
            return MeshStorage::size();
        }

        //- Return const access to global points
        const pointField& points() const
        {
            return MeshStorage::points();
        }

        //- Return const access to the faces
        const List<FaceType>& faces() const
        {
            return static_cast<const List<FaceType> &>(*this);
        }

        const surfGroupList& patches() const
        {
            return patches_;
        }

        surfGroupList& patches()
        {
            return patches_;
        }

    // Edit

        //- Move points
        virtual void movePoints(const pointField&);

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar&);

        //- Triangulate the surface, return the number of added faces.
        //  The patch list will be adjusted accordingly.
        label triangulate();

        //- Join the faces by removing duplicate points.
        //  Returns true if any points merged
        bool stitchFaces(const scalar tol=SMALL, const bool verbose=false);

        //- Check/fix duplicate/degenerate faces
        void checkFaces(const bool verbose);

        //- Remove invalid faces
        void cleanup(const bool verbose);

        //- Return new surface.
        // Returns pointMap, faceMap from subsetMeshMap
        meshedSurface subsetMesh
        (
            const boolList& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Transfer the contents of the argument and annull the argument
        void transfer(meshedSurface&);

        //- Transfer the contents of the argument and annull the argument
        void transfer(keyedSurface&);


    // Read

        //- Read from file. Chooses reader based on extension
        bool read
        (
            const fileName&,
            const word& ext,
            const bool triangulate=false
        );


    // Write

        //- Write to Ostream in simple FOAM format
        virtual void write(Ostream&) const;

        //- Generic write routine. Chooses writer based on extension.
        virtual void write(const fileName& fName) const
        {
            write(fName, *this);
        }

        //- Write to database
        void write(const Time&) const;

        //- Write some statistics
        void writeStats(Ostream&) const;


    // Member operators

        void operator=(const meshedSurface&);

    // Ostream Operator

        // friend Istream& operator>>(Istream&, meshedSurface&);
        friend Ostream& operator<<(Ostream&, const meshedSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
