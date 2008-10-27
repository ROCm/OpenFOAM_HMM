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
    Foam::keyedSurface

Description
    A surface geometry mesh, in which the patch information is conveyed by
    the 'key' attached to each face.

    This form of surface description is particularly useful for reading in
    surface meshes from third-party formats (eg, obj, stl, gts, etc.). It
    can also be particularly useful for situations in which the surface
    many be adjusted in an arbitrary manner without worrying about needed
    to adjust the patch information (eg, surface refinement).

See Also
    The Foam::meshedSurface - which is organized as a surface mesh, but
    with independent patch information.

SourceFiles
    keyedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef keyedSurface_H
#define keyedSurface_H

#include "pointField.H"
#include "PrimitivePatchExtra.H"
#include "boolList.H"
#include "geometricSurfacePatchList.H"
#include "surfacePatchList.H"
#include "face.H"
#include "Keyed.H"
#include "xfer.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;
class IFstream;
class meshedSurface;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                      Class keyedSurface Declaration
\*---------------------------------------------------------------------------*/

class keyedSurface
:
    public PrimitivePatchExtra<Keyed<face>, ::Foam::List, pointField, point>
{
    friend class meshedSurface;

protected:

    // Protected Member Data

        //- Typedef if this type has not already been defined
        typedef Keyed<face> keyedFace;

        //- Typedef for similar code in keyedSurface and meshedSurface
        typedef keyedFace FaceType;

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

    // Private data

        //- Patch information (face ordering nFaces/startFace only used
        //  during reading and writing)
        geometricSurfacePatchList geoPatches_;

    // Private member functions

        //- Read in Foam format
        bool read(Istream&);

        //- Generic read routine. Chooses reader based on extension.
        //? bool read(const fileName&, const word& ext);

protected:

    // Static protected functions

        //- Read non-comment line
        static string getLineNoComment(IFstream&);

    // Protected Member functions

        //- Sets default patch names based on the maximum patch number
        void setPatches(const label maxPatch);

        //- Finds maximum patch number and sets default patch names
        void setPatches();

        //- Sets patch names from hashed values (id -> name)
        void setPatches
        (
            const Map<word>& regionNames,
            const label maxPatch = -1
        );

        //- Sets patch names from hashed values (name -> id)
        void setPatches(const HashTable<label>& groupToPatch);

        //- Join the faces by removing duplicate points.
        //  Returns true if any points merged
        bool stitchFaces(const scalar tol=SMALL, const bool verbose=false);


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
        ClassName("keyedSurface");

        //- Return the default geometric patch type (usually "empty")
        static word defaultGeometricType;

    // Static

        //- Name of keyedSurface directory to use.
        static fileName triSurfInstance(const Time&);

        //- Name of keyedSurface directory to use.
        static fileName triSurfName(const Time&);


    // Constructors

        //- Construct null
        keyedSurface();

        //- Construct from components (points, labelled faces, patches).
        keyedSurface
        (
            const pointField&,
            const List<keyedFace>&,
            const geometricSurfacePatchList&
        );

        //- Construct by transferring components
        //  (points, labelled faces, patches.
        keyedSurface
        (
            const xfer<pointField>&,
            const xfer<List<keyedFace> >&,
            const xfer<geometricSurfacePatchList>&
        );

        //- Construct by transferring points and labelled faces
        //  With region names per map or set to default.
        keyedSurface
        (
            const xfer<pointField>&,
            const xfer<List<keyedFace> >&,
            const Map<word>& regionNames = Map<word>::null()
        );

        //- Construct by transferring points and labelled faces
        //  with patch-names from hash
        keyedSurface
        (
            const xfer<pointField>&,
            const xfer<List<keyedFace> >&,
            const HashTable<label>& labelToRegion
        );

        //- Construct by transferring points, copying unlabelled faces.
        //  Set single default patch.
        keyedSurface
        (
            const xfer<pointField>&,
            const List<face>&
        );

        //- Construct from a boundary mesh with local points/faces
        keyedSurface
        (
            const polyBoundaryMesh&,
            const bool globalPoints=false
        );

        //- Construct from a meshedSurface
        keyedSurface(const meshedSurface&);

        //- Construct from file name (uses extension to determine type)
        keyedSurface
        (
            const fileName&,
            const bool triangulate=false
        );

        //- Construct from file name (uses extension to determine type)
        keyedSurface
        (
            const fileName&,
            const word&,
            const bool triangulate=false
        );

        //- Construct from Istream
        keyedSurface(Istream&);

        //- Construct from objectRegistry
        keyedSurface(const Time&);

        //- Construct as copy
        keyedSurface(const keyedSurface&);

        //- Construct by transferring the contents from a keyedSurface
        keyedSurface(const xfer<keyedSurface>&);

        //- Construct by transferring the contents from a meshedSurface
        keyedSurface(const xfer<meshedSurface>&);


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            keyedSurface,
            fileExtension,
            (
                const fileName& fName,
                const bool triangulate
            ),
            (fName, triangulate)
        );

    // Selectors

        //- Select constructed from filename (implicit extension)
        static autoPtr<keyedSurface> New
        (
            const fileName& fName,
            const bool triangulate=false
        );

        //- Select constructed from filename (explicit extension)
        static autoPtr<keyedSurface> New
        (
            const fileName& fName,
            const word& ext,
            const bool triangulate=false
        );

    // Destructor

        ~keyedSurface();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            keyedSurface,
            write,
            fileExtension,
            (
                const fileName& fName,
                const keyedSurface& surf
            ),
            (fName, surf)
        );

        //- Write to file
        static void write(const fileName&, const keyedSurface&);

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
                return faces().size();
            }

            //- The surface size is the number of faces
            label size() const
            {
                return nFaces();
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

            const geometricSurfacePatchList& geoPatches() const
            {
                return geoPatches_;
            }

            geometricSurfacePatchList& geoPatches()
            {
                return geoPatches_;
            }

            //- Determine the sort order from the face regions.
            //  Returns patch list and sets faceMap to index within faceLst
            static surfacePatchList sortedRegions
            (
                const List<keyedFace>& faceLst,
                const Map<word>& patchNames,
                labelList& faceMap
            );

            //- Determine the sort order from the region list.
            //  Returns patch list and sets faceMap to indices within faceLst
            static surfacePatchList sortedRegions
            (
                const List<label>& regionLst,
                const Map<word>& patchNames,
                labelList& faceMap
            );

            //- Sort faces according to region.
            //  Returns patch list and sets faceMap to index within faces()
            surfacePatchList sortedRegions(labelList& faceMap) const;

        // Edit

            //- Move points
            virtual void movePoints(const pointField&);

            //- Scale points. A non-positive factor is ignored
            virtual void scalePoints(const scalar&);

            void checkFaces(const bool verbose);

            //- Remove invalid faces
            void cleanup(const bool verbose);

            //- Return new surface. Returns pointMap, faceMap from
            //  subsetMeshMap
            keyedSurface subsetMesh
            (
                const boolList& include,
                labelList& pointMap,
                labelList& faceMap
            ) const;

            //- Transfer the contents of the argument and annull the argument
            void transfer(keyedSurface&);

            //- Transfer the contents of the argument and annull the argument
            void transfer(meshedSurface&);

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

    // Friend Functions

    // Member operators

        void operator=(const keyedSurface&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const keyedSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
