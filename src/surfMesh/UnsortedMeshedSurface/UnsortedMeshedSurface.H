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
    Foam::UnsortedMeshedSurface

Description
    A surface geometry mesh, in which the patch information is conveyed by
    the 'region' associated with each face.

    This form of surface description is particularly useful for reading in
    surface meshes from third-party formats (eg, obj, stl, gts, etc.). It
    can also be particularly useful for situations in which the surface
    many be adjusted in an arbitrary manner without worrying about needed
    to adjust the patch information (eg, surface refinement).

See Also
    The Foam::meshedSurface - which is organized as a surface mesh, but
    with independent patch information.

SourceFiles
    UnsortedMeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef UnsortedMeshedSurface_H
#define UnsortedMeshedSurface_H

#include "pointField.H"
#include "PrimitivePatchExtra.H"
#include "boolList.H"
#include "surfPatchIdentifierList.H"
#include "surfGroupList.H"
#include "surfaceFormatsCore.H"
#include "xfer.H"
#include "face.H"
#include "triFace.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Time;
class IFstream;
template<class Face> class UnsortedMeshedSurface;
template<class Face> class MeshedSurface;

class polyBoundaryMesh;

template<class Face>
    Ostream& operator<<(Ostream&, const UnsortedMeshedSurface<Face>&);


/*---------------------------------------------------------------------------*\
                      Class UnsortedMeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class UnsortedMeshedSurface
:
    public PrimitivePatchExtra<Face, ::Foam::List, pointField, point>,
    public fileFormats::surfaceFormatsCore
{
    friend class MeshedSurface<Face>;

protected:

    // Protected Member Data

        //- Typedef for type holding the region (patch) informationm
        typedef surfPatchIdentifier PatchRegionType;

private:

    //- Private typedefs for convenience
        typedef UnsortedMeshedSurface<Face> ThisType;
        typedef PrimitivePatchExtra
        <
            Face,
            ::Foam::List,
            pointField,
            point
        >
        ParentType;

    // Private Member Data

        //- The regions associated with the faces
        labelList regions_;

        //- Patch information (face ordering nFaces/startFace only used
        //  during reading and writing)
        List<PatchRegionType> patches_;

    // Private member functions

        //- Disable setSize with value
        void setSize(const label, const Face&);

        //- Read OpenFOAM Surface format
        bool read(Istream&);

protected:

    // Static protected functions

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

        //- Return non-const access to global points
        pointField& points()
        {
            return const_cast<pointField&>(ParentType::points());
        }

        //- Return non-const access to the faces
        List<Face>& faces()
        {
            return static_cast<List<Face> &>(*this);
        }

        //- Return non-const access to the region Ids
        List<label>& regions()
        {
            return regions_;
        }

public:

        //- Runtime type information
        TypeName("UnsortedMeshedSurface");

    // Static

        //- Can we read this file format?
        static bool canRead(const word& ext, const bool verbose=false);

        //- Can we write this file format?
        static bool canWrite(const word& ext, const bool verbose=false);

    // Constructors

        //- Construct null
        UnsortedMeshedSurface();

        //- Construct by transferring components
        //  (points, faces, region ids, patches).
        UnsortedMeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const xfer<List<label> >& regionIds,
            const xfer<surfPatchIdentifierList>&
        );

        //- Construct by transferring points, faces and regionIds
        //  with region names per map or set to default.
        UnsortedMeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const xfer<List<label> >& regionIds,
            const Map<word>& regionNames = Map<word>::null()
        );

        //- Construct by transferring points, faces and regionIds
        //  with patch-names from hash
        UnsortedMeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const xfer<List<label> >& regionIds,
            const HashTable<label>& labelToRegion
        );

        //- Construct by transferring points, faces.
        //  Set single default patch.
        UnsortedMeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&
        );

#if 0
        //- Construct from a boundary mesh with local points/faces
        UnsortedMeshedSurface
        (
            const polyBoundaryMesh&,
            const bool globalPoints=false
        );
#endif

        //- Construct from a meshedSurface
        UnsortedMeshedSurface(const MeshedSurface<Face>&);

        //- Construct by transferring the contents from a UnsortedMeshedSurface
        UnsortedMeshedSurface(const xfer<UnsortedMeshedSurface<Face> >&);

        //- Construct by transferring the contents from a meshedSurface
        UnsortedMeshedSurface(const xfer<MeshedSurface<Face> >&);

        //- Construct from file name (uses extension to determine type)
        UnsortedMeshedSurface(const fileName&);

        //- Construct from file name (uses extension to determine type)
        UnsortedMeshedSurface(const fileName&, const word&);

        //- Construct from Istream
        UnsortedMeshedSurface(Istream&);

        //- Construct from objectRegistry
        UnsortedMeshedSurface(const Time&);

        //- Construct as copy
        UnsortedMeshedSurface(const UnsortedMeshedSurface<Face>&);

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            UnsortedMeshedSurface,
            fileExtension,
            (
                const fileName& fName
            ),
            (fName)
        );

    // Selectors

        //- Select constructed from filename (implicit extension)
        static autoPtr<UnsortedMeshedSurface> New
        (
            const fileName&
        );

        //- Select constructed from filename (explicit extension)
        static autoPtr<UnsortedMeshedSurface> New
        (
            const fileName&,
            const word& ext
        );

    // Destructor

        virtual ~UnsortedMeshedSurface();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            UnsortedMeshedSurface,
            write,
            fileExtension,
            (
                const fileName& fName,
                const UnsortedMeshedSurface<Face>& surf
            ),
            (fName, surf)
        );

        //- Write to file
        static void write(const fileName&, const UnsortedMeshedSurface<Face>&);


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
            return ParentType::size();
        }

        //- The surface size is the number of faces
        label size() const
        {
            return ParentType::size();
        }

        //- Reset size of face and region list
        void setSize(const label);

        //- Return const access to global points
        const pointField& points() const
        {
            return ParentType::points();
        }

        //- Return const access to the faces
        const List<Face>& faces() const
        {
            return static_cast<const List<Face> &>(*this);
        }

        //- Return const access to the regions
        const List<label>& regions() const
        {
            return regions_;
        }

        const surfPatchIdentifierList& patches() const
        {
            return patches_;
        }

        surfPatchIdentifierList& patches()
        {
            return patches_;
        }

        //- Sort faces according to region.
        //  Returns patch list and sets faceMap to index within faces()
        surfGroupList sortedRegions(labelList& faceMap) const;

    // Edit

        //- Clear all storage
        virtual void clear();

        //- Move points
        virtual void movePoints(const pointField&);

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar&);

        //- Remove invalid faces
        void cleanup(const bool verbose);

        //- Check/fix duplicate/degenerate faces
        void checkFaces(const bool verbose);

        //- Join the faces by removing duplicate points.
        //  Returns true if any points merged
        bool stitchFaces(const scalar tol=SMALL, const bool verbose=false);

        //- Triangulate the surface, return the number of added faces.
        label triangulate();

        //- Return new surface. Returns pointMap, faceMap from
        //  subsetMeshMap
        UnsortedMeshedSurface subsetMesh
        (
            const UList<bool>& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Transfer the contents of the argument and annull the argument
        void transfer(UnsortedMeshedSurface<Face>&);

        //- Transfer the contents of the argument and annull the argument
        void transfer(MeshedSurface<Face>&);

    // Read

        //- Read from file. Chooses reader based on explicit extension
        bool read(const fileName&, const word& ext);

        //- Read from file. Chooses reader based on detected extension
        virtual bool read(const fileName&);

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

        void operator=(const UnsortedMeshedSurface<Face>&);


    // Ostream Operator

        friend Ostream& operator<<
        <Face>
        (
            Ostream&,
            const UnsortedMeshedSurface<Face>&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "UnsortedMeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
