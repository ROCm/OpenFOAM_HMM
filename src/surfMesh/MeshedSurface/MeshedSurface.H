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
    Foam::MeshedSurface

Description
    A surface geometry mesh with patch information, not to be confused
    with a similarily named surfaceMesh, which actually refers to
    the cell faces of a volume mesh!

   The MeshedSurface is intended to surfaces from a variety of sources.
   - A set of points and faces without any patch information.
   - A set of points and faces with randomly sorted patch information.
     This could arise, for example, from reading external file formats
     such as STL, etc.

SourceFiles
    MeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef MeshedSurface_H
#define MeshedSurface_H

#include "PrimitiveMeshedSurface.H"
#include "surfGroupList.H"
#include "surfaceFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Time;
template<class Face> class MeshedSurface;
template<class Face> class UnsortedMeshedSurface;

class polyBoundaryMesh;
class surfMesh;

template<class Face>
    Ostream& operator<<(Ostream&, const MeshedSurface<Face>&);

/*---------------------------------------------------------------------------*\
                      Class MeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class MeshedSurface
:
    public PrimitiveMeshedSurface<Face>,
    public fileFormats::surfaceFormatsCore
{
    friend class UnsortedMeshedSurface<Face>;

private:

    //- Private typedefs for convenience
        typedef PrimitiveMeshedSurface<Face> ParentType;

    // Private Member Data

        //- Patch information (face ordering nFaces/startFace only used
        //  during reading and writing)
        List<surfGroup> patches_;

    // Private member functions

        //- set a single patch
        void onePatch();

        //- basic sanity check on patches
        void checkPatches();

        //- Sort faces by regionIds and set patches
        void sortFacesByRegion(const UList<label>&, const Map<word>&);

        //- Set new regions from faceMap
        void remapRegions(List<label>& faceMap);

        //- Read OpenFOAM Surface format
        bool read(Istream&);

public:

        //- Runtime type information
        ClassName("MeshedSurface");

    // Static

        //- Can we read this file format?
        static bool canRead(const fileName&, const bool verbose=false);

        //- Can we read this file format?
        static bool canReadType(const word& ext, const bool verbose=false);

        //- Can we write this file format?
        static bool canWriteType(const word& ext, const bool verbose=false);


    // Constructors

        //- Construct null
        MeshedSurface();

        //- Construct by transferring components
        //  (points, faces and patches).
        MeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const xfer<surfGroupList>&
        );

        //- Construct from points, faces, and patch information
        MeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const UList<label>& patchSizes,
            const UList<word>& patchNames,
            const UList<word>& patchTypes
        );

        //- Construct by transferring points, faces.
        //  Sort faces according to regionIds with the names of the regions,
        //  or jsut set a single default patch.
        MeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const UList<label>& regionIds = UList<label>::null(),
            const Map<word>& regionNames = Map<word>::null()
        );

        //- Construct by transferring points, faces.
        //  Sort faces according to regionIds with patch-names from hash
        MeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const UList<label>& regionIds,
            const HashTable<label>& nameToRegionMapping
        );

        //- Construct from a boundary mesh with local points/faces
        MeshedSurface
        (
            const polyBoundaryMesh&,
            const bool globalPoints=false
        );

        //- Construct from a UnsortedMeshedSurface
        MeshedSurface(const UnsortedMeshedSurface<Face>&);

        //- Construct from a surfMesh
        MeshedSurface(const surfMesh&);

        //- Construct by transferring the contents from a UnsortedMeshedSurface
        MeshedSurface(const xfer<UnsortedMeshedSurface<Face> >&);

        //- Construct by transferring the contents from a MeshedSurface
        MeshedSurface(const xfer<MeshedSurface<Face> >&);

        //- Construct from file name (uses extension to determine type)
        MeshedSurface(const fileName&);

        //- Construct from file name (uses extension to determine type)
        MeshedSurface(const fileName&, const word& ext);

        //- Construct from Istream
        MeshedSurface(Istream&);

        //- Construct from objectRegistry
        MeshedSurface(const Time&);

        //- Construct as copy
        MeshedSurface(const MeshedSurface&);

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            MeshedSurface,
            fileExtension,
            (
                const fileName& fName
            ),
            (fName)
        );

    // Selectors

        //- Select constructed from filename (explicit extension)
        static autoPtr<MeshedSurface> New
        (
            const fileName&,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<MeshedSurface> New(const fileName&);

    // Destructor

        virtual ~MeshedSurface();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            MeshedSurface,
            write,
            fileExtension,
            (
                const fileName& fName,
                const MeshedSurface<Face>& surf
            ),
            (fName, surf)
        );

        //- Write to file
        static void write(const fileName&, const MeshedSurface<Face>&);


    // Member Functions

    // Access

        //- The surface size is the number of faces
        label size() const
        {
            return ParentType::size();
        }

        const List<surfGroup>& patches() const
        {
            return patches_;
        }

        void addPatches(const UList<surfGroup>&);


    // Edit

        //- Clear all storage
        virtual void clear();

        //- Remove invalid faces
        void cleanup(const bool verbose);

        //- Check/fix duplicate/degenerate faces
        virtual bool checkFaces(const bool verbose);

        //- Join the faces by removing duplicate points.
        //  Returns true if any points merged
        virtual bool stitchFaces
        (
            const scalar tol=SMALL,
            const bool verbose=false
        );

        //- Triangulate the surface, return the number of added faces.
        //  The patch list will be adjusted accordingly.
        virtual label triangulate();

        //- Return new surface.
        //  Returns return pointMap, faceMap from subsetMeshMap
        MeshedSurface subsetMesh
        (
            const UList<bool>& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return new surface.
        MeshedSurface subsetMesh
        (
            const UList<bool>& include
        ) const;

        //- Transfer the contents of the argument and annull the argument
        void transfer(MeshedSurface<Face>&);

        //- Transfer the contents of the argument and annull the argument
        void transfer(UnsortedMeshedSurface<Face>&);


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


    // Member operators

        void operator=(const MeshedSurface<Face>&);

    // Ostream Operator

        friend Ostream& operator<<
        <Face>
        (
            Ostream&,
            const MeshedSurface<Face>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialization for holding triangulated information
template<>
inline label MeshedSurface<triFace>::triangulate()
{
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
