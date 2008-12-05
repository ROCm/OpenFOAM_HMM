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
    the cell faces of a volume mesh.

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
#include "HashSet.H"

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
        typedef UnsortedMeshedSurface<Face>  SiblingType;

    // Private Member Data

        //- Patch information (face ordering nFaces/startFace only used
        //  during reading and writing)
        List<surfGroup> patches_;

    // Private member functions

        //- Read OpenFOAM Surface format
        bool read(Istream&);

protected:

    // Protected Member functions

        //- basic sanity check on patches
        void checkPatches();

        //- sort faces by regions and store sorted faces
        void sortFacesAndStore
        (
            const xfer<List<Face> >& unsortedFaces,
            const xfer<List<label> >& regionIds,
            const bool sorted
        );

        //- Set new regions from faceMap
        virtual void remapFaces(const UList<label>& faceMap);

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

        static wordHashSet readTypes();
        static wordHashSet writeTypes();

    // Constructors

        //- Construct null
        MeshedSurface();

        //- Construct by transferring components (points, faces, patches).
        MeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const xfer<surfGroupList>&
        );

        //- Construct by transferring points, faces.
        //  Use patch information, or set single default patch.
        MeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&,
            const UList<label>& patchSizes = UList<label>::null(),
            const UList<word>& patchNames = UList<word>::null()
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
                const fileName& name
            ),
            (name)
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
                const fileName& name,
                const MeshedSurface<Face>& surf
            ),
            (name, surf)
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

        //- set a single patch, optionally with a specific name
        void onePatch(const word& name = word::null);

        //- add patches
        void addPatches
        (
            const UList<surfGroup>&,
            const bool cullEmpty=false
        );

        //- add patches
        void addPatches
        (
            const UList<label>& sizes,
            const UList<word>& names,
            const bool cullEmpty=false
        );

        //- add patches
        void addPatches
        (
            const UList<label>& sizes,
            const bool cullEmpty=false
        );

    // Edit

        //- Clear all storage
        virtual void clear();

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
        virtual void write(const fileName& name) const
        {
            write(name, *this);
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
