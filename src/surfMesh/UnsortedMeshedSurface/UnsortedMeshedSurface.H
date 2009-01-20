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

#include "BasicMeshedSurface.H"
#include "surfPatchIdentifierList.H"
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
    public BasicMeshedSurface<Face>,
    public fileFormats::surfaceFormatsCore
{
    friend class MeshedSurface<Face>;

private:

    //- Typedefs for convenience
        typedef BasicMeshedSurface<Face> ParentType;
        typedef MeshedSurface<Face> SiblingType;

    // Private Member Data

        //- The regions associated with the faces
        labelList regions_;

        //- Patch information (face ordering nFaces/startFace only used
        //  during reading and writing)
        List<surfPatchIdentifier> patches_;

    // Private member functions

        //- Disable resize with value
        void resize(const label, const Face&);

        //- Disable setSize with value
        void setSize(const label, const Face&);

        //- Read OpenFOAM Surface format
        bool read(Istream&);

protected:

    // Protected Member functions

        //- Return non-const access to the faces
        List<label>& storedRegions()
        {
            return regions_;
        }

        //- Return non-const access to the patches
        List<surfPatchIdentifier>& storedPatches()
        {
            return patches_;
        }

        //- Set new regions from faceMap
        virtual void remapFaces(const UList<label>& faceMap);

public:

        //- Runtime type information
        TypeName("UnsortedMeshedSurface");

    // Static

        //- Can we read this file format?
        static bool canReadType(const word& ext, const bool verbose=false);

        //- Can we read this file format?
        static bool canRead(const fileName&, const bool verbose=false);

        //- Can we write this file format?
        static bool canWriteType(const word& ext, const bool verbose=false);

        static wordHashSet readTypes();
        static wordHashSet writeTypes();

    // Constructors

        //- Construct null
        UnsortedMeshedSurface();

        //- Construct by transferring components
        //  (points, faces, region ids, patches).
        UnsortedMeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face> >&,
            const Xfer<List<label> >& regionIds,
            const Xfer<surfPatchIdentifierList>&
        );

        //- Construct by transferring points, faces.
        //  Use patch information, or set single default patch
        UnsortedMeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face> >&,
            const UList<label>& patchSizes = UList<label>::null(),
            const UList<word>& patchNames = UList<word>::null()
        );

        //- Construct from a boundary mesh with local points/faces
        UnsortedMeshedSurface
        (
            const polyBoundaryMesh&,
            const bool globalPoints=false
        );

        //- Construct from a meshedSurface
        UnsortedMeshedSurface(const MeshedSurface<Face>&);

        //- Construct by transferring the contents from a UnsortedMeshedSurface
        UnsortedMeshedSurface(const Xfer<UnsortedMeshedSurface<Face> >&);

        //- Construct by transferring the contents from a meshedSurface
        UnsortedMeshedSurface(const Xfer<MeshedSurface<Face> >&);

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
                const fileName& name
            ),
            (name)
        );

    // Selectors

        //- Select constructed from filename (explicit extension)
        static autoPtr<UnsortedMeshedSurface> New
        (
            const fileName&,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<UnsortedMeshedSurface> New(const fileName&);

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
                const fileName& name,
                const UnsortedMeshedSurface<Face>& surf
            ),
            (name, surf)
        );

        //- Write to file
        static void write(const fileName&, const UnsortedMeshedSurface<Face>&);


    // Member Functions

    // Access

        //- The surface size is the number of faces
        label size() const
        {
            return ParentType::size();
        }

        //- Reset size of face and region list
        void setSize(const label);

        //- Return const access to the regions
        const List<label>& regions() const
        {
            return regions_;
        }

        //- Return const access to the patches
        const List<surfPatchIdentifier>& patches() const
        {
            return patches_;
        }

        //- Sort faces according to region.
        //  Returns patch list and sets faceMap to index within faces()
        List<surfGroup> sortedRegions(labelList& faceMap) const;

        //- Set regions to 0 and set a single patch
        //  Optionally with a specific name
        void onePatch(const word& name = word::null);

        //- Set regions and patches
        void setPatches(const surfGroupList&);

        //- Set regions and patches
        void setPatches(const UList<label>& sizes, const UList<word>& names);

        //- Set regions and patches with default names
        void setPatches(const UList<label>& sizes);


    // Edit

        //- Clear all storage
        virtual void clear();

        //- Return new surface.
        //  Returns return pointMap, faceMap from subsetMeshMap
        UnsortedMeshedSurface subsetMesh
        (
            const UList<bool>& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return new surface.
        UnsortedMeshedSurface subsetMesh
        (
            const UList<bool>& include
        ) const;

        //- Transfer components (points, faces, region ids).
        virtual void reset
        (
            const Xfer<pointField>&,
            const Xfer<List<Face> >&,
            const Xfer<List<label> >& regionIds = Xfer<List<label> >::null()
        );

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
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }

        //- Write to database
        void write(const Time&) const;


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
