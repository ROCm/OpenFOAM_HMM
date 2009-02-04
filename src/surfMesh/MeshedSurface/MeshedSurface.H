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
    Foam::MeshedSurface

Description
    A surface geometry mesh with zone information, not to be confused
    with a similarily named surfaceMesh, which actually refers to
    the cell faces of a volume mesh.

   The MeshedSurface is intended for surfaces from a variety of sources.
   - A set of points and faces without any surface zone information.
   - A set of points and faces with randomly ordered zone information.
     This could arise, for example, from reading external file formats
     such as STL, etc.

SourceFiles
    MeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef MeshedSurface_H
#define MeshedSurface_H

#include "BasicMeshedSurface.H"
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Time;
class surfMesh;
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
    public BasicMeshedSurface<Face>,
    public fileFormats::surfaceFormatsCore
{
    // friends despite different faces
    template<class Face2>
    friend class MeshedSurface;

    // friends despite different faces
    template<class Face2>
    friend class UnsortedMeshedSurface;

    friend class surfMesh;

private:

    //- Private typedefs for convenience
        typedef BasicMeshedSurface<Face> ParentType;
        typedef UnsortedMeshedSurface<Face>  SiblingType;

    // Private Member Data

        //- Zone information
        // (face ordering nFaces/startFace only used during reading/writing)
        List<surfZone> zones_;

    // Private member functions

        //- Read OpenFOAM Surface format
        bool read(Istream&);

        //- Transfer points/zones and transcribe face -> triFace
        void transcribe(MeshedSurface<face>&);

protected:

    // Protected Member functions

        //- basic sanity check on zones
        void checkZones();

        //- Non-const access to the zones
        surfZoneList& storedZones()
        {
            return zones_;
        }

        //- sort faces by zones and store sorted faces
        void sortFacesAndStore
        (
            const Xfer<List<Face> >& unsortedFaces,
            const Xfer<List<label> >& zoneIds,
            const bool sorted
        );

        //- Set new zones from faceMap
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

        //- Construct by transferring components (points, faces, zones).
        MeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face> >&,
            const Xfer<surfZoneList>&
        );

        //- Construct by transferring points, faces.
        //  Use zone information, or set single default zone.
        MeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face> >&,
            const UList<label>& zoneSizes = UList<label>::null(),
            const UList<word>& zoneNames = UList<word>::null()
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
        MeshedSurface(const Xfer<UnsortedMeshedSurface<Face> >&);

        //- Construct by transferring the contents from a MeshedSurface
        MeshedSurface(const Xfer<MeshedSurface<Face> >&);

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

        const List<surfZone>& zones() const
        {
            return zones_;
        }

        //- set a single zone, optionally with a specific name
        void oneZone(const word& name = word::null);

        //- Add zones
        void addZones
        (
            const UList<surfZone>&,
            const bool cullEmpty=false
        );

        //- Add zones
        void addZones
        (
            const UList<label>& sizes,
            const UList<word>& names,
            const bool cullEmpty=false
        );

        //- Add zones
        void addZones
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
            const labelHashSet& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return new surface.
        MeshedSurface subsetMesh
        (
            const labelHashSet& include
        ) const;

        //- Transfer the contents of the argument and annull the argument
        void transfer(MeshedSurface<Face>&);

        //- Transfer the contents of the argument and annull the argument
        void transfer(UnsortedMeshedSurface<Face>&);

        //- Transfer contents to the Xfer container
        Xfer< MeshedSurface<Face> > xfer();

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
