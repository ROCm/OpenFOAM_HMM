/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::UnsortedMeshedSurface

Description
    A surface geometry mesh, in which the surface zone information is
    conveyed by the 'zoneId' associated with each face.

    This form of surface description is particularly useful for reading in
    surface meshes from third-party formats (eg, obj, stl, gts, etc.). It
    can also be particularly useful for situations in which the surface
    many be adjusted in an arbitrary manner without worrying about needed
    to adjust the zone information (eg, surface refinement).

See also
    The Foam::MeshedSurface - which is organized as a surface mesh, but
    with independent zone information.

SourceFiles
    UnsortedMeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef UnsortedMeshedSurface_H
#define UnsortedMeshedSurface_H

#include "MeshedSurface.H"
#include "surfZoneIdentifierList.H"
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "HashSet.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
//
class Time;
class Istream;
class Ostream;

template<class Face> class MeshedSurface;
template<class Face> class MeshedSurfaceProxy;
template<class Face> class UnsortedMeshedSurface;

template<class Face>
Istream& operator>>(Istream&, UnsortedMeshedSurface<Face>&);
template<class Face>
Ostream& operator<<(Ostream&, const UnsortedMeshedSurface<Face>&);

/*---------------------------------------------------------------------------*\
                   Class UnsortedMeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class UnsortedMeshedSurface
:
    public MeshedSurface<Face>
{
    // Friends, regardless of face representations
    template<class Face2> friend class MeshedSurface;
    template<class Face2> friend class UnsortedMeshedSurface;

    // Friendship with surfMesh is needed for transferring
    friend class surfMesh;

private:

    // Private typedefs for convenience

        typedef MeshedSurface<Face>       ParentType;
        typedef MeshedSurface<Face>       FriendType;
        typedef MeshedSurfaceProxy<Face>  ProxyType;


    // Private Member Data

        //- The zone Id associated with each face
        labelList zoneIds_;

        //- Zone information (face ordering nFaces/startFace only used
        //  during reading and writing)
        List<surfZoneIdentifier> zoneToc_;


    // Private Member Functions

        //- Disable resize with value
        void resize(const label, const Face&) = delete;

        //- Disable setSize with value
        void setSize(const label, const Face&) = delete;

        //- Read/construct from Istream
        Istream& read(Istream&);

        //- Write to Ostream
        Ostream& write(Ostream&) const;

protected:

    // Protected Member functions

        //- Return non-const access to the zone Ids
        labelList& storedZoneIds()
        {
            return zoneIds_;
        }

        //- Return non-const access to the zone table-of-contents
        List<surfZoneIdentifier>& storedZoneToc()
        {
            return zoneToc_;
        }

        //- Set new zones from faceMap
        virtual void remapFaces(const labelUList& faceMap);


public:

    // Public typedefs

        //- Face type used
        typedef Face FaceType;

        //- Runtime type information
        TypeName("UnsortedMeshedSurface");


    // Static

        //- Can we read this file format type?
        static bool canReadType(const word& ext, bool verbose=false);

        //- Can we read this file format?
        static bool canRead(const fileName& name, bool verbose=false);

        //- Can we write this file format type?
        static bool canWriteType(const word& ext, bool verbose=false);

        //- Known readable file-types
        static wordHashSet readTypes();

        //- Known writable file-types
        static wordHashSet writeTypes();


    // Constructors

        //- Construct null
        UnsortedMeshedSurface();

        //- Construct as copy
        UnsortedMeshedSurface(const UnsortedMeshedSurface<Face>& surf);

        //- Construct from a MeshedSurface
        UnsortedMeshedSurface(const MeshedSurface<Face>& surf);

        //- Move construct from a UnsortedMeshedSurface
        UnsortedMeshedSurface(UnsortedMeshedSurface<Face>&& surf);

        //- Move construct from MeshedSurface
        UnsortedMeshedSurface(MeshedSurface<Face>&& surf);

        //- Move construct from components (points, faces, zone ids).
        //  Zone information is fairly lightweight and is copied.
        UnsortedMeshedSurface
        (
            pointField&& pointLst,
            List<Face>&& faceLst,
            List<label>&& zoneIds,
            UList<surfZoneIdentifier>& tocInfo
        );

        //- Construct from file name (uses extension to determine type)
        UnsortedMeshedSurface(const fileName& name);

        //- Construct from file name (uses extension to determine type)
        UnsortedMeshedSurface(const fileName& name, const word& ext);

        //- Construct from Istream
        UnsortedMeshedSurface(Istream& is);

        //- Construct from objectRegistry and a named surface
        UnsortedMeshedSurface
        (
            const Time& t,
            const word& surfName = word::null
        );


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
            const fileName& name,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<UnsortedMeshedSurface> New(const fileName& name);


    //- Destructor
    virtual ~UnsortedMeshedSurface() = default;


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            UnsortedMeshedSurface,
            write,
            fileExtension,
            (
                const fileName& name,
                const UnsortedMeshedSurface<Face>& surf,
                const dictionary& options
            ),
            (name, surf, options)
        );

        //- Write to file, select based on its extension
        static void write
        (
            const fileName& name,
            const UnsortedMeshedSurface<Face>& surf,
            const dictionary& options = dictionary::null
        );

        //- Write to file, selected based on given extension
        static void write
        (
            const fileName& name,
            const word& ext,
            const UnsortedMeshedSurface<Face>& surf,
            const dictionary& options
        );


    // Member Functions

    // Access

        //- The surface size is the number of faces
        label size() const
        {
            return ParentType::size();
        }

        //- Reset size of face and zone list
        void setSize(const label);

        //- Return const access to the zone ids
        virtual const labelList& zoneIds() const
        {
            return zoneIds_;
        }

        //- Return const access to the zone table-of-contents
        const List<surfZoneIdentifier>& zoneToc() const
        {
            return zoneToc_;
        }

        //- Sort faces according to zoneIds
        //  Returns a surfZoneList and sets faceMap to index within faces()
        //  (i.e. map from original,unsorted to sorted)
        surfZoneList sortedZones(labelList& faceMap) const;

        //- Set zones to 0 and set a single zone
        void setOneZone();

        //- Set zone ids and zones
        void setZones(const surfZoneList& zoneLst);

        //- Set zone ids and zones
        void setZones(const labelUList& sizes, const UList<word>& names);

        //- Set zone ids and zones with default names
        void setZones(const labelUList& sizes);


    // Edit

        //- Clear all storage
        virtual void clear();

        //- Return new surface.
        //  \param[in] include the faces to select
        //  \param[out] pointMap from subsetMeshMap
        //  \param[out] faceMap from subsetMeshMap
        template<class BoolListType>
        UnsortedMeshedSurface subsetMesh
        (
            const BoolListType& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return new surface
        UnsortedMeshedSurface subsetMesh(const bitSet& include) const;

        //- Return new surface
        UnsortedMeshedSurface subsetMesh(const labelHashSet& include) const;


        //- Swap contents - disabled
        void swap(MeshedSurface<Face>& surf) = delete;

        //- Swap contents
        void swap(UnsortedMeshedSurface<Face>& surf);

        //- Transfer the contents of the argument and annul the argument
        void transfer(UnsortedMeshedSurface<Face>& surf);

        //- Transfer the contents of the argument and annul the argument
        void transfer(MeshedSurface<Face>& surf);

        //- Release (clear) stored zoneIds and return for reuse
        autoPtr<labelList> releaseZoneIds();


    // Read

        //- Read from file. Chooses reader based on explicit extension
        bool read(const fileName& name, const word& ext);

        //- Read from file. Chooses reader based on detected extension
        virtual bool read(const fileName& name);


    // Write

        //- Generic write routine. Chooses writer based on its extension.
        virtual void write
        (
            const fileName& name,
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, *this, options);
        }

        //- Write to database
        void write
        (
            const Time& t,
            const word& surfName = word::null
        ) const;


    // Member Operators

        //- Copy assignment
        void operator=(const UnsortedMeshedSurface<Face>& surf);

        //- Move assignment
        void operator=(UnsortedMeshedSurface<Face>&& surf);

        //- Conversion operator to MeshedSurfaceProxy
        operator MeshedSurfaceProxy<Face>() const;


    // IOstream Operators

        //- Read UnsortedMeshedSurface from Istream.
        //  Avoid using to read/write file content (fragile).
        friend Istream& operator>> <Face>
        (
            Istream& is,
            UnsortedMeshedSurface<Face>& surf
        );

        //- Write UnsortedMeshedSurface to Ostream.
        //  Avoid using to read/write file content (fragile).
        friend Ostream& operator<< <Face>
        (
            Ostream& os,
            const UnsortedMeshedSurface<Face>& surf
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UnsortedMeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
