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
    Foam::MeshedSurface

Description
    A surface geometry mesh with zone information, not to be confused with
    the similarly named surfaceMesh, which actually refers to the cell faces
    of a volume mesh.

    A MeshedSurface can have zero or more surface zones (roughly equivalent
    to faceZones for a polyMesh). If surface zones are defined, they must
    be contiguous and cover all of the faces.

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

#include "PrimitivePatch.H"
#include "PatchTools.H"
#include "pointField.H"
#include "face.H"
#include "labelledTri.H"
#include "bitSet.H"
#include "HashSet.H"
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class Time;
class surfMesh;
class polyBoundaryMesh;
class Istream;
class Ostream;

template<class Face> class MeshedSurface;
template<class Face> class MeshedSurfaceProxy;
template<class Face> class UnsortedMeshedSurface;

template<class Face>
Istream& operator>>(Istream&, MeshedSurface<Face>&);
template<class Face>
Ostream& operator<<(Ostream&, const MeshedSurface<Face>&);

/*---------------------------------------------------------------------------*\
                      Class MeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class MeshedSurface
:
    public PrimitivePatch<Face, ::Foam::List, pointField, point>,
    public fileFormats::surfaceFormatsCore
{
    // Friends, regardless of face representations
    template<class Face2> friend class MeshedSurface;
    template<class Face2> friend class UnsortedMeshedSurface;

    // Friendship with surfMesh is needed for transferring
    friend class surfMesh;


private:

    // Private typedefs for convenience

        typedef PrimitivePatch
        <
            Face,
            ::Foam::List,
            pointField,
            point
        >
        ParentType;

        typedef UnsortedMeshedSurface<Face>  FriendType;
        typedef MeshedSurfaceProxy<Face>     ProxyType;


    // Private Member Data

        //- Zone information
        // (face ordering nFaces/startFace only used during reading/writing)
        surfZoneList zones_;


    // Private Member Functions

        //- Read/construct from Istream
        Istream& read(Istream& is);

        //- Write to Ostream
        Ostream& write(Ostream& os) const;


protected:

    // Protected Member functions

        //- Transfer points/zones from 'face' to other other shapes.
        //  Eg, transcribe face to triFace, or face -> labelledTri, including
        //  any addZonesToFaces adjustment.
        //  No general form, only specializations.
        void transcribe(MeshedSurface<face>& surf);

        //- Basic sanity check on zones
        void checkZones();

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

        //- Non-const access to the zones
        surfZoneList& storedZones()
        {
            return zones_;
        }

        //- Sort faces by zones and store sorted faces
        void sortFacesAndStore
        (
            DynamicList<Face>& unsortedFaces,
            DynamicList<label>& zoneIds,
            const bool sorted
        );

        //- Set new zones from faceMap
        virtual void remapFaces(const labelUList& faceMap);


public:

    // Public typedefs

        //- Face type used
        typedef Face FaceType;

        //- Runtime type information
        ClassName("MeshedSurface");


    // Static

        //- Can we read this file format?
        static bool canRead(const fileName& name, bool verbose=false);

        //- Can we read this file format?
        static bool canReadType(const word& ext, bool verbose=false);

        //- Can we write this file format?
        static bool canWriteType(const word& ext, bool verbose=false);

        //- Known readable file-types
        static wordHashSet readTypes();

        //- Known writable file-types
        static wordHashSet writeTypes();


    // Constructors

        //- Construct null, an empty surface
        MeshedSurface();

        //- Copy construct
        MeshedSurface(const MeshedSurface& surf);

        //- Copy construct from an UnsortedMeshedSurface
        MeshedSurface(const UnsortedMeshedSurface<Face>& surf);

        //- Move construct
        MeshedSurface(MeshedSurface&& surf);

        //- Move construct from an UnsortedMeshedSurface
        MeshedSurface(UnsortedMeshedSurface<Face>&& surf);

        //- Copy construct from components (points, faces, zones).
        MeshedSurface
        (
            const pointField& pointLst,
            const UList<Face>& faceLst,
            const UList<surfZone>& zoneLst
        );

        //- Move construct from components (points, faces).
        //  Zone information is fairly lightweight and is copied.
        MeshedSurface
        (
            pointField&& pointLst,
            List<Face>&& faceLst,
            const UList<surfZone>& zoneLst
        );

        //- Copy construct from components (points, faces).
        //  Use zone information if available
        MeshedSurface
        (
            const pointField& pointLst,
            const UList<Face>& faceLst,
            const labelUList& zoneSizes = labelUList(),
            const UList<word>& zoneNames = UList<word>()
        );

        //- Move construct from components (points, faces).
        //  Use zone information if available
        MeshedSurface
        (
            pointField&& pointLst,
            List<Face>&& faceLst,
            const labelUList& zoneSizes = labelUList(),
            const UList<word>& zoneNames = UList<word>()
        );

        //- Construct from a boundary mesh with local points/faces
        MeshedSurface
        (
            const polyBoundaryMesh& bMesh,
            const bool globalPoints = false
        );

        //- Construct from a surfMesh
        MeshedSurface(const surfMesh& mesh);

        //- Construct from file name (uses extension to determine type)
        MeshedSurface(const fileName& name);

        //- Construct from file name (uses extension to determine type)
        MeshedSurface(const fileName& name, const word& ext);

        //- Construct from Istream
        MeshedSurface(Istream& is);

        //- Construct from database
        MeshedSurface(const Time& runTime, const word& surfName = word::null);


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
            const fileName& name,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<MeshedSurface> New(const fileName& name);


    //- Destructor
    virtual ~MeshedSurface();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            UnsortedMeshedSurface,
            write,
            fileExtension,
            (
                const fileName& name,
                const MeshedSurface<Face>& surf,
                const dictionary& options
            ),
            (name, surf, options)
        );

        //- Write to file, selecting writer based on its extension
        static void write
        (
            const fileName& name,
            const MeshedSurface<Face>& surf,
            const dictionary& options = dictionary::null
        );

        //- Write to file, selecting writer based on the given extension
        static void write
        (
            const fileName& name,
            const word& ext,
            const MeshedSurface<Face>& surf,
            const dictionary& options = dictionary::null
        );


    // Member Functions

    // Access

        //- The surface size is the number of faces
        label size() const
        {
            return ParentType::size();
        }

        //- Return const access to the faces
        inline const List<Face>& surfFaces() const
        {
            return static_cast<const List<Face>&>(*this);
        }

        //- Const access to the surface zones.
        //  If zones are defined, they must be contiguous and cover the
        //  entire surface
        const surfZoneList& surfZones() const
        {
            return zones_;
        }

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


    // Edit

        //- Clear all storage
        virtual void clear();


        //- Add surface zones
        virtual void addZones
        (
            const UList<surfZone>&,
            const bool cullEmpty=false
        );

        //- Add surface zones
        virtual void addZones
        (
            const labelUList& sizes,
            const UList<word>& names,
            const bool cullEmpty=false
        );

        //- Add surface zones
        virtual void addZones
        (
            const labelUList& sizes,
            const bool cullEmpty=false
        );

        //- Propagate zone information on face regions.
        //  Normally a no-op, only used by the labelledTri specialization.
        //  Specializations return true, others return false.
        bool addZonesToFaces();


        //- Remove surface zones
        virtual void removeZones();


        //- Move points
        virtual void movePoints(const pointField& newPoints);

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar scaleFactor);

        //- Remove invalid faces
        virtual void cleanup(const bool verbose);

        virtual bool stitchFaces
        (
            const scalar tol=SMALL,
            const bool verbose=false
        );

        virtual bool checkFaces
        (
            const bool verbose=false
        );

        //- Count number of triangles.
        virtual label nTriangles() const;

        //- Count number of triangles, returning a face map of original ids.
        //  The faceMap is zero-sized when no triangulation would be needed.
        virtual label nTriangles(List<label>& faceMap) const;

        //- Triangulate in-place, returning the number of triangles added.
        virtual label triangulate();

        //- Triangulate in-place, returning the number of triangles added
        //  and setting a map of original face Ids.
        //  The faceMap is zero-sized when no triangulation was done.
        virtual label triangulate(List<label>& faceMap);

        //- Return new surface.
        //
        //  \param[in] include the faces to select
        //  \param[out] pointMap from subsetMeshMap
        //  \param[out] faceMap from subsetMeshMap
        template<class BoolListType>
        MeshedSurface subsetMesh
        (
            const BoolListType& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return new surface
        MeshedSurface subsetMesh(const bitSet& include) const;

        //- Return new surface
        MeshedSurface subsetMesh(const labelHashSet& include) const;

        //- Swap contents
        void swap(MeshedSurface<Face>& surf);

        //- Transfer the components
        void transfer(pointField& pointLst, List<Face>& faceLst);

        //- Transfer the contents of the argument and annul the argument
        void transfer(MeshedSurface<Face>& surf);

        //- Transfer the contents of the argument and annul the argument
        void transfer(UnsortedMeshedSurface<Face>& surf);

        //- Release (clear) geometry and return for reuse
        autoPtr<MeshedSurface<Face>> releaseGeom();

        //- Swap the stored faces
        void swapFaces(List<Face>& faces);

        //- Swap the stored points
        void swapPoints(pointField& points);

        //- Swap the stored zones
        void swapZones(surfZoneList& zones);


    // Read

        //- Read from file. Chooses reader based on explicit extension
        bool read(const fileName& name, const word& ext);

        //- Read from file. Chooses reader based on detected extension
        virtual bool read(const fileName& name);


    // Write

        void writeStats(Ostream& os) const;

        //- Generic write routine. Chooses writer based on extension.
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


    // Member operators

        //- Copy assignment
        void operator=(const MeshedSurface<Face>& surf);

        //- Move assignment
        void operator=(MeshedSurface<Face>&& surf);

        //- Conversion operator to MeshedSurfaceProxy
        operator MeshedSurfaceProxy<Face>() const;


    // IOstream Operators

        //- Read MeshedSurface from Istream.
        //  Avoid using to read/write file content (fragile).
        friend Istream& operator>> <Face>
        (
            Istream& is,
            MeshedSurface<Face>& surf
        );


        //- Write MeshedSurface to Ostream.
        //  Avoid using to read/write file content (fragile).
        friend Ostream& operator<< <Face>
        (
            Ostream& os,
            const MeshedSurface<Face>& surf
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialization for labelledTri.
template<>
bool MeshedSurface<labelledTri>::addZonesToFaces();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
