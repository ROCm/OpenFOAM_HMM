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
    Foam::surfMesh

Description
    A surface mesh consisting of general polygon faces.

SourceFiles
    surfMesh.C
    surfMeshClear.C
    surfMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef surfMesh_H
#define surfMesh_H

#include "surfaceRegistry.H"
#include "PrimitivePatch.H"
#include "pointField.H"
#include "faceList.H"
#include "pointIOField.H"
#include "faceIOList.H"
#include "labelIOList.H"
#include "surfZoneIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Face>
class MeshedSurface;

/*---------------------------------------------------------------------------*\
                      Class surfMeshAllocator Declaration
\*---------------------------------------------------------------------------*/

//- A helper class for storing points and faces
class surfMeshAllocator
{
protected:
    pointIOField storedPoints_;
    faceIOList   storedFaces_;

    surfMeshAllocator(const IOobject& ioPoints, const IOobject& ioFaces)
    :
        storedPoints_(ioPoints),
        storedFaces_(ioFaces)
    {}

    surfMeshAllocator
    (
        const IOobject& ioPoints,
        const pointField& points,
        const IOobject& ioFaces,
        const faceList& faces
    )
    :
        storedPoints_(ioPoints, points),
        storedFaces_(ioFaces, faces)
    {}

    surfMeshAllocator
    (
        const IOobject& ioPoints,
        const Xfer<pointField>& points,
        const IOobject& ioFaces,
        const Xfer<faceList>& faces
    )
    :
        storedPoints_(ioPoints, points),
        storedFaces_(ioFaces, faces)
    {}

};



/*---------------------------------------------------------------------------*\
                          Class surfMesh Declaration
\*---------------------------------------------------------------------------*/

class surfMesh
:
    public surfaceRegistry,
    public surfMeshAllocator,
    public PrimitivePatch<face, ::Foam::UList, const pointField&, point>
{

public:

    // Public data types

        //- Enumeration defining the state of the mesh after a read update.
        //  Used for post-processing applications, where the mesh
        //  needs to update based on the files written in time
        //  directores
        enum readUpdateState
        {
            UNCHANGED,
            POINTS_MOVED,
            TOPO_CHANGE,
            TOPO_PATCH_CHANGE
        };


private:

    // Private typedefs

        typedef PrimitivePatch
        <
            face,
            ::Foam::UList,
            const pointField&,
            point
        >
        MeshReference;

    // Permanent data

        // Zoning information

           //- Face zones
           surfZoneIOList surfZones_;


    // Private member functions

        //- Set a single zone
        void oneZone();

        //- Disallow construct as copy
        surfMesh(const surfMesh&);

        //- Disallow default bitwise assignment
        void operator=(const surfMesh&);

protected:

    // Protected Member Functions

        //- Non-const access to global points
        pointIOField& storedPoints()
        {
            return surfMeshAllocator::storedPoints_;
        }

        //- Non-const access to the faces
        faceIOList& storedFaces()
        {
            return surfMeshAllocator::storedFaces_;
        }

        //- Non-const access to the zones
        surfZoneIOList& storedZones()
        {
            return surfZones_;
        }

public:

    // Public typedefs

    typedef surfMesh Mesh;

    //- Placeholder only, but do not remove - it is needed for GeoMesh
    typedef bool BoundaryMesh;


    //- Runtime type information
    TypeName("surfMesh");

    //- Return the mesh sub-directory name (normally "surfMesh")
    static word meshSubDir;

    // Constructors

        //- Construct from IOobject, with alternative surface name
        explicit surfMesh(const IOobject&, const word& surfName="");

        //- Construct by transferring components (points, faces) without zones.
        //  surfZones are added using addZones() member function
        surfMesh
        (
            const IOobject&,
            const Xfer<pointField>&,
            const Xfer<faceList>&,
            const word& surfName=""
        );

        //- Construct copy/move from MeshedSurface
        surfMesh
        (
            const IOobject&,
            const Xfer< MeshedSurface<face> >& surf,
            const word& surfName=""
        );

    // Destructor

        virtual ~surfMesh();


    // Member Functions

        // Database

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshDir() const;

            //- Return the current instance directory for points
            //  Used in the consruction of geometric mesh data dependent
            //  on points
            const fileName& pointsInstance() const;

            //- Return the current instance directory for faces
            const fileName& facesInstance() const;

            //- Set the instance for mesh files
            void setInstance(const fileName&);


        // Access

            //- Return raw points
            virtual label nPoints() const;

            //- Return raw faces
            virtual label nFaces() const;

            //- Return number of faces
            virtual label size() const
            {
                return nFaces();
            }


            //- Return raw points
            virtual const pointField& points() const;

            //- Return raw faces
            virtual const faceList& faces() const;

            //- Return surface zones
            const surfZoneList& surfZones() const
            {
                return surfZones_;
            }

            //- Return non-const access to the zones
            surfZoneList& surfZones()
            {
                return surfZones_;
            }

            //- Check the surface zone definitions
            void checkZones();

            //- Add surface zones patches
            void addZones
            (
                const List<surfZone>&,
                const bool validate = true
            );

            //- Update the mesh based on the mesh files saved in
            //  time directories
            virtual readUpdateState readUpdate();

            //- Update the mesh corresponding to given map

            //- Remove surface zones
            void removeZones();

            //- Rename surface
            virtual void rename(const word&);

            //- Reset mesh primitive data.
            void resetPrimitives
            (
                const Xfer<pointField>& points,
                const Xfer<faceList>& faces,
                const Xfer<surfZoneList>& zones,
                const bool validate = true
            );


            //- Transfer the contents of the argument and annull the argument
            void transfer(MeshedSurface<face>&);

        //  Storage management

            //- Clear geometry
            void clearGeom();

            //- Clear addressing
            void clearAddressing();

            //- Clear all geometry and addressing unnecessary for CFD
            void clearOut();

            //- Clear primitive data (points, faces and cells)
            void clearPrimitives();

            //- Remove all files from mesh instance
            void removeFiles(const fileName& instanceDir) const;

            //- Remove all files from mesh instance()
            void removeFiles() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
