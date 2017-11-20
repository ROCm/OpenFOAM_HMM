/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::MeshedSurfaceProxy

Description
    A proxy for writing MeshedSurface, UnsortedMeshedSurface and surfMesh
    to various file formats.

SourceFiles
    MeshedSurfaceProxy.C
    MeshedSurfaceProxys.C

\*---------------------------------------------------------------------------*/

#ifndef MeshedSurfaceProxy_H
#define MeshedSurfaceProxy_H

#include "pointField.H"
#include "labelledTri.H"
#include "HashSet.H"
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Face> class MeshedSurface;

/*---------------------------------------------------------------------------*\
                     Class MeshedSurfaceProxy Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class MeshedSurfaceProxy
:
    public fileFormats::surfaceFormatsCore
{
    // Private Member Data

        const pointField& points_;

        const UList<Face>& faces_;

        const UList<surfZone>& zones_;

        const UList<label>& faceMap_;

public:

    // Public typedefs

        //- Face type used
        typedef Face FaceType;


    // Static

        //- Runtime type information
        ClassName("MeshedSurfaceProxy");

        //- The file format types that can be written via MeshedSurfaceProxy
        static wordHashSet writeTypes();

        //- Can this file format type be written via MeshedSurfaceProxy?
        static bool canWriteType(const word& ext, const bool verbose=false);


    // Constructors

        //- Construct from component references
        MeshedSurfaceProxy
        (
            const pointField& pointLst,
            const UList<Face>& faceLst,
            const UList<surfZone>& zoneLst = List<surfZone>(),
            const labelUList& faceMap = List<label>()
        );


    //- Destructor
    virtual ~MeshedSurfaceProxy() = default;


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            MeshedSurfaceProxy,
            write,
            fileExtension,
            (
                const fileName& name,
                const MeshedSurfaceProxy<Face>& surf,
                const dictionary& options
            ),
            (name, surf, options)
        );

        //- Write to file, select based on its extension
        static void write
        (
            const fileName& name,
            const MeshedSurfaceProxy& surf,
            const dictionary& options
        );

        //- Write to file, selected based on given extension
        static void write
        (
            const fileName& name,
            const word& ext,
            const MeshedSurfaceProxy& surf,
            const dictionary& options
        );


    // Member Functions

        // Access

            //- The surface size is the number of faces
            inline label size() const
            {
                return faces_.size();
            }

            //- Return const access to the points
            inline const pointField& points() const
            {
                return points_;
            }

            //- Return const access to the faces
            inline const UList<Face>& surfFaces() const
            {
                return faces_;
            }

            //- Const access to the surface zones.
            //  If zones are defined, they must be contiguous and cover the
            //  entire surface
            inline const UList<surfZone>& surfZones() const
            {
                return zones_;
            }

            //- Const access to the faceMap, zero-sized when unused
            inline const labelUList& faceMap() const
            {
                return faceMap_;
            }

            //- Use faceMap?
            inline bool useFaceMap() const
            {
                return faceMap_.size() == faces_.size();
            }

            //- Count number of triangles.
            inline label nTriangles() const;


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

            //- Generic write routine. Chooses writer based on extension.
            virtual void write
            (
                const fileName& name,
                const word& ext,
                const dictionary& options = dictionary::null
            ) const
            {
                write(name, ext, *this, options);
            }

            //- Write to database
            virtual void write
            (
                const Time& t,
                const word& surfName = word::null
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeshedSurfaceProxy.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
