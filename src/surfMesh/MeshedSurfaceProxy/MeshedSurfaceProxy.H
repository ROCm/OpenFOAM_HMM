/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Face> class MeshedSurface;

/*---------------------------------------------------------------------------*\
                     Class MeshedSurfaceProxy Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class MeshedSurfaceProxy
:
    public fileFormats::surfaceFormatsCore
{
    // Private Data

        const pointField& points_;

        const UList<Face>& faces_;

        const UList<surfZone>& zones_;

        const UList<label>& faceMap_;

public:

    // Public Typedefs

        //- The face type
        typedef Face face_type;

        //- The point type
        typedef point point_type;


    //- Declare type-name (with debug switch)
    ClassName("MeshedSurfaceProxy");


    // Static Functions

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
            const labelUList& faceMap = labelUList::null()
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
                IOstreamOption streamOpt,
                const dictionary& options
            ),
            (name, surf, streamOpt, options)
        );

        //- Write to file, select based on its extension
        static void write
        (
            const fileName& name,
            const MeshedSurfaceProxy& surf,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );

        //- Write to file, selected based on given extension
        static void write
        (
            const fileName& name,
            const word& ext,
            const MeshedSurfaceProxy& surf,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );


    // Member Functions

        // Access

            //- The surface size is the number of faces
            label size() const
            {
                return faces_.size();
            }

            //- Return const access to the points
            const pointField& points() const
            {
                return points_;
            }

            //- Return const access to the faces
            const UList<Face>& surfFaces() const
            {
                return faces_;
            }

            //- Const access to the surface zones.
            //  If zones are defined, they must be contiguous and cover the
            //  entire surface
            const UList<surfZone>& surfZones() const
            {
                return zones_;
            }

            //- Const access to the faceMap, zero-sized when unused
            const labelUList& faceMap() const
            {
                return faceMap_;
            }

            //- Can/should use faceMap?
            bool useFaceMap() const
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
                IOstreamOption streamOpt = IOstreamOption(),
                const dictionary& options = dictionary::null
            ) const
            {
                write(name, *this, streamOpt, options);
            }

            //- Generic write routine. Chooses writer based on extension.
            virtual void write
            (
                const fileName& name,
                const word& ext,
                IOstreamOption streamOpt = IOstreamOption(),
                const dictionary& options = dictionary::null
            ) const
            {
                write(name, ext, *this, streamOpt, options);
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
