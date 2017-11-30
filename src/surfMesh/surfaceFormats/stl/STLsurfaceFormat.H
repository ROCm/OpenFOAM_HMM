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
    Foam::fileFormats::STLsurfaceFormat

Description
    Provide a means of reading/writing STL files (ASCII and BINARY).

Note
    For efficiency, the zones are sorted before creating the faces.
    The class is thus derived from MeshedSurface.

SourceFiles
    STLsurfaceFormat.C
    STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLsurfaceFormat_H
#define STLsurfaceFormat_H

#include "STLReader.H"
#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::STLsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class STLsurfaceFormat
:
    public MeshedSurface<Face>,
    public STLCore
{
    // Private Member Functions

        //- Write Face (ASCII)
        static inline void writeShell
        (
            Ostream& os,
            const UList<point>& pts,
            const Face& f
        );

        //- Write Face (BINARY)
        static inline void writeShell
        (
            ostream& os,
            const UList<point>& pts,
            const Face& f,
            const label zoneI
        );


        //- Disallow default bitwise copy construct
        STLsurfaceFormat(const STLsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const STLsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct from file name
        STLsurfaceFormat(const fileName& filename);


    //- Destructor
    virtual ~STLsurfaceFormat() = default;


    // Static Member Functions

        //- Write surface mesh components by proxy (as ASCII)
        static void writeAscii
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf
        );

        //- Write surface mesh components by proxy (as BINARY)
        static void writeBinary
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf
        );

        //- Write UnsortedMeshedSurface (as ASCII) sorted by zone
        static void writeAscii
        (
            const fileName& filename,
            const UnsortedMeshedSurface<Face>& surf
        );

        //- Write UnsortedMeshedSurface (as BINARY) unsorted by zone
        static void writeBinary
        (
            const fileName& filename,
            const UnsortedMeshedSurface<Face>& surf
        );

        //- Write surface mesh components by proxy
        //  as ASCII or BINARY or dependent on the extension
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            const STLFormat format
        );

        //- Write UnsortedMeshedSurface
        //  as ASCII or BINARY or dependent on the extension
        static void write
        (
            const fileName& filename,
            const UnsortedMeshedSurface<Face>& surf,
            const STLFormat format
        );


        //- Write surface mesh components by proxy
        //  as ASCII or BINARY, depending on the extension
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            const dictionary& options = dictionary::null
        );

        //- Write UnsortedMeshedSurface
        //  as ASCII or BINARY, depending on the extension
        static void write
        (
            const fileName& filename,
            const UnsortedMeshedSurface<Face>& surf,
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Read from file
        virtual bool read(const fileName& filename);

        //- Write object
        virtual void write
        (
            const fileName& name,
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this), options);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "STLsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
