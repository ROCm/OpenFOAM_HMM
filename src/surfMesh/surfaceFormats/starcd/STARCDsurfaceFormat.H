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
    Foam::fileFormats::STARCDsurfaceFormat

Description
    Read/write the surface shells from pro-STAR vrt/cel files.

Note
    Uses the extension @a .inp (input) to denote the format.

See Also
    Foam::meshReaders::STARCD

SourceFiles
    STARCDsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDsurfaceFormat_H
#define STARCDsurfaceFormat_H

#include "IFstream.H"
#include "Ostream.H"
#include "OFstream.H"
#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class STARCDsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class STARCDsurfaceFormat
:
    public UnsortedMeshedSurface<Face>
{
    //- Private typedefs for convenience
    typedef STARCDsurfaceFormat<Face>   ThisType;
    typedef UnsortedMeshedSurface<Face> ParentType;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        STARCDsurfaceFormat(const STARCDsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const STARCDsurfaceFormat<Face>&);


        static bool readHeader(IFstream&, const word&);

        static void writeHeader(Ostream&, const char* filetype);

        static void writePoints(Ostream&, const pointField&);

        static inline void writeShell
        (
            Ostream&,
            const Face&,
            const label cellId,
            const label cellTableId
        );

public:

    // Constructors

        //- Construct null
        STARCDsurfaceFormat();

        //- Construct from file name
        STARCDsurfaceFormat(const fileName&);

    // Selectors

        //- Read file and return keyedSurface
        static autoPtr<UnsortedMeshedSurface<Face> > New
        (
            const fileName& fName
        )
        {
            return autoPtr<UnsortedMeshedSurface<Face> >
            (
                new STARCDsurfaceFormat<Face>(fName)
            );
        }

    // Destructor
        virtual ~STARCDsurfaceFormat()
        {}


    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);

        //- Write UnsortedMeshedSurface
        static void write(const fileName&, const UnsortedMeshedSurface<Face>&);

        //- Write MeshedSurface
        static void write(const fileName&, const MeshedSurface<Face>&);

        //- Write to Ostream as one large file - not really useful
        virtual void write(Ostream&) const
        {
            notImplemented("fileFormats::STARCDsurfaceFormat::write(Ostream&)");
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "STARCDsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
