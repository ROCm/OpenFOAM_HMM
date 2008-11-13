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
    Foam::fileFormats::AC3DsurfaceFormat

Description
    Provide a means of reading/writing AC3D format.

    http://www.inivis.com/ac3d/man/ac3dfileformat.html

Note
   Since the faces are already organized as patches, the reader could be
   optimized for this, but at expense of losing a common reading approach.

SourceFiles
    AC3DfileFormat.C

\*---------------------------------------------------------------------------*/

#ifndef AC3DsurfaceFormat_H
#define AC3DsurfaceFormat_H

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
                      Class AC3DfileFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class AC3DsurfaceFormat
:
    public UnsortedMeshedSurface<Face>
{
    //- Private typedefs for convenience
    typedef AC3DsurfaceFormat<Face>     ThisType;
    typedef UnsortedMeshedSurface<Face> ParentType;

    // Private Member Functions

        //- Read a type via IStringStream
        template<class Type>
        static Type parse(const string&);

        //- Read cmd, args from IFstream
        static bool readCmd(IFstream&, string& cmd, string& args);

        //- Cue up to cmd, reading args
        static bool cueTo(IFstream&, const string& cmd, string& args);

        //- Cue up to cmd, reading args or exit with a FatalError
        //  returns the command args
        static string cueToOrDie
        (
            IFstream&,
            const string& cmd,
            const string& errorMsg=string::null
        );

        //- Write header with materials
        static void writeHeader(Ostream&, const List<surfGroup>&);

        //- Disallow default bitwise copy construct
        AC3DsurfaceFormat(const AC3DsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const AC3DsurfaceFormat<Face>&);

public:

    // Constructors

        //- Construct null
        AC3DsurfaceFormat();

        //- Construct from file name
        AC3DsurfaceFormat(const fileName&);

    // Selectors

        //- Read file and return surface
        static autoPtr<UnsortedMeshedSurface<Face> > New
        (
            const fileName& fName
        )
        {
            return autoPtr<UnsortedMeshedSurface<Face> >
            (
                new AC3DsurfaceFormat<Face>(fName)
            );
        }

    // Destructor

        virtual ~AC3DsurfaceFormat()
        {}


    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);

        //- Write UnsortedMeshedSurface
        static void write
        (
            Ostream&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface
        static void write
        (
            const fileName& fName,
            const UnsortedMeshedSurface<Face>& surf
        )
        {
            write(OFstream(fName)(), surf);
        }

        //- Write MeshedSurface
        static void write
        (
            Ostream&,
            const MeshedSurface<Face>&
        );

        //- Write MeshedSurface
        static void write
        (
            const fileName& fName,
            const MeshedSurface<Face>& surf
        )
        {
            write(OFstream(fName)(), surf);
        }

        //- Write object
        virtual void write(Ostream& os) const
        {
            write(os, *this);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "AC3DsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
