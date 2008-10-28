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
    Foam::fileFormats::AC3DfileFormat

Description
    Provide a means of reading/writing AC3D format.

    http://www.inivis.com/ac3d/man/ac3dfileformat.html

Note
   Since the faces are already organized as patches, the reader could be
   optimized for this, but at expense of losing a common reading approach.

SourceFiles
    AC3DfileFormat.C

\*---------------------------------------------------------------------------*/

#ifndef AC3DfileFormat_H
#define AC3DfileFormat_H

#include "Ostream.H"
#include "OFstream.H"
#include "keyedSurface.H"
#include "meshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

class AC3DfileFormat;

/*---------------------------------------------------------------------------*\
                      Class AC3DfileFormat Declaration
\*---------------------------------------------------------------------------*/

class AC3DfileFormat
:
    public keyedSurface
{

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
        static void writeHeader(Ostream&, const List<surfacePatch>&);

        //- Disallow default bitwise copy construct
        AC3DfileFormat(const AC3DfileFormat&);

        //- Disallow default bitwise assignment
        void operator=(const AC3DfileFormat&);

public:

    // Constructors

        //- Construct null
        AC3DfileFormat();

        //- Construct from file name
        AC3DfileFormat(const fileName&, const bool triangulate=false);

    // Selectors

        //- Read file and return keyedSurface
        static autoPtr<keyedSurface> New
        (
            const fileName& fName,
            const bool triangulate=false
        )
        {
            return autoPtr<keyedSurface>
            (
                new AC3DfileFormat(fName, triangulate)
            );
        }

    // Destructor

    // Member Functions

        // Write

        //- Write keyedSurface
        static void write(Ostream&, const keyedSurface&);

        //- Write keyedSurface
        static void write(const fileName& fName, const keyedSurface& surf)
        {
            write(OFstream(fName)(), surf);
        }

        //- Write meshedSurface
        static void write(Ostream&, const meshedSurface&);

        //- Write meshedSurface
        static void write(const fileName& fName, const meshedSurface& surf)
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
