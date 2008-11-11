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
    Foam::fileFormats::STLfileFormat

Description
    Provide a means of reading/writing STL files (ASCII and binary).

SourceFiles
    STLfileFormat.C
    STLfileFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLfileFormat_H
#define STLfileFormat_H

#include "STLtriangle.H"
#include "IFstream.H"
#include "Ostream.H"
#include "OFstream.H"
#include "keyedSurface.H"
#include "meshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                       Class STLfileFormat Declaration
\*---------------------------------------------------------------------------*/

class STLfileFormat
:
    public keyedSurface
{
    // Private data

        //- The number of bytes in the STL binary header
        static const int headerSize=80;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        STLfileFormat(const STLfileFormat&);

        //- Disallow default bitwise assignment
        void operator=(const STLfileFormat&);

        bool readASCII(IFstream&, const off_t);
        bool readBINARY(IFstream&, const off_t);

        static bool detectBINARY(const fileName&);

        //- write face in ASCII
        static inline void writeShell
        (
            Ostream&,
            const pointField&,
            const face&,
            const vector& norm
        );

        //- write face in BINARY
        static inline void writeShell
        (
            ostream&,
            const pointField&,
            const face&,
            const vector&,
            const label patchI
        );


        //- Write keySurface
        static void writeASCII(Ostream&, const keyedSurface&);

        //- Write keyedSurface
        static void writeBINARY(ostream&, const keyedSurface&);

        //- Write meshedSurface
        static void writeASCII(Ostream&, const meshedSurface&);

        //- Write meshedSurface
        static void writeBINARY(ostream&, const meshedSurface&);

public:

    // Constructors

        //- Construct null
        STLfileFormat();

        //- Construct from file name
        STLfileFormat(const fileName&, const bool triangulate=true);

    // Selectors

        //- Read file and return keyedSurface
        static autoPtr<keyedSurface> New
        (
            const fileName& fName,
            const bool triangulate=true
        )
        {
            return autoPtr<keyedSurface>
            (
                new STLfileFormat(fName,true)
            );
        }

    // Destructor

    // Member Functions

        //- Read from file
        virtual bool read(const fileName&, const bool triangulate);

        //- Write keyedSurface
        //  The ASCII output is sorted by patch, the binary output is unsorted
        static void write(Ostream&, const keyedSurface&);

        //- Write keyedSurface
        //  The ASCII output is sorted by patch, the binary output is unsorted
        static void write(const fileName&, const keyedSurface&);

        //- Write meshedSurface
        static void write(Ostream&, const meshedSurface&);

        //- Write meshedSurface
        static void write(const fileName&, const meshedSurface&);

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
