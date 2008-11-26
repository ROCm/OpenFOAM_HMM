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
    Foam::fileFormats::STLsurfaceFormat

Description
    Provide a means of reading/writing STL files (ASCII and binary).

Note
    For efficiency, the regions are sorted before creating the faces.
    The class is thus derived from MeshedSurface.

SourceFiles
    STLsurfaceFormat.C
    STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLsurfaceFormat_H
#define STLsurfaceFormat_H

#include "STLsurfaceFormatCore.H"
#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class STLsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class STLsurfaceFormat
:
    public MeshedSurface<Face>
{
    // Private Member Functions

        //- Write Face in ASCII
        static inline void writeShell
        (
            Ostream&,
            const pointField&,
            const Face&,
            const vector& norm
        );

        //- Write Face in BINARY
        static inline void writeShell
        (
            ostream&,
            const pointField&,
            const Face&,
            const vector&,
            const label patchI
        );


        //- Write UnsortedMeshedSurface
        static void writeASCII(Ostream&, const UnsortedMeshedSurface<Face>&);

        //- Write UnsortedMeshedSurface
        static void writeBINARY(ostream&, const UnsortedMeshedSurface<Face>&);

        //- Write MeshedSurface
        static void writeASCII(Ostream&, const MeshedSurface<Face>&);

        //- Write MeshedSurface
        static void writeBINARY(ostream&, const MeshedSurface<Face>&);


        //- Disallow default bitwise copy construct
        STLsurfaceFormat(const STLsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const STLsurfaceFormat<Face>&);
public:

    // Constructors

        //- Construct from file name
        STLsurfaceFormat(const fileName&);

    // Selectors

        //- Read file and return surface
        static autoPtr<MeshedSurface<Face> > New
        (
            const fileName& fName
        )
        {
            return autoPtr<MeshedSurface<Face> >
            (
                new STLsurfaceFormat<Face>(fName)
            );
        }

    // Destructor
        virtual ~STLsurfaceFormat()
        {}

    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);

        //- Write MeshedSurface (as ASCII)
        static void write
        (
            Ostream&,
            const MeshedSurface<Face>&
        );

        //- Write MeshedSurface
        static void write
        (
            const fileName&,
            const MeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface (as ASCII) sorted by region
        static void write
        (
            Ostream&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface
        //  ASCII output is sorted by region; binary output is unsorted
        static void write
        (
            const fileName&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write object (only ASCII output)
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
#   include "STLsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
