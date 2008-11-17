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
    Foam::fileFormats::GTSsurfaceFormat

Description
    Provide a means of reading/writing GTS format.

SourceFiles
    GTSsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef GTSsurfaceFormat_H
#define GTSsurfaceFormat_H

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
                     Class GTSsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class GTSsurfaceFormat
:
    public UnsortedMeshedSurface<Face>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        GTSsurfaceFormat(const GTSsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const GTSsurfaceFormat<Face>&);

public:

    // Constructors

        //- Construct from file name
        GTSsurfaceFormat(const fileName&);

    // Selectors

        //- Read file and return surface
        static autoPtr<UnsortedMeshedSurface<Face> > New
        (
            const fileName& fName
        )
        {
            return autoPtr<UnsortedMeshedSurface<Face> >
            (
                new GTSsurfaceFormat<Face>(fName)
            );
        }

    // Destructor
        virtual ~GTSsurfaceFormat()
        {}

    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);

        //- Write UnsortedMeshedSurface
        //  The output is sorted by regions
        static void write
        (
            Ostream&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface
        //  The output is sorted by regions
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
#   include "GTSsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
