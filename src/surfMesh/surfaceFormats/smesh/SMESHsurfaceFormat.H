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
    Foam::fileFormats::SMESHsurfaceFormat

Description
    Provide a means of writing tetgen SMESH format.

    Tetgen http://tetgen.berlios.de

See Also
    File format information:
    http://tetgen.berlios.de/fformats.smesh.html

SourceFiles
    SMESHsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef SMESHsurfaceFormat_H
#define SMESHsurfaceFormat_H

#include "IFstream.H"
#include "OFstream.H"
#include "Ostream.H"
#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "SMESHsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class SMESHsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class SMESHsurfaceFormat
:
    public MeshedSurface<Face>,
    public SMESHsurfaceFormatCore
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        SMESHsurfaceFormat(const SMESHsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const SMESHsurfaceFormat<Face>&);

public:

    // Constructors

        //- Construct null
        SMESHsurfaceFormat();

    // Selectors

    // Destructor

        virtual ~SMESHsurfaceFormat()
        {}


    // Member Functions

        //- Write MeshedSurface
        static void write
        (
            Ostream&,
            const MeshedSurface<Face>&
        );

        //- Write MeshedSurface
        static void write
        (
            const fileName& name,
            const MeshedSurface<Face>& surf
        )
        {
            write(OFstream(name)(), surf);
        }

        //- Write UnsortedMeshedSurface
        //  The output is sorted by region.
        static void write
        (
            Ostream&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface
        //  The output is sorted by region.
        static void write
        (
            const fileName& name,
            const UnsortedMeshedSurface<Face>& surf
        )
        {
            write(OFstream(name)(), surf);
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
#   include "SMESHsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
