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
    Foam::fileFormats::VTKsurfaceFormat

Description
    Provide a means of writing VTK legacy format.
    The output is never sorted by zone.

SourceFiles
    VTKsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTKsurfaceFormat_H
#define VTKsurfaceFormat_H

#include "Ostream.H"
#include "OFstream.H"
#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "VTKsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class VTKsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class VTKsurfaceFormat
:
    public MeshedSurface<Face>,
    public VTKsurfaceFormatCore
{
    // Private Member Functions

        //- Write header information about number of polygon points
        static void writeHeaderPolygons(Ostream&, const UList<Face>&);

        //- Disallow default bitwise copy construct
        VTKsurfaceFormat(const VTKsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const VTKsurfaceFormat<Face>&);

public:

    // Constructors

        //- Construct null
        VTKsurfaceFormat();

    // Selectors

    // Destructor

        virtual ~VTKsurfaceFormat()
        {}

    // Member Functions

        // Write

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
        //  The output remains unsorted
        static void write
        (
            Ostream&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface
        //  The output remains unsorted
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
#   include "VTKsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
