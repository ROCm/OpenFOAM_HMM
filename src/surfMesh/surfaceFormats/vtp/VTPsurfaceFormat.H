/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::fileFormats::VTPsurfaceFormat

Description
    Write surfaces in VTP (xml) format.

    The output is never sorted by zone.

SourceFiles
    VTPsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTPsurfaceFormat_H
#define VTPsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "VTPsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::VTPsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class VTPsurfaceFormat
:
    public MeshedSurface<Face>,
    public VTPsurfaceFormatCore
{
    // Private Member Functions

        //- Write polygons
        static void writePolys
        (
            vtk::formatter& format,
            const UList<Face>& faces
        );


        //- Disallow default bitwise copy construct
        VTPsurfaceFormat(const VTPsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const VTPsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct null
        VTPsurfaceFormat() = default;


    //- Destructor
    virtual ~VTPsurfaceFormat() = default;


    // Static Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            const dictionary& options = dictionary::null
        );

        //- Write UnsortedMeshedSurface, the output remains unsorted
        static void write
        (
            const fileName& filename,
            const UnsortedMeshedSurface<Face>& surf,
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Write object file
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
    #include "VTPsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
