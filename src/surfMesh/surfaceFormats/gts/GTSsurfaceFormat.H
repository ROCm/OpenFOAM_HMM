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
    Foam::fileFormats::GTSsurfaceFormat

Description
    Read/write GTS format.

    The output is never sorted by zone and is only written if it consists
    entirely of triangles.

SourceFiles
    GTSsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef GTSsurfaceFormat_H
#define GTSsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::GTSsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class GTSsurfaceFormat
:
    public UnsortedMeshedSurface<Face>
{
    // Private Member Functions

        //- Check and raise FatalError if output is not triangulated
        //  Triangulatinng on-the-fly is otherwise too annoying
        static bool checkIfTriangulated(const UList<Face>& faceLst);

        //- Disallow default bitwise copy construct
        GTSsurfaceFormat(const GTSsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const GTSsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct from file name
        GTSsurfaceFormat(const fileName& filename);


    //- Destructor
    virtual ~GTSsurfaceFormat() = default;


    // Static Member Functions

        //- Write MeshedSurface
        static void write
        (
            const fileName& filename,
            const MeshedSurface<Face>& surf,
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

        //- Read from file
        virtual bool read(const fileName& filename);

        //- Write object
        virtual void write
        (
            const fileName& name,
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, *this, options);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GTSsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
