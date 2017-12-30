/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::fileFormats::STARCDsurfaceFormat

Description
    Read/write the surface shells from PROSTAR vrt/cel files.

Note
    Uses the extension \a .inp (input) to denote the format.

See also
    Foam::meshReaders::STARCD

SourceFiles
    STARCDsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDsurfaceFormat_H
#define STARCDsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "STARCDsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
              Class fileFormats::STARCDsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class STARCDsurfaceFormat
:
    public MeshedSurface<Face>,
    public STARCDsurfaceFormatCore
{
    // Private Member Functions

        static inline void writeShell
        (
            Ostream& os,
            const Face& f,
            const label cellId,
            const label cellTableId
        );

        //- Disallow default bitwise copy construct
        STARCDsurfaceFormat(const STARCDsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const STARCDsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct from file name
        STARCDsurfaceFormat(const fileName& filename);


    //- Destructor
    virtual ~STARCDsurfaceFormat() = default;


    // Static Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
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
            write(name, MeshedSurfaceProxy<Face>(*this), options);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "STARCDsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
