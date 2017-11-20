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
    Foam::fileFormats::OFFsurfaceFormat

Description
    Read/write Geomview OFF polyList format.

See also
    The <a href="http://www.geoview.org">Geoview</a>
    file format information:
    http://www.geomview.org/docs/html/OFF.html#OFF

Note
    When reading, the optional \a colorspec is ignored.
    When writing, it is set to the zone number (integer).

SourceFiles
    OFFsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef OFFsurfaceFormat_H
#define OFFsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::OFFsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class OFFsurfaceFormat
:
    public MeshedSurface<Face>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        OFFsurfaceFormat(const OFFsurfaceFormat&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const OFFsurfaceFormat&) = delete;


public:

    // Constructors

        //- Construct from file name
        OFFsurfaceFormat(const fileName& filename);


    //- Destructor
    virtual ~OFFsurfaceFormat() = default;


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
    #include "OFFsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
