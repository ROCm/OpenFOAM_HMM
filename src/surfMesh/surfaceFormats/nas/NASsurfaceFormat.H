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
    Foam::fileFormats::NASsurfaceFormat

Description
    Nastran surface reader.

    - Uses the Ansa "$ANSA_NAME" or the Hypermesh "$HMNAME COMP" extensions
      to obtain zone names.
    - Handles Nastran short, long formats and comma-separated free format.
    - Properly handles the Nastran compact floating point notation: \n
    \verbatim
        GRID          28        10.20269-.030265-2.358-8
    \endverbatim

SourceFiles
    NASsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef NASsurfaceFormat_H
#define NASsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "NASCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::NASsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class NASsurfaceFormat
:
    public MeshedSurface<Face>,
    public NASCore
{
    // Private Member Functions

        //- Output CTRIA3 or CQUAD4
        inline static label writeShell
        (
            Ostream& os,
            const Face& f,
            const label groupId,
            label elementId
        );

        //- Disallow default bitwise copy construct
        NASsurfaceFormat(const NASsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const NASsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct from file name
        NASsurfaceFormat(const fileName& filename);


    //- Destructor
    virtual ~NASsurfaceFormat() = default;


    // Static Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Read from a file
        virtual bool read(const fileName& filename);

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
    #include "NASsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
