/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::fileFormats::FLMAsurfaceFormat

Description
    Provide a means of writing AVL/FIRE FLMA format.

SourceFiles
    FLMAsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef FLMAsurfaceFormat_H
#define FLMAsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "FIRECore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
               Class fileFormats::FLMAsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class FLMAsurfaceFormat
:
    public MeshedSurface<Face>,
    public FIRECore
{
    // Private Member Functions

        static inline void writeShell(OSstream& os, const Face& f);
        static inline void writeType(OSstream& os, const Face& f);

        //- Disallow default bitwise copy construct
        FLMAsurfaceFormat(const FLMAsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const FLMAsurfaceFormat<Face>&) = delete;


protected:

    // Protected Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            OSstream& os,
            const MeshedSurfaceProxy<Face>& surf
        );


        //- Write surface mesh components by proxy with/without compression
        static void write
        (
            bool compress,
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf
        );


public:

    // Constructors

        //- Construct null
        FLMAsurfaceFormat() = default;


    //- Destructor
    virtual ~FLMAsurfaceFormat() = default;


    // Static Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Write flma file
        virtual void write
        (
            const fileName& name,
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this), options);
        }

};


/*---------------------------------------------------------------------------*\
                     Class FLMAZsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class FLMAZsurfaceFormat
:
    public FLMAsurfaceFormat<Face>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        FLMAZsurfaceFormat(const FLMAZsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const FLMAZsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct null
        FLMAZsurfaceFormat() = default;


    //- Destructor
    virtual ~FLMAZsurfaceFormat() = default;


    // Static Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Write flmaz file
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
#   include "FLMAsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
