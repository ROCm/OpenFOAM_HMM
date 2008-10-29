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
    Foam::fileFormats::GTSfileFormat

Description
    Provide a means of reading/writing GTS format.

SourceFiles
    GTSfileFormat.C

\*---------------------------------------------------------------------------*/

#ifndef GTSfileFormat_H
#define GTSfileFormat_H

#include "Ostream.H"
#include "OFstream.H"
#include "keyedSurface.H"
#include "meshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

class GTSfileFormat;

/*---------------------------------------------------------------------------*\
                      Class GTSfileFormat Declaration
\*---------------------------------------------------------------------------*/

class GTSfileFormat
:
    public keyedSurface
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GTSfileFormat(const GTSfileFormat&);

        //- Disallow default bitwise assignment
        void operator=(const GTSfileFormat&);

public:

    // Constructors

        //- Construct null
        GTSfileFormat();

        //- Construct from file name
        GTSfileFormat(const fileName&, const bool triangulate=true);

    // Selectors

        //- Read file and return keyedSurface
        static autoPtr<keyedSurface> New
        (
            const fileName& fName,
            const bool triangulate=true
        )
        {
            return autoPtr<keyedSurface>
            (
                new GTSfileFormat(fName, triangulate)
            );
        }

    // Destructor

    // Member Functions

        // Write

        //- Write keyedSurface
        static void write(Ostream&, const keyedSurface&);

        //- Write keyedSurface
        static void write(const fileName& fName, const keyedSurface& surf)
        {
            write(OFstream(fName)(), surf);
        }

        //- Write meshedSurface
        static void write(Ostream&, const meshedSurface&);

        //- Write meshedSurface
        static void write(const fileName& fName, const meshedSurface& surf)
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
