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
    Foam::fileFormats::NASfileFormat

Description
    Nastran surface reader.

    - Uses the Ansa "$ANSA_NAME" or the Hypermesh "$HMNAME COMP" extensions
      to obtain patch names.
    - Handles Nastran short and long formats, but not free format.
    - Properly handles the Nastran compact floating point notation: \n
    @verbatim
        GRID          28        10.20269-.030265-2.358-8
    @endverbatim

SourceFiles
    NASfileFormat.C

\*---------------------------------------------------------------------------*/

#ifndef NASfileFormat_H
#define NASfileFormat_H

#include "Ostream.H"
#include "OFstream.H"
#include "keyedSurface.H"
#include "meshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                       Class NASfileFormat Declaration
\*---------------------------------------------------------------------------*/

class NASfileFormat
:
    public keyedSurface
{

    // Private Member Functions

        //- Do weird things to extract number
        static scalar parseNASCoord(const string&);

        //- Disallow default bitwise copy construct
        NASfileFormat(const NASfileFormat&);

        //- Disallow default bitwise assignment
        void operator=(const NASfileFormat&);

public:

    // Constructors

        //- Construct null
        NASfileFormat();

        //- Construct from file name
        NASfileFormat(const fileName&, const bool triangulate=false);

    // Selectors

        //- Read file and return keyedSurface
        static autoPtr<keyedSurface> New
        (
            const fileName& fName,
            const bool triangulate=false
        )
        {
            return autoPtr<keyedSurface>
            (
                new NASfileFormat(fName,triangulate)
            );
        }

    // Destructor

    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
