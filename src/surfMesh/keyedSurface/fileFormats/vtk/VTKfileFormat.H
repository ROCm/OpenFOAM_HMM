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
    Foam::fileFormats::VTKfileFormat

Description
    Provide a means of writing VTK legacy format.

SourceFiles
    VTKfileFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTKfileFormat_H
#define VTKfileFormat_H

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
                       Class VTKfileFormat Declaration
\*---------------------------------------------------------------------------*/

class VTKfileFormat
:
    public keyedSurface
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        VTKfileFormat(const VTKfileFormat&);

        //- Disallow default bitwise assignment
        void operator=(const VTKfileFormat&);

        static void writeHead(Ostream&, const pointField&, const List<face>&);
        static void writeTail(Ostream&, const List<surfacePatch>&);

public:

    // Constructors

        //- Construct null
        VTKfileFormat();

    // Selectors

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
