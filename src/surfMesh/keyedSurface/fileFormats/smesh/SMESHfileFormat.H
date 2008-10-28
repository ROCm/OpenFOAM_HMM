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
    Foam::fileFormats::SMESHfileFormat

Description
    Provide a means of writing tetgen SMESH format.

    Tetgen http://tetgen.berlios.de

See Also
    File format information:
    http://tetgen.berlios.de/fformats.smesh.html

SourceFiles
    SMESHfileFormat.C

\*---------------------------------------------------------------------------*/

#ifndef SMESHfileFormat_H
#define SMESHfileFormat_H

#include "IFstream.H"
#include "OFstream.H"
#include "Ostream.H"
#include "keyedSurface.H"
#include "meshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class SMESHfileFormat Declaration
\*---------------------------------------------------------------------------*/

class SMESHfileFormat
:
    public keyedSurface
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        SMESHfileFormat(const SMESHfileFormat&);

        //- Disallow default bitwise assignment
        void operator=(const SMESHfileFormat&);

        static void writeHead(Ostream&, const pointField&, const List<face>&);
        static void writeTail(Ostream&);

public:

    // Constructors

        //- Construct null
        SMESHfileFormat();

    // Selectors

    // Destructor

    // Member Functions

        //- Write keyedSurface
        static void write(Ostream&, const keyedSurface&);

        //- Write meshedSurface
        static void write(Ostream&, const meshedSurface&);

        //- Write keyedSurface
        static void write(const fileName& fName, const keyedSurface& surf)
        {
            write(OFstream(fName)(), surf);
        }

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
