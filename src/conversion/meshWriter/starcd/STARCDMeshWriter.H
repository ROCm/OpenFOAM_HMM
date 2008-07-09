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
    Foam::meshWriters::STARCD

Description
    Writes polyMesh in pro-STAR (v4) bnd/cel/vrt format

    Alternatively, extracts the surface of the FOAM mesh into
    pro-STAR (v4) .cel/.vrt/ format.
    This can be useful, for example, for surface morphing in an external
    package.

    The cellTableId and cellTable information are used (if available).
    Otherwise the cellZones are used (if available).

SourceFiles
    STARCDMeshWriter.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDMeshWriter_H
#define STARCDMeshWriter_H

#include "meshWriter.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace meshWriters
{

/*---------------------------------------------------------------------------*\
                    Class meshWriters::STARCD Declaration
\*---------------------------------------------------------------------------*/

class STARCD
:
    public meshWriter
{
    // Private Data
        static const char* defaultBoundaryName;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        STARCD(const STARCD&);

        //- Disallow default bitwise assignment
        void operator=(const STARCD&);

        //- pro-STAR 4+ header format
        static void writeHeader(Ostream&, const char* filetype);

        //- write points
        void writePoints(const fileName& baseName) const;

        //- write cells
        void writeCells(const fileName& baseName) const;

        //- write boundary
        void writeBoundary(const fileName& baseName) const;

        void getCellTable();

        label findDefaultBoundary() const;

public:

    // Static data members

        //- Face addressing from OpenFOAM faces -> pro-STAR faces
        static const label foamToStarFaceAddr[4][6];


    // Constructors

        //- open a file for writing
        STARCD
        (
            const polyMesh&,
            const scalar scaleFactor = 1.0
        );

    // Selectors

    // Destructor

        virtual ~STARCD();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


        //- Remove STAR-CD files for the baseName
        void rmFiles(const fileName& baseName) const;

        //- write volume mesh
        virtual bool write
        (
            const fileName& meshName = fileName::null
        ) const;

        //- write surface mesh with optional triangulation
        virtual bool writeSurface
        (
            const fileName& meshName = fileName::null,
            const bool& triangulate = false
        ) const;


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshWriters

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
