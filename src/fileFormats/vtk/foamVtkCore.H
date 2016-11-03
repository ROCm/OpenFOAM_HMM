/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::foamVtkCore

Description
    Core routines for dealing with VTK files.

SourceFiles
    foamVtkCore.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkCore_H
#define foamVtkCore_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                  Class fileFormats::foamVtkCore Declaration
\*---------------------------------------------------------------------------*/

class foamVtkCore
{
public:

    // Public Data, Declarations

        //- The context when outputting a VTK file (XML or legacy).
        enum OutputContext
        {
            INLINE,     //<! Generate header and inline data
            HEADER,     //<! Generate header only
            APPEND      //<! Generate append-data
        };


        //- Equivalent to enumeration in "vtkCellType.h"
        enum vtkTypes
        {
            VTK_EMPTY_CELL       = 0,
            VTK_VERTEX           = 1,
            VTK_POLY_VERTEX      = 2,
            VTK_LINE             = 3,
            VTK_POLY_LINE        = 4,
            VTK_TRIANGLE         = 5,
            VTK_TRIANGLE_STRIP   = 6,
            VTK_POLYGON          = 7,
            VTK_PIXEL            = 8,
            VTK_QUAD             = 9,
            VTK_TETRA            = 10,
            VTK_VOXEL            = 11,
            VTK_HEXAHEDRON       = 12,
            VTK_WEDGE            = 13,
            VTK_PYRAMID          = 14,
            VTK_PENTAGONAL_PRISM = 15,
            VTK_HEXAGONAL_PRISM  = 16,
            VTK_POLYHEDRON       = 42
        };


protected:

    // Constructors

        //- Construct null
        foamVtkCore();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
