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

Namespace
    Foam::vtk

Description
    Namespace for handling VTK input/output.

SourceFiles
    foamVtkCore.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkCore_H
#define foamVtkCore_H

#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // Enumerations

    //- The context when outputting a VTK file (XML or legacy).
    enum OutputContext
    {
        INLINE,     //<! Generate header and inline data
        HEADER,     //<! Generate header only
        APPEND      //<! Generate append-data
    };


    //- The output format type for file contents.
    //  Upper bits for output type, lower bits for the format itself.
    enum class formatType
    {
        /** XML inline ASCII, using the asciiFormatter */
        INLINE_ASCII  = 0,
        /** XML inline base64, using the base64Formatter */
        INLINE_BASE64 = 0x01,
        /** XML append base64, using the appendBase64Formatter */
        APPEND_BASE64 = 0x11,
        /** XML append raw binary, using the appendRawFormatter */
        APPEND_BINARY = 0x12,
        /** Legacy ASCII, using the legacyAsciiFormatter */
        LEGACY_ASCII  = 0x20,
        /** Legacy raw binary, using the legacyRawFormatter */
        LEGACY_BINARY = 0x22,
    };

    //- Equivalent to enumeration in "vtkCellType.h"
    enum cellType
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


    //- Some common XML tags for vtk files
    enum class fileTag
    {
        VTK_FILE,               //!< "VTKFile"
        DATA_ARRAY,             //!< "DataArray"
        PIECE,                  //!< "Piece"
        DATA_SET,               //!< "DataSet"
        POINTS,                 //!< "Points"
        CELLS,                  //!< "Cells"
        POLYS,                  //!< "Polys"
        VERTS,                  //!< "Verts"
        LINES,                  //!< "Lines"
        CELL_DATA,              //!< "CellData"
        POINT_DATA,             //!< "PointData"
        POLY_DATA,              //!< "PolyData"
        UNSTRUCTURED_GRID,      //!< "UnstructuredGrid"
    };

    //- Strings corresponding to the vtk xml tags
    extern const Foam::Enum<fileTag> fileTagNames;

    //- Some common XML attributes for vtk files
    enum class fileAttr
    {
        OFFSET,                 //!< "offset"
        NUMBER_OF_COMPONENTS,   //!< "NumberOfComponents"
        NUMBER_OF_POINTS,       //!< "NumberOfPoints"
        NUMBER_OF_CELLS,        //!< "NumberOfCells"
        NUMBER_OF_POLYS,        //!< "NumberOfPolys"
        NUMBER_OF_VERTS,        //!< "NumberOfVerts"
        NUMBER_OF_LINES,        //!< "NumberOfLines"
    };

    //- Strings corresponding to the vtk xml attributes
    extern const Foam::Enum<fileAttr> fileAttrNames;

    //- Some common names for XML data arrays
    enum class dataArrayAttr
    {
        POINTS,                 //!< "Points"
        OFFSETS,                //!< "offsets"
        CONNECTIVITY,           //!< "connectivity"
        TYPES,                  //!< "types"
        FACES,                  //!< "faces"
        FACEOFFSETS,            //!< "faceoffsets"
    };

    //- Strings corresponding to the vtk xml attributes
    extern const Foam::Enum<dataArrayAttr> dataArrayAttrNames;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
