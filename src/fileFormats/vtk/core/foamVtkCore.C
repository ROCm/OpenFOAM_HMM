/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "foamVtkCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::vtk::fileTag
>
Foam::vtk::fileTagNames
{
    { fileTag::VTK_FILE, "VTKFile" },
    { fileTag::DATA_ARRAY, "DataArray" },
    { fileTag::PIECE, "Piece" },
    { fileTag::DATA_SET, "DataSet" },
    { fileTag::POINTS, "Points" },
    { fileTag::CELLS, "Cells" },
    { fileTag::POLYS, "Polys" },
    { fileTag::VERTS, "Verts" },
    { fileTag::LINES, "Lines" },
    { fileTag::CELL_DATA, "CellData" },
    { fileTag::POINT_DATA, "PointData" },
    { fileTag::POLY_DATA, "PolyData" },
    { fileTag::UNSTRUCTURED_GRID, "UnstructuredGrid" },
};


const Foam::Enum
<
    Foam::vtk::fileAttr
>
Foam::vtk::fileAttrNames
{
    { fileAttr::OFFSET, "offset" },
    { fileAttr::NUMBER_OF_COMPONENTS, "NumberOfComponents" },
    { fileAttr::NUMBER_OF_POINTS, "NumberOfPoints" },
    { fileAttr::NUMBER_OF_CELLS, "NumberOfCells" },
    { fileAttr::NUMBER_OF_POLYS, "NumberOfPolys" },
    { fileAttr::NUMBER_OF_VERTS, "NumberOfVerts" },
    { fileAttr::NUMBER_OF_LINES, "NumberOfLines" },
};


const Foam::Enum
<
    Foam::vtk::dataArrayAttr
>
Foam::vtk::dataArrayAttrNames
{
    { dataArrayAttr::POINTS, "Points" },
    { dataArrayAttr::OFFSETS, "offsets" },
    { dataArrayAttr::CONNECTIVITY, "connectivity" },
    { dataArrayAttr::TYPES, "types" },
    { dataArrayAttr::FACES, "faces" },
    { dataArrayAttr::FACEOFFSETS, "faceoffsets" },
};


// ************************************************************************* //
