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

InNamespace
    Foam::foam

Description
    Enumerations for commonly used VTK file items.

SourceFiles
    foamVtkFileEnums.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkFileEnums_H
#define foamVtkFileEnums_H

#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Some common XML tags for vtk files
enum class vtkFileTag
{
    VTK_FILE,           //!< "VTKFile"
    DATA_ARRAY,         //!< "DataArray"
    PIECE,              //!< "Piece"
    DATA_SET,           //!< "DataSet"
    POINTS,             //!< "Points"
    CELLS,              //!< "Cells"
    POLYS,              //!< "Polys"
    VERTS,              //!< "Verts"
    CELL_DATA,          //!< "CellData"
    POINT_DATA,         //!< "PointData"
    POLY_DATA,          //!< "PolyData"
    UNSTRUCTURED_GRID,  //!< "UnstructuredGrid"
};

//- Strings corresponding to the elements
extern const Foam::Enum<vtkFileTag> vtkFileTagNames;


//- Some common XML attributes for vtk files
enum class vtkFileAttr
{
    OFFSET,             //!< "offset"
    NUMBER_OF_POINTS,   //!< "NumberOfPoints"
    NUMBER_OF_CELLS,    //!< "NumberOfCells"
    NUMBER_OF_POLYS,    //!< "NumberOfPolys"
    NUMBER_OF_VERTS,    //!< "NumberOfVerts"
};

//- Strings corresponding to the elements
extern const Foam::Enum<vtkFileAttr> vtkFileAttrNames;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
