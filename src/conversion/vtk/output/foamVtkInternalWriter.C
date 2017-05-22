/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "foamVtkInternalWriter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkOutput::internalWriter::internalWriter
(
    const fvMesh& mesh,
    enum foamVtkOutput::formatType fmtType,
    const foamVtkCells& cells,
    const fileName& outputName
)
:
    mesh_(mesh),
    format_(),
    vtkCells_(cells),
    os_(outputName.c_str())
{
    format_ = foamVtkOutput::newFormatter(os_, fmtType);

    // Write header
    foamVtkOutput::legacy::fileHeader(format(), mesh.time().caseName())
        << "DATASET UNSTRUCTURED_GRID" << nl;

    //------------------------------------------------------------------
    //
    // Write topology
    //
    //------------------------------------------------------------------

    os_ << "POINTS " << vtkCells_.nFieldPoints() << " float" << nl;
    foamVtkOutput::writeList(format(), mesh.points());

    const pointField& ctrs = mesh.cellCentres();
    foamVtkOutput::writeList(format(), ctrs, vtkCells_.addPointCellLabels());

    format().flush();

    //
    // Write cells
    //

    const List<uint8_t>& cellTypes = vtkCells_.cellTypes();
    const labelList& vertLabels = vtkCells_.vertLabels();

    os_ << "CELLS " << vtkCells_.nFieldCells() << ' '
        << vertLabels.size() << nl;

    foamVtkOutput::writeList(format(), vertLabels);
    format().flush();

    os_ << "CELL_TYPES " << cellTypes.size() << nl;

    // No nComponents for char, so cannot use foamVtkOutput::writeList
    forAll(cellTypes, i)
    {
        format().write(cellTypes[i]);
    }
    format().flush();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::foamVtkOutput::internalWriter::beginCellData(label nFields)
{
    foamVtkOutput::legacy::cellDataHeader
    (
        os(),
        vtkCells_.nFieldCells(),
        nFields
    );
}


void Foam::foamVtkOutput::internalWriter::endCellData()
{}


void Foam::foamVtkOutput::internalWriter::beginPointData(label nFields)
{
    foamVtkOutput::legacy::pointDataHeader
    (
        os(),
        vtkCells_.nFieldPoints(),
        nFields
    );
}


void Foam::foamVtkOutput::internalWriter::endPointData()
{}


void Foam::foamVtkOutput::internalWriter::writeCellIDs()
{
    const labelList& cellMap = vtkCells_.cellMap();

    // Cell ids first
    os_ << "cellID 1 " << vtkCells_.nFieldCells() << " int" << nl;

    foamVtkOutput::writeList(format(), cellMap);
    format().flush();
}


// ************************************************************************* //
