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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::foamVtkOutput::internalWriter::beginPiece()
{
    if (!legacy_)
    {
        format()
            .openTag(vtkFileTag::PIECE)
            ( "NumberOfPoints", vtkCells_.nFieldPoints() )
            ( "NumberOfCells",  vtkCells_.nFieldCells() )
            .closeTag();
    }
}


void Foam::foamVtkOutput::internalWriter::writePoints()
{
    // payload size
    const uint64_t payLoad = (vtkCells_.nFieldPoints() * 3 * sizeof(float));

    if (legacy_)
    {
        legacy::beginPoints(os_, vtkCells_.nFieldPoints());
    }
    else
    {
        format()
            .tag(vtkFileTag::POINTS)
            .openDataArray<float,3>(vtkFileTag::POINTS)
            .closeTag();
    }

    format().writeSize(payLoad);

    foamVtkOutput::writeList(format(), mesh_.points());
    foamVtkOutput::writeList
    (
        format(),
        mesh_.cellCentres(),
        vtkCells_.addPointCellLabels()
    );

    format().flush();

    if (!legacy_)
    {
        format()
            .endDataArray()
            .endTag(vtkFileTag::POINTS);
    }
}


void Foam::foamVtkOutput::internalWriter::writeCellsLegacy()
{
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


void Foam::foamVtkOutput::internalWriter::writeCells()
{
    format().tag(vtkFileTag::CELLS);

    //
    // 'connectivity'
    //
    {
        const labelList& vertLabels = vtkCells_.vertLabels();
        const uint64_t payLoad = vertLabels.size() * sizeof(label);

        format().openDataArray<label>("connectivity")
            .closeTag();

        format().writeSize(payLoad);
        foamVtkOutput::writeList(format(), vertLabels);
        format().flush();

        format().endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        const labelList& vertOffsets = vtkCells_.vertOffsets();
        const uint64_t payLoad = vertOffsets.size() * sizeof(label);

        format().openDataArray<label>("offsets")
            .closeTag();

        format().writeSize(payLoad);
        foamVtkOutput::writeList(format(), vertOffsets);
        format().flush();

        format().endDataArray();
    }


    //
    // 'types' (cell types)
    //
    {
        const List<uint8_t>& cellTypes = vtkCells_.cellTypes();
        const uint64_t payLoad = cellTypes.size() * sizeof(uint8_t);

        format().openDataArray<uint8_t>("types")
            .closeTag();

        format().writeSize(payLoad);
        forAll(cellTypes, i)
        {
            // No nComponents for char, cannot use foamVtkOutput::writeList here
            format().write(cellTypes[i]);
        }
        format().flush();

        format().endDataArray();
    }


    //
    // can quit here if there are NO face streams
    //
    if (vtkCells_.faceLabels().empty())
    {
        format().endTag(vtkFileTag::CELLS);

        return;
    }


    // --------------------------------------------------

    //
    // 'faces' (face streams)
    //
    {
        const labelList& faceLabels = vtkCells_.faceLabels();
        const uint64_t payLoad = faceLabels.size() * sizeof(label);

        format().openDataArray<label>("faces")
            .closeTag();

        format().writeSize(payLoad);
        foamVtkOutput::writeList(format(), faceLabels);
        format().flush();

        format().endDataArray();
    }


    // 'faceoffsets' (face stream offsets)
    // -1 to indicate that the cell is a primitive type that does not
    // have a face stream
    {
        const labelList& faceOffsets = vtkCells_.faceOffsets();
        const uint64_t payLoad = faceOffsets.size() * sizeof(label);

        format().openDataArray<label>("faceoffsets")
            .closeTag();

        format().writeSize(payLoad);
        foamVtkOutput::writeList(format(), faceOffsets);
        format().flush();

        format().endDataArray();
    }

    format().endTag(vtkFileTag::CELLS);
}


void Foam::foamVtkOutput::internalWriter::writeMesh()
{
    writePoints();
    if (legacy_)
    {
        writeCellsLegacy();
    }
    else
    {
        writeCells();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkOutput::internalWriter::internalWriter
(
    const fvMesh& mesh,
    const foamVtkCells& cells,
    const fileName& baseName,
    const foamVtkOutput::outputOptions outOpts
)
:
    mesh_(mesh),
    legacy_(outOpts.legacy()),
    format_(),
    vtkCells_(cells),
    os_()
{
    outputOptions opts(outOpts);
    opts.append(false);  // No append supported

    os_.open((baseName + (legacy_ ? ".vtk" : ".vtu")).c_str());
    format_ = opts.newFormatter(os_);

    const auto& title = mesh_.time().caseName();

    if (legacy_)
    {
        legacy::fileHeader(format(), title, vtkFileTag::UNSTRUCTURED_GRID);
    }
    else
    {
        // XML (inline)

        format()
            .xmlHeader()
            .xmlComment(title)
            .beginVTKFile(vtkFileTag::UNSTRUCTURED_GRID, "0.1");
    }

    beginPiece();
    writeMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkOutput::internalWriter::~internalWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::foamVtkOutput::internalWriter::beginCellData(label nFields)
{
    if (legacy_)
    {
        legacy::dataHeader
        (
            os(),
            vtkFileTag::CELL_DATA,
            vtkCells_.nFieldCells(),
            nFields
        );
    }
    else
    {
        format().tag(vtkFileTag::CELL_DATA);
    }
}


void Foam::foamVtkOutput::internalWriter::endCellData()
{
    if (!legacy_)
    {
        format().endTag(vtkFileTag::CELL_DATA);
    }
}


void Foam::foamVtkOutput::internalWriter::beginPointData(label nFields)
{
    if (legacy_)
    {
        legacy::dataHeader
        (
            os(),
            vtkFileTag::POINT_DATA,
            vtkCells_.nFieldPoints(),
            nFields
        );
    }
    else
    {
        format().tag(vtkFileTag::POINT_DATA);
    }
}


void Foam::foamVtkOutput::internalWriter::endPointData()
{
    if (!legacy_)
    {
        format().endTag(vtkFileTag::POINT_DATA);
    }
}


void Foam::foamVtkOutput::internalWriter::writeFooter()
{
    if (!legacy_)
    {
        // slight cheat. </Piece> too
        format().endTag(vtkFileTag::PIECE);

        format().endTag(vtkFileTag::UNSTRUCTURED_GRID)
            .endVTKFile();
    }
}


void Foam::foamVtkOutput::internalWriter::writeCellIDs()
{
    // Cell ids first
    const labelList& cellMap = vtkCells_.cellMap();
    const uint64_t payLoad = vtkCells_.nFieldCells() * sizeof(label);

    if (legacy_)
    {
        os_ << "cellID 1 " << vtkCells_.nFieldCells() << " int" << nl;
    }
    else
    {
        format().openDataArray<label>("cellID")
            .closeTag();
    }

    format().writeSize(payLoad);

    foamVtkOutput::writeList(format(), cellMap);
    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


// ************************************************************************* //
