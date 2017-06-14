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

void Foam::vtk::internalWriter::beginPiece()
{
    if (!legacy_)
    {
        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, vtuCells_.nFieldPoints())
            .xmlAttr(vtk::fileAttr::NUMBER_OF_CELLS,  vtuCells_.nFieldCells())
            .closeTag();
    }
}


void Foam::vtk::internalWriter::writePoints()
{
    // payload size
    const uint64_t payLoad = (vtuCells_.nFieldPoints() * 3 * sizeof(float));

    if (legacy_)
    {
        legacy::beginPoints(os_, vtuCells_.nFieldPoints());
    }
    else
    {
        format()
            .tag(vtk::fileTag::POINTS)
            .openDataArray<float,3>(vtk::dataArrayAttr::POINTS)
            .closeTag();
    }

    format().writeSize(payLoad);

    vtk::writeList(format(), mesh_.points());
    vtk::writeList
    (
        format(),
        mesh_.cellCentres(),
        vtuCells_.addPointCellLabels()
    );

    format().flush();

    if (!legacy_)
    {
        format()
            .endDataArray()
            .endTag(vtk::fileTag::POINTS);
    }
}


void Foam::vtk::internalWriter::writeCellsLegacy()
{
    const List<uint8_t>& cellTypes = vtuCells_.cellTypes();
    const labelList& vertLabels = vtuCells_.vertLabels();

    os_ << "CELLS " << vtuCells_.nFieldCells() << ' '
        << vertLabels.size() << nl;

    vtk::writeList(format(), vertLabels);
    format().flush();

    os_ << "CELL_TYPES " << cellTypes.size() << nl;

    // No nComponents for char, so cannot use vtk::writeList
    forAll(cellTypes, i)
    {
        format().write(cellTypes[i]);
    }
    format().flush();
}


void Foam::vtk::internalWriter::writeCells()
{
    format().tag(vtk::fileTag::CELLS);

    //
    // 'connectivity'
    //
    {
        const labelList& vertLabels = vtuCells_.vertLabels();
        const uint64_t payLoad = vertLabels.size() * sizeof(label);

        format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
            .closeTag();

        format().writeSize(payLoad);
        vtk::writeList(format(), vertLabels);
        format().flush();

        format().endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        const labelList& vertOffsets = vtuCells_.vertOffsets();
        const uint64_t payLoad = vertOffsets.size() * sizeof(label);

        format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
            .closeTag();

        format().writeSize(payLoad);
        vtk::writeList(format(), vertOffsets);
        format().flush();

        format().endDataArray();
    }


    //
    // 'types' (cell types)
    //
    {
        const List<uint8_t>& cellTypes = vtuCells_.cellTypes();
        const uint64_t payLoad = cellTypes.size() * sizeof(uint8_t);

        format().openDataArray<uint8_t>(vtk::dataArrayAttr::TYPES)
            .closeTag();

        format().writeSize(payLoad);
        forAll(cellTypes, i)
        {
            // No nComponents for char, cannot use vtk::writeList here
            format().write(cellTypes[i]);
        }
        format().flush();

        format().endDataArray();
    }


    //
    // can quit here if there are NO face streams
    //
    if (vtuCells_.faceLabels().empty())
    {
        format().endTag(vtk::fileTag::CELLS);

        return;
    }


    // --------------------------------------------------

    //
    // 'faces' (face streams)
    //
    {
        const labelList& faceLabels = vtuCells_.faceLabels();
        const uint64_t payLoad = faceLabels.size() * sizeof(label);

        format().openDataArray<label>(vtk::dataArrayAttr::FACES)
            .closeTag();

        format().writeSize(payLoad);
        vtk::writeList(format(), faceLabels);
        format().flush();

        format().endDataArray();
    }


    // 'faceoffsets' (face stream offsets)
    // -1 to indicate that the cell is a primitive type that does not
    // have a face stream
    {
        const labelList& faceOffsets = vtuCells_.faceOffsets();
        const uint64_t payLoad = faceOffsets.size() * sizeof(label);

        format().openDataArray<label>(vtk::dataArrayAttr::FACEOFFSETS)
            .closeTag();

        format().writeSize(payLoad);
        vtk::writeList(format(), faceOffsets);
        format().flush();

        format().endDataArray();
    }

    format().endTag(vtk::fileTag::CELLS);
}


void Foam::vtk::internalWriter::writeMesh()
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

Foam::vtk::internalWriter::internalWriter
(
    const fvMesh& mesh,
    const vtk::vtuCells& cells,
    const fileName& baseName,
    const vtk::outputOptions outOpts
)
:
    mesh_(mesh),
    legacy_(outOpts.legacy()),
    format_(),
    vtuCells_(cells),
    os_()
{
    outputOptions opts(outOpts);
    opts.append(false);  // No append supported

    os_.open((baseName + (legacy_ ? ".vtk" : ".vtu")).c_str());
    format_ = opts.newFormatter(os_);

    const auto& title = mesh_.time().caseName();

    if (legacy_)
    {
        legacy::fileHeader(format(), title, vtk::fileTag::UNSTRUCTURED_GRID);
    }
    else
    {
        // XML (inline)

        format()
            .xmlHeader()
            .xmlComment(title)
            .beginVTKFile(vtk::fileTag::UNSTRUCTURED_GRID, "0.1");
    }

    beginPiece();
    writeMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::internalWriter::~internalWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::vtk::internalWriter::beginCellData(label nFields)
{
    if (legacy_)
    {
        legacy::dataHeader
        (
            os(),
            vtk::fileTag::CELL_DATA,
            vtuCells_.nFieldCells(),
            nFields
        );
    }
    else
    {
        format().tag(vtk::fileTag::CELL_DATA);
    }
}


void Foam::vtk::internalWriter::endCellData()
{
    if (!legacy_)
    {
        format().endTag(vtk::fileTag::CELL_DATA);
    }
}


void Foam::vtk::internalWriter::beginPointData(label nFields)
{
    if (legacy_)
    {
        legacy::dataHeader
        (
            os(),
            vtk::fileTag::POINT_DATA,
            vtuCells_.nFieldPoints(),
            nFields
        );
    }
    else
    {
        format().tag(vtk::fileTag::POINT_DATA);
    }
}


void Foam::vtk::internalWriter::endPointData()
{
    if (!legacy_)
    {
        format().endTag(vtk::fileTag::POINT_DATA);
    }
}


void Foam::vtk::internalWriter::writeFooter()
{
    if (!legacy_)
    {
        // slight cheat. </Piece> too
        format().endTag(vtk::fileTag::PIECE);

        format().endTag(vtk::fileTag::UNSTRUCTURED_GRID)
            .endVTKFile();
    }
}


void Foam::vtk::internalWriter::writeCellIDs()
{
    // Cell ids first
    const labelList& cellMap = vtuCells_.cellMap();
    const uint64_t payLoad = vtuCells_.nFieldCells() * sizeof(label);

    if (legacy_)
    {
        os_ << "cellID 1 " << vtuCells_.nFieldCells() << " int" << nl;
    }
    else
    {
        format().openDataArray<label>("cellID")
            .closeTag();
    }

    format().writeSize(payLoad);

    vtk::writeList(format(), cellMap);
    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


// ************************************************************************* //
