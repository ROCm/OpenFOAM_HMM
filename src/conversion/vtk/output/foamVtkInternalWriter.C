/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "foamVtkInternalWriter.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::vtk::internalWriter::debug = 0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::internalWriter::beginPiece()
{
    // Basic sizes

    numberOfPoints_ = vtuCells_.nFieldPoints(); // With addPointCellLabels
    numberOfCells_  = vtuCells_.nFieldCells();  // With decomposed cells

    if (parallel_)
    {
        reduce(numberOfPoints_, sumOp<label>());
        reduce(numberOfCells_,  sumOp<label>());
    }

    // Nothing else to do for legacy
    if (legacy()) return;

    DebugInFunction
        << "nPoints=" << numberOfPoints_ << " nCells=" << numberOfCells_ << nl;

    if (format_)
    {
        format()
            .tag
            (
                vtk::fileTag::PIECE,
                vtk::fileAttr::NUMBER_OF_POINTS, numberOfPoints_,
                vtk::fileAttr::NUMBER_OF_CELLS,  numberOfCells_
            );
    }
}


void Foam::vtk::internalWriter::writePoints()
{
    if (format_)
    {
        if (legacy())
        {
            legacy::beginPoints(os_, numberOfPoints_);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float,3>(numberOfPoints_);

            format()
                .tag(vtk::fileTag::POINTS)
                .beginDataArray<float,3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);
        }
    }


    if (parallel_)
    {
        vtk::writeListsParallel
        (
            format_.ref(),
            mesh_.points(),
            mesh_.cellCentres(),
            vtuCells_.addPointCellLabels()
        );
    }
    else
    {
        vtk::writeLists
        (
            format(),
            mesh_.points(),
            mesh_.cellCentres(),
            vtuCells_.addPointCellLabels()
        );
    }


    if (format_)
    {
        format().flush();
        format().endDataArray();

        if (!legacy())
        {
            format()
                .endTag(vtk::fileTag::POINTS);
        }
    }
}


void Foam::vtk::internalWriter::writeCellsLegacy(const label pointOffset)
{
    const List<uint8_t>& cellTypes = vtuCells_.cellTypes();
    const labelList& vertLabels = vtuCells_.vertLabels();

    label nCells = cellTypes.size();
    label nVerts = vertLabels.size();

    if (parallel_)
    {
        reduce(nCells, sumOp<label>());
        reduce(nVerts, sumOp<label>());
    }

    if (nCells != numberOfCells_)
    {
        FatalErrorInFunction
            << "Expecting " << numberOfCells_
            << " cells, but found " << nCells
            << exit(FatalError);
    }


    // CELLS
    {
        if (format_)
        {
            os_ << nl
                << "CELLS " << nCells << ' ' << nVerts << nl;
        }

        if (parallel_)
        {
            vtk::writeListParallel
            (
                format_.ref(),
                vtk::vtuSizing::copyVertLabelsLegacy
                (
                    vertLabels,
                    pointOffset
                )
            );
        }
        else
        {
            vtk::writeList(format(), vertLabels);
        }

        if (format_)
        {
            format().flush();
        }
    }


    // CELL_TYPES
    {
        if (format_)
        {
            os_ << nl
                << "CELL_TYPES " << nCells << nl;
        }

        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), cellTypes);
        }
        else
        {
            vtk::writeList(format(), cellTypes);
        }

        if (format_)
        {
            format().flush();
        }
    }
}


void Foam::vtk::internalWriter::writeCellsConnectivity(const label pointOffset)
{
    //
    // 'connectivity'
    //
    {
        const labelList& vertLabels = vtuCells_.vertLabels();
        label nVerts = vertLabels.size();

        if (parallel_)
        {
            reduce(nVerts, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nVerts);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);
        }

        if (parallel_)
        {
            vtk::writeListParallel
            (
                format_.ref(),
                vtk::vtuSizing::copyVertLabelsXml
                (
                    vertLabels,
                    pointOffset
                )
            );
        }
        else
        {
            vtk::writeList(format(), vertLabels);
        }

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        const labelList& vertOffsets = vtuCells_.vertOffsets();
        label nOffs = vertOffsets.size();

        if (parallel_)
        {
            reduce(nOffs, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad =
                vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }

        if (parallel_)
        {
            // processor-local connectivity offsets
            const globalIndex procOffset
            (
                vertOffsets.empty() ? 0 : vertOffsets.last()
            );

            vtk::writeListParallel(format_.ref(), vertOffsets, procOffset);
        }
        else
        {
            vtk::writeList(format(), vertOffsets);
        }

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }


    //
    // 'types' (cell types)
    //
    {
        const List<uint8_t>& cellTypes = vtuCells_.cellTypes();
        label nCells = cellTypes.size();

        if (parallel_)
        {
            reduce(nCells, sumOp<label>());
        }

        if (nCells != numberOfCells_)
        {
            FatalErrorInFunction
                << "Expecting " << numberOfCells_
                << " cells, but found " << nCells
                << exit(FatalError);
        }

        if (format_)
        {
            const uint64_t payLoad =
                vtk::sizeofData<uint8_t>(nCells);

            format().beginDataArray<uint8_t>(vtk::dataArrayAttr::TYPES);
            format().writeSize(payLoad);
        }

        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), cellTypes);
        }
        else
        {
            vtk::writeList(format(), cellTypes);
        }

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }
}


void Foam::vtk::internalWriter::writeCellsFaces(const label pointOffset)
{
    label nFaceLabels = vtuCells_.faceLabels().size();

    if (parallel_)
    {
        reduce(nFaceLabels, sumOp<label>());
    }

    // Can quit now if there are NO face streams
    if (!nFaceLabels)
    {
        return;
    }

    // --------------------------------------------------

    //
    // 'faces' (face streams)
    //
    const labelList& faceLabels = vtuCells_.faceLabels();

    {
        // Already have nFaceLabels (above)

        if (format_)
        {
            const uint64_t payLoad =
                vtk::sizeofData<label>(nFaceLabels);

            format().beginDataArray<label>(vtk::dataArrayAttr::FACES);
            format().writeSize(payLoad);
        }


        if (parallel_)
        {
            vtk::writeListParallel
            (
                format_.ref(),
                vtk::vtuSizing::copyFaceLabelsXml
                (
                    faceLabels,
                    pointOffset
                )
            );
        }
        else
        {
            vtk::writeList(format(), faceLabels);
        }


        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }

    // 'faceoffsets' (face stream offsets)
    // -1 to indicate that the cell is a primitive type that does not
    // have a face stream

    // If the processor-local mesh has any polyhedrals, we have a list with
    // the faceoffsets and we just need to renumber.
    // If the processor-local mesh has NO polyhedrals (but others do), we
    // need to generate a list of -1 for that processor.
    //
    // Result: A face offset value for each cell.
    {
        const labelList& faceOffsets = vtuCells_.faceOffsets();
        const label nLocalCells = vtuCells_.cellTypes().size();

        label nCells = nLocalCells;

        if (parallel_)
        {
            reduce(nCells, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad =
                vtk::sizeofData<label>(nCells);

            format().beginDataArray<label>(vtk::dataArrayAttr::FACEOFFSETS);
            format().writeSize(payLoad);
        }


        if (parallel_)
        {
            const List<uint8_t>& cellTypes = vtuCells_.cellTypes();
            const label nLocalCells = cellTypes.size();

            const globalIndex procOffset(faceLabels.size());

            labelList faceOffsetsRenumber;

            if (faceOffsets.size()) // Or check procOffset.localSize()
            {
                faceOffsetsRenumber =
                    vtk::vtuSizing::copyFaceOffsetsXml
                    (
                        faceOffsets,
                        procOffset.localStart()
                    );
            }
            else
            {
                faceOffsetsRenumber.resize(nLocalCells, -1);
            }

            vtk::writeListParallel(format_.ref(), faceOffsetsRenumber);
        }
        else
        {
            vtk::writeList(format(), faceOffsets);
        }


        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::internalWriter::internalWriter
(
    const fvMesh& mesh,
    const vtk::vtuCells& cells,
    const vtk::outputOptions opts
)
:
    vtk::fileWriter(vtk::fileTag::UNSTRUCTURED_GRID, opts),
    mesh_(mesh),
    vtuCells_(cells),
    numberOfPoints_(0),
    numberOfCells_(0)
{
    // We do not currently support append mode
    opts_.append(false);
}


Foam::vtk::internalWriter::internalWriter
(
    const fvMesh& mesh,
    const vtk::vtuCells& cells,
    const fileName& file,
    bool parallel
)
:
    internalWriter(mesh, cells)
{
    open(file, parallel);
}


Foam::vtk::internalWriter::internalWriter
(
    const fvMesh& mesh,
    const vtk::vtuCells& cells,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    internalWriter(mesh, cells, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::internalWriter::beginFile(std::string title)
{
    if (title.size())
    {
        return vtk::fileWriter::beginFile(title);
    }

    // Provide default title

    DebugInFunction
        << "case=" << mesh_.time().caseName()
        << " region=" << mesh_.name()
        << " time=" << mesh_.time().timeName()
        << " index=" << mesh_.time().timeIndex() << endl;


    if (legacy())
    {
        return vtk::fileWriter::beginFile
        (
            mesh_.time().globalCaseName()
        );
    }


    // XML (inline)

    return vtk::fileWriter::beginFile
    (
        "case='" + mesh_.time().globalCaseName()
      + "' region='" + mesh_.name()
      + "' time='" + mesh_.time().timeName()
      + "' index='" + Foam::name(mesh_.time().timeIndex())
      + "'"
    );
}


bool Foam::vtk::internalWriter::writeGeometry()
{
    enter_Piece();

    beginPiece();

    writePoints();

    // Include addPointCellLabels for the point offsets
    const label pointOffset =
    (
        parallel_ ? globalIndex(vtuCells_.nFieldPoints()).localStart() : 0
    );

    if (legacy())
    {
        writeCellsLegacy(pointOffset);
        return true;
    }

    if (format_)
    {
        format().tag(vtk::fileTag::CELLS);
    }

    writeCellsConnectivity(pointOffset);
    writeCellsFaces(pointOffset);

    if (format_)
    {
        format().endTag(vtk::fileTag::CELLS);
    }

    return true;
}


bool Foam::vtk::internalWriter::beginCellData(label nFields)
{
    return enter_CellData(numberOfCells_, nFields);
}


bool Foam::vtk::internalWriter::beginPointData(label nFields)
{
    return enter_PointData(numberOfPoints_, nFields);
}


void Foam::vtk::internalWriter::writeCellIDs()
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") for cellID field" << nl << endl
            << exit(FatalError);
    }

    const labelList& cellMap = vtuCells_.cellMap();

    if (format_)
    {
        if (legacy())
        {
            vtk::legacy::intField<1>(format(), "cellID", numberOfCells_);
        }
        else
        {
            const uint64_t payLoad = vtk::sizeofData<label>(numberOfCells_);

            format().beginDataArray<label>("cellID");
            format().writeSize(payLoad);
        }
    }


    if (parallel_)
    {
        // With decomposed cells for the cell offsets
        const globalIndex globalCellOffset(vtuCells_.nFieldCells());

        vtk::writeListParallel(format_.ref(), cellMap, globalCellOffset);
    }
    else
    {
        vtk::writeList(format(), cellMap);
    }

    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


bool Foam::vtk::internalWriter::writeProcIDs()
{
    if (!parallel_)
    {
        // Disabled in serial output (meaningless)
        return false;
    }

    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") for procID field" << nl << endl
            << exit(FatalError);
    }

    const globalIndex procMaps(vtuCells_.nFieldCells());

    bool good = false;

    if (Pstream::master())
    {
        const label nCells = procMaps.size();

        if (format_)
        {
            if (legacy())
            {
                vtk::legacy::intField<1>(format(), "procID", nCells);
            }
            else
            {
                const uint64_t payLoad =
                    vtk::sizeofData<label>(nCells);

                format().beginDataArray<label>("procID");
                format().writeSize(payLoad);
            }
        }

        // Per-processor ids
        for (label proci=0; proci < Pstream::nProcs(); ++proci)
        {
            label len = procMaps.localSize(proci);

            while (len--)
            {
                format().write(proci);
            }
        }

        format().flush();
        format().endDataArray();

        good = true;
    }

    // MPI barrier
    return returnReduce(good, orOp<bool>());
}


void Foam::vtk::internalWriter::writePointIDs()
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::POINT_DATA]
            << ") for pointID field" << nl << endl
            << exit(FatalError);
    }

    if (format_)
    {
        if (legacy())
        {
            vtk::legacy::intField<1>(format(), "pointID", numberOfPoints_);
        }
        else
        {
            const uint64_t payLoad = vtk::sizeofData<label>(numberOfPoints_);

            format().beginDataArray<label>("pointID");
            format().writeSize(payLoad);
        }
    }


    // Point offset for regular mesh points (without decomposed)
    const label pointOffset =
    (
        parallel_ ? globalIndex(vtuCells_.nPoints()).localStart() : 0
    );

    // Cell offset for *regular* mesh cells (without decomposed)
    const label cellOffset =
    (
        parallel_ ? globalIndex(vtuCells_.nCells()).localStart() : 0
    );


    labelList pointIds = identity(vtuCells_.nFieldPoints(), pointOffset);

    // The pointID for added points is the cellID, tag as a negative number
    label pointi = vtuCells_.nPoints();
    for (const label celli : vtuCells_.addPointCellLabels())
    {
        pointIds[pointi] = (-1 - celli - cellOffset);
        ++pointi;
    }

    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), pointIds);
    }
    else
    {
        vtk::writeList(format(), pointIds);
    }

    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


// ************************************************************************* //
