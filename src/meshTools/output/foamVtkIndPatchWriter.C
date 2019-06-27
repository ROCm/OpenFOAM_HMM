/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "foamVtkIndPatchWriter.H"
#include "foamVtkOutput.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::indirectPatchWriter::beginPiece()
{
    // Basic sizes
    nLocalPoints_ = pp_.nPoints();
    nLocalFaces_  = pp_.size();
    nLocalVerts_  = 0;

    for (const face& f : pp_)
    {
        nLocalVerts_ += f.size();
    }

    numberOfPoints_ = nLocalPoints_;
    numberOfCells_  = nLocalFaces_;

    if (parallel_)
    {
        reduce(numberOfPoints_, sumOp<label>());
        reduce(numberOfCells_,  sumOp<label>());
    }


    // Nothing else to do for legacy
    if (legacy()) return;

    if (format_)
    {
        format().tag
        (
            vtk::fileTag::PIECE,
            vtk::fileAttr::NUMBER_OF_POINTS, numberOfPoints_,
            vtk::fileAttr::NUMBER_OF_POLYS,  numberOfCells_
        );
    }
}


void Foam::vtk::indirectPatchWriter::writePoints()
{
    if (format_)
    {
        if (legacy())
        {
            legacy::beginPoints(os_, numberOfPoints_);
        }
        else
        {
            const uint64_t payLoad = vtk::sizeofData<float, 3>(numberOfPoints_);

            format()
                .tag(vtk::fileTag::POINTS)
                .beginDataArray<float,3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);
        }
    }


    if (parallel_ ? Pstream::master() : true)
    {
        {
            vtk::writeList(format(), pp_.localPoints());
        }
    }


    if (parallel_)
    {
        if (Pstream::master())
        {
            pointField recv;

            // Receive each point field and write
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                ++slave
            )
            {
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);

                {
                    fromSlave >> recv;

                    vtk::writeList(format(), recv);
                }
            }
        }
        else
        {
            // Send to master
            OPstream toMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            {
                toMaster << pp_.localPoints();
            }
        }
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


void Foam::vtk::indirectPatchWriter::writePolysLegacy(const label pointOffset)
{
    // Connectivity count without additional storage (done internally)

    label nFaces = nLocalFaces_;
    label nVerts = nLocalVerts_;

    if (parallel_)
    {
        reduce(nFaces, sumOp<label>());
        reduce(nVerts, sumOp<label>());
    }

    if (nFaces != numberOfCells_)
    {
        FatalErrorInFunction
            << "Expecting " << numberOfCells_
            << " faces, but found " << nFaces
            << exit(FatalError);
    }

    legacy::beginPolys(os_, nFaces, nVerts);

    labelList vertLabels(nLocalFaces_ + nLocalVerts_);

    {
        // Legacy: size + connectivity together
        // [nPts, id1, id2, ..., nPts, id1, id2, ...]

        auto iter = vertLabels.begin();

        label off = pointOffset;

        {
            for (const face& f : pp_.localFaces())
            {
                *iter = f.size();       // The size prefix
                ++iter;

                for (const label pfi : f)
                {
                    *iter = pfi + off;  // Face vertex label
                    ++iter;
                }
            }
            // off += pp_.nPoints();
        }
    }


    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), vertLabels);
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


void Foam::vtk::indirectPatchWriter::writePolys(const label pointOffset)
{
    if (format_)
    {
        format().tag(vtk::fileTag::POLYS);
    }

    //
    // 'connectivity'
    //
    {
        labelList vertLabels(nLocalVerts_);

        label nVerts = nLocalVerts_;

        if (parallel_)
        {
            reduce(nVerts, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nVerts);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad * sizeof(label));
        }

        {
            // XML: connectivity only
            // [id1, id2, ..., id1, id2, ...]

            auto iter = vertLabels.begin();

            label off = pointOffset;

            {
                for (const face& f : pp_.localFaces())
                {
                    for (const label pfi : f)
                    {
                        *iter = pfi + off;  // Face vertex label
                        ++iter;
                    }
                }
                // off += pp_.nPoints();
            }
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertLabels);
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
        labelList vertOffsets(nLocalFaces_);
        label nOffs = vertOffsets.size();

        if (parallel_)
        {
            reduce(nOffs, sumOp<label>());
        }

        if (format_)
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }


        // processor-local connectivity offsets
        label off =
        (
            parallel_ ? globalIndex(nLocalVerts_).localStart() : 0
        );

        auto iter = vertOffsets.begin();

        {
            for (const face& f : pp_)
            {
                off += f.size();   // End offset
                *iter = off;
                ++iter;
            }
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertOffsets);
        }
        else
        {
            vtk::writeList(format_.ref(), vertOffsets);
        }


        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }

    if (format_)
    {
        format().endTag(vtk::fileTag::POLYS);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::indirectPatchWriter::indirectPatchWriter
(
    const indirectPrimitivePatch& pp,
    const vtk::outputOptions opts
)
:
    vtk::fileWriter(vtk::fileTag::POLY_DATA, opts),
    pp_(pp),
    numberOfPoints_(0),
    numberOfCells_(0),
    nLocalPoints_(0),
    nLocalFaces_(0),
    nLocalVerts_(0)
{
    // We do not currently support append mode
    opts_.append(false);
}


Foam::vtk::indirectPatchWriter::indirectPatchWriter
(
    const indirectPrimitivePatch& pp,
    const fileName& file,
    bool parallel
)
:
    indirectPatchWriter(pp)
{
    open(file, parallel);
}


Foam::vtk::indirectPatchWriter::indirectPatchWriter
(
    const indirectPrimitivePatch& pp,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    indirectPatchWriter(pp, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::indirectPatchWriter::beginFile(std::string title)
{
    if (title.size())
    {
        return vtk::fileWriter::beginFile(title);
    }

    // Provide default title
    return vtk::fileWriter::beginFile("surfaces");
}


bool Foam::vtk::indirectPatchWriter::writeGeometry()
{
    enter_Piece();

    beginPiece();

    writePoints();

    const label pointOffset =
    (
        parallel_ ? globalIndex(nLocalPoints_).localStart() : 0
    );

    if (legacy())
    {
        writePolysLegacy(pointOffset);
    }
    else
    {
        writePolys(pointOffset);
    }

    return true;
}


bool Foam::vtk::indirectPatchWriter::beginCellData(label nFields)
{
    return enter_CellData(numberOfCells_, nFields);
}


bool Foam::vtk::indirectPatchWriter::beginPointData(label nFields)
{
    return enter_PointData(numberOfPoints_, nFields);
}


// ************************************************************************* //
