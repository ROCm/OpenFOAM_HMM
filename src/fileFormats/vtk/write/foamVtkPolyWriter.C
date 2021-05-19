/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "foamVtkPolyWriter.H"
#include "foamVtkOutput.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::polyWriter::beginPiece
(
    const pointField& points,
    const faceList& faces
)
{
    // Basic sizes
    nLocalPoints_ = points.size();
    nLocalFaces_  = faces.size();
    nLocalVerts_  = 0;

    for (const face& f : faces)
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


void Foam::vtk::polyWriter::writePoints
(
    const pointField& points
)
{
    this->beginPoints(numberOfPoints_);

    if (parallel_ ? Pstream::master() : true)
    {
        {
            vtk::writeList(format(), points);
        }
    }

    if (parallel_)
    {
        if (Pstream::master())
        {
            pointField recv;

            // Receive each point field and write
            for (const int subproci : Pstream::subProcs())
            {
                IPstream fromProc(Pstream::commsTypes::blocking, subproci);

                {
                    fromProc >> recv;

                    vtk::writeList(format(), recv);
                }
            }
        }
        else
        {
            // Send
            OPstream toProc
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            {
                toProc << points;
            }
        }
    }


    this->endPoints();
}


void Foam::vtk::polyWriter::writePolysLegacy
(
    const faceList& faces,
    const label pointOffset
)
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
            for (const face& f : faces)
            {
                *iter = f.size();       // The size prefix
                ++iter;

                for (const label pfi : f)
                {
                    *iter = pfi + off;  // Face vertex label
                    ++iter;
                }
            }
            // off += points.size();
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


void Foam::vtk::polyWriter::writePolys
(
    const faceList& faces,
    const label pointOffset
)
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
                for (const face& f : faces)
                {
                    for (const label pfi : f)
                    {
                        *iter = pfi + off;  // Face vertex label
                        ++iter;
                    }
                }
                // off += points.size();
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
            for (const face& f : faces)
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

Foam::vtk::polyWriter::polyWriter
(
    const vtk::outputOptions opts
)
:
    vtk::fileWriter(vtk::fileTag::POLY_DATA, opts),
    numberOfPoints_(0),
    numberOfCells_(0),
    nLocalPoints_(0),
    nLocalFaces_(0),
    nLocalVerts_(0)
{
    // We do not currently support append mode
    opts_.append(false);
}


Foam::vtk::polyWriter::polyWriter
(
    const fileName& file,
    bool parallel
)
:
    // Default parameter fails for gcc-4.8.5, thus specify format here
    polyWriter(vtk::formatType::INLINE_BASE64)
{
    open(file, parallel);
}


Foam::vtk::polyWriter::polyWriter
(
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    polyWriter(opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::polyWriter::writeGeometry()
{
    FatalErrorInFunction
        << "Method was not overloaded, called without a geometry!!" << nl
        << "    Indicates a programming error" << nl << endl
        << abort(FatalError);

    return false;
}


bool Foam::vtk::polyWriter::writePolyGeometry
(
    const pointField& points,
    const faceList& faces
)
{
    enter_Piece();

    beginPiece(points, faces);

    writePoints(points);

    const label pointOffset =
    (
        parallel_ ? globalIndex(nLocalPoints_).localStart() : 0
    );

    if (legacy())
    {
        writePolysLegacy(faces, pointOffset);
    }
    else
    {
        writePolys(faces, pointOffset);
    }

    return true;
}


bool Foam::vtk::polyWriter::beginCellData(label nFields)
{
    return enter_CellData(numberOfCells_, nFields);
}


bool Foam::vtk::polyWriter::beginPointData(label nFields)
{
    return enter_PointData(numberOfPoints_, nFields);
}


// ************************************************************************* //
