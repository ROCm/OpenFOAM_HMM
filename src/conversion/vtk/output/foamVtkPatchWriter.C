/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "foamVtkPatchWriter.H"
#include "foamVtkOutput.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::patchWriter::beginPiece()
{
    if (!legacy_)
    {
        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, nPoints_)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POLYS,  nFaces_)
            .closeTag();
    }
}


void Foam::vtk::patchWriter::writePoints()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const uint64_t payLoad = (nPoints_*3*sizeof(float));

    if (legacy_)
    {
        legacy::beginPoints(os_, nPoints_);
    }
    else
    {
        format().tag(vtk::fileTag::POINTS)
            .openDataArray<float, 3>(vtk::dataArrayAttr::POINTS)
            .closeTag();
    }

    format().writeSize(payLoad);

    for (const label patchId : patchIDs_)
    {
        const polyPatch& pp = patches[patchId];

        vtk::writeList(format(), pp.localPoints());
    }
    format().flush();

    if (!legacy_)
    {
        format()
            .endDataArray()
            .endTag(vtk::fileTag::POINTS);
    }
}


void Foam::vtk::patchWriter::writePolysLegacy()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // connectivity count without additional storage (done internally)
    uint64_t nConnectivity = 0;
    for (const label patchId : patchIDs_)
    {
        const polyPatch& pp = patches[patchId];

        forAll(pp, facei)
        {
            nConnectivity += pp[facei].size();
        }
    }

    legacy::beginPolys(os_, nFaces_, nConnectivity);


    // legacy: size + connectivity together
    // [nPts, id1, id2, ..., nPts, id1, id2, ...]

    label off = 0;
    for (const label patchId : patchIDs_)
    {
        const polyPatch& pp = patches[patchId];

        forAll(pp, facei)
        {
            const face& f = pp.localFaces()[facei];

            format().write(f.size());  // The size prefix
            forAll(f, fi)
            {
                format().write(off + f[fi]);
            }
        }
        off += pp.nPoints();
    }

    format().flush();
}


void Foam::vtk::patchWriter::writePolys()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    //
    // 'connectivity'
    //

    format().tag(vtk::fileTag::POLYS);

    //
    // 'connectivity'
    //
    {
        // payload count
        uint64_t payLoad = 0;
        for (const label patchId : patchIDs_)
        {
            const polyPatch& pp = patches[patchId];

            forAll(pp, facei)
            {
                payLoad += pp[facei].size();
            }
        }

        format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
            .closeTag();

        // payload size
        format().writeSize(payLoad * sizeof(label));

        label off = 0;
        for (const label patchId : patchIDs_)
        {
            const polyPatch& pp = patches[patchId];

            forAll(pp, facei)
            {
                const face& f = pp.localFaces()[facei];
                forAll(f, fi)
                {
                    format().write(off + f[fi]);
                }
            }

            off += pp.nPoints();
        }

        format().flush();

        format()
            .endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        format()
            .openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
            .closeTag();

        // payload size
        format().writeSize(nFaces_ * sizeof(label));

        label off = 0;
        for (const label patchId : patchIDs_)
        {
            const polyPatch& pp = patches[patchId];

            forAll(pp, facei)
            {
                off += pp[facei].size();

                format().write(off);
            }
        }

        format().flush();
        format().endDataArray();
    }

    format().endTag(vtk::fileTag::POLYS);
}


void Foam::vtk::patchWriter::writeMesh()
{
    writePoints();
    if (legacy_)
    {
        writePolysLegacy();
    }
    else
    {
        writePolys();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::patchWriter::patchWriter
(
    const fvMesh& mesh,
    const fileName& baseName,
    const vtk::outputOptions outOpts,
    const bool nearCellValue,
    const labelList& patchIDs
)
:
    mesh_(mesh),
    legacy_(outOpts.legacy()),
    format_(),
    nearCellValue_(nearCellValue),
    patchIDs_(patchIDs),
    os_(),
    nPoints_(0),
    nFaces_(0)
{
    outputOptions opts(outOpts);
    opts.append(false);  // No append supported

    os_.open((baseName + (legacy_ ? ".vtk" : ".vtp")).c_str());
    format_ = opts.newFormatter(os_);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const word& title =
    (
        patchIDs_.size() == 1
      ? patches[patchIDs_.first()].name()
      : "patches"
    );

    // Basic sizes
    nPoints_ = nFaces_ = 0;
    for (const label patchId : patchIDs_)
    {
        const polyPatch& pp = patches[patchId];

        nPoints_ += pp.nPoints();
        nFaces_  += pp.size();
    }


    if (legacy_)
    {
        legacy::fileHeader(format(), title, vtk::fileTag::POLY_DATA);
    }
    else
    {
        // XML (inline)

        format()
            .xmlHeader()
            .xmlComment(title)
            .beginVTKFile(vtk::fileTag::POLY_DATA, "0.1");
    }

    beginPiece();
    writeMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::patchWriter::~patchWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::patchWriter::beginCellData(label nFields)
{
    if (legacy_)
    {
        legacy::dataHeader
        (
            os(),
            vtk::fileTag::CELL_DATA,
            nFaces_,
            nFields
        );
    }
    else
    {
        format().tag(vtk::fileTag::CELL_DATA);
    }
}


void Foam::vtk::patchWriter::endCellData()
{
    if (!legacy_)
    {
        format().endTag(vtk::fileTag::CELL_DATA);
    }
}


void Foam::vtk::patchWriter::beginPointData(label nFields)
{
    if (legacy_)
    {
        legacy::dataHeader
        (
            os(),
            vtk::fileTag::POINT_DATA,
            nPoints_,
            nFields
        );
    }
    else
    {
        format().tag(vtk::fileTag::POINT_DATA);
    }
}


void Foam::vtk::patchWriter::endPointData()
{
    if (!legacy_)
    {
        format().endTag(vtk::fileTag::POINT_DATA);
    }
}


void Foam::vtk::patchWriter::writeFooter()
{
    if (!legacy_)
    {
        // slight cheat. </Piece> too
        format().endTag(vtk::fileTag::PIECE);

        format().endTag(vtk::fileTag::POLY_DATA)
            .endVTKFile();
    }
}


void Foam::vtk::patchWriter::writePatchIDs()
{
    // Patch ids first
    const uint64_t payLoad = nFaces_ * sizeof(label);

    if (legacy_)
    {
        legacy::intField(os_, "patchID", 1, nFaces_);
    }
    else
    {
        format().openDataArray<label>("patchID")
            .closeTag();
    }

    format().writeSize(payLoad);

    for (const label patchId : patchIDs_)
    {
        const label sz = mesh_.boundaryMesh()[patchId].size();

        for (label facei = 0; facei < sz; ++facei)
        {
            format().write(patchId);
        }
    }
    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


// ************************************************************************* //
