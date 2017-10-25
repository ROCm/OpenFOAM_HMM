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

#include "foamVtkSurfaceMeshWriter.H"
#include "foamVtkOutput.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::surfaceMeshWriter::beginPiece()
{
    if (!legacy_)
    {
        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, pp_.nPoints())
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POLYS,  pp_.size())
            .closeTag();
    }
}


void Foam::vtk::surfaceMeshWriter::writePoints()
{
    const uint64_t payLoad = (pp_.nPoints()*3*sizeof(float));

    if (legacy_)
    {
        legacy::beginPoints(os_, pp_.nPoints());
    }
    else
    {
        format().tag(vtk::fileTag::POINTS)
            .openDataArray<float, 3>(vtk::dataArrayAttr::POINTS)
            .closeTag();
    }

    format().writeSize(payLoad);

    vtk::writeList(format(), pp_.localPoints());
    format().flush();

    if (!legacy_)
    {
        format()
            .endDataArray()
            .endTag(vtk::fileTag::POINTS);
    }
}


void Foam::vtk::surfaceMeshWriter::writePolysLegacy()
{
    // connectivity count without additional storage (done internally)
    uint64_t nConnectivity = 0;
    forAll(pp_, facei)
    {
        nConnectivity += pp_[facei].size();
    }

    legacy::beginPolys(os_, pp_.size(), nConnectivity);


    // legacy: size + connectivity together
    // [nPts, id1, id2, ..., nPts, id1, id2, ...]

    forAll(pp_, facei)
    {
        const face& f = pp_.localFaces()[facei];

        format().write(f.size());  // The size prefix
        vtk::writeList(format(), f);
    }

    format().flush();
}


void Foam::vtk::surfaceMeshWriter::writePolys()
{
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
        forAll(pp_, facei)
        {
            payLoad += pp_[facei].size();
        }

        format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
            .closeTag();

        // payload size
        format().writeSize(payLoad * sizeof(label));

        forAll(pp_, facei)
        {
            const face& f = pp_.localFaces()[facei];
            vtk::writeList(format(), f);
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
        format().writeSize(pp_.size() * sizeof(label));

        label off = 0;
        forAll(pp_, facei)
        {
            off += pp_[facei].size();

            format().write(off);
        }

        format().flush();
        format().endDataArray();
    }

    format().endTag(vtk::fileTag::POLYS);
}


void Foam::vtk::surfaceMeshWriter::writeMesh()
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

Foam::vtk::surfaceMeshWriter::surfaceMeshWriter
(
    const indirectPrimitivePatch& pp,
    const word& name,
    const fileName& baseName,
    const vtk::outputOptions outOpts
)
:
    pp_(pp),
    legacy_(outOpts.legacy()),
    format_(),
    os_()
{
    outputOptions opts(outOpts);
    opts.append(false);  // No append supported

    os_.open((baseName + (legacy_ ? ".vtk" : ".vtp")).c_str());
    format_ = opts.newFormatter(os_);

    if (legacy_)
    {
        legacy::fileHeader(format(), name, vtk::fileTag::POLY_DATA);
    }
    else
    {
        // XML (inline)

        format()
            .xmlHeader()
            .xmlComment(name)
            .beginVTKFile(vtk::fileTag::POLY_DATA, "0.1");

    }

    beginPiece();
    writeMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::surfaceMeshWriter::~surfaceMeshWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::surfaceMeshWriter::beginCellData(label nFields)
{
    if (legacy_)
    {
        legacy::dataHeader(os(), vtk::fileTag::CELL_DATA, pp_.size(), nFields);
    }
    else
    {
        format().tag(vtk::fileTag::CELL_DATA);
    }
}


void Foam::vtk::surfaceMeshWriter::endCellData()
{
    if (!legacy_)
    {
        format().endTag(vtk::fileTag::CELL_DATA);
    }
}


void Foam::vtk::surfaceMeshWriter::writeFooter()
{
    if (!legacy_)
    {
        // slight cheat. </Piece> too
        format().endTag(vtk::fileTag::PIECE);

        format().endTag(vtk::fileTag::POLY_DATA)
            .endVTKFile();
    }
}


// ************************************************************************* //
