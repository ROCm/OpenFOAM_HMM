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

void Foam::foamVtkOutput::surfaceMeshWriter::beginPiece()
{
    if (!legacy_)
    {
        format()
            .openTag(vtkFileTag::PIECE)
            ( "NumberOfPoints", pp_.nPoints() )
            ( "NumberOfPolys",  pp_.size() )
            .closeTag();
    }
}


void Foam::foamVtkOutput::surfaceMeshWriter::writePoints()
{
    // payload count
    const uint64_t payLoad = (pp_.nPoints()*3*sizeof(float));

    if (legacy_)
    {
        legacy::beginPoints(os_, pp_.nPoints());
    }
    else
    {
        format().tag(vtkFileTag::POINTS)
            .openDataArray<float, 3>(vtkFileTag::POINTS)
            .closeTag();
    }

    format().writeSize(payLoad);

    foamVtkOutput::writeList(format(), pp_.localPoints());
    format().flush();

    if (!legacy_)
    {
        format()
            .endDataArray()
            .endTag(vtkFileTag::POINTS);
    }
}


void Foam::foamVtkOutput::surfaceMeshWriter::writePolysLegacy()
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
        foamVtkOutput::writeList(format(), f);
    }

    format().flush();
}


void Foam::foamVtkOutput::surfaceMeshWriter::writePolys()
{
    //
    // 'connectivity'
    //

    format().tag(vtkFileTag::POLYS);

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

        format().openDataArray<label>("connectivity")
            .closeTag();

        // payload size
        format().writeSize(payLoad * sizeof(label));

        forAll(pp_, facei)
        {
            const face& f = pp_.localFaces()[facei];
            foamVtkOutput::writeList(format(), f);
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
            .openDataArray<label>("offsets")
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

    format().endTag(vtkFileTag::POLYS);
}


void Foam::foamVtkOutput::surfaceMeshWriter::writeMesh()
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

Foam::foamVtkOutput::surfaceMeshWriter::surfaceMeshWriter
(
    const indirectPrimitivePatch& pp,
    const word& name,
    const fileName& baseName,
    const foamVtkOutput::outputOptions outOpts
)
:
    pp_(pp),
    legacy_(outOpts.legacy()),
    format_(),
    os_()
{
    outputOptions opts(outOpts);
    opts.legacy(true);  // No append supported

    os_.open((baseName + (legacy_ ? ".vtk" : ".vtp")).c_str());
    format_ = opts.newFormatter(os_);

    if (legacy_)
    {
        legacy::fileHeader(format(), name, vtkFileTag::POLY_DATA);
    }
    else
    {
        // XML (inline)

        format()
            .xmlHeader()
            .xmlComment(name)
            .beginVTKFile(vtkFileTag::POLY_DATA, "0.1");

    }

    beginPiece();
    writeMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkOutput::surfaceMeshWriter::~surfaceMeshWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamVtkOutput::surfaceMeshWriter::beginCellData(label nFields)
{
    if (legacy_)
    {
        legacy::dataHeader(os(), vtkFileTag::CELL_DATA, pp_.size(), nFields);
    }
    else
    {
        format().tag(vtkFileTag::CELL_DATA);
    }
}


void Foam::foamVtkOutput::surfaceMeshWriter::endCellData()
{
    if (!legacy_)
    {
        format().endTag(vtkFileTag::CELL_DATA);
    }
}


void Foam::foamVtkOutput::surfaceMeshWriter::writeFooter()
{
    if (!legacy_)
    {
        // slight cheat. </Piece> too
        format().endTag(vtkFileTag::PIECE);

        format().endTag(vtkFileTag::POLY_DATA)
            .endVTKFile();
    }
}


// ************************************************************************* //
