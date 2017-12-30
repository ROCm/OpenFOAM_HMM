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

#include "foamVtkLagrangianWriter.H"
#include "Cloud.H"
#include "passiveParticle.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::lagrangianWriter::beginPiece()
{
    if (!legacy_)
    {
        if (useVerts_)
        {
            format()
                .openTag(vtk::fileTag::PIECE)
                .xmlAttr(fileAttr::NUMBER_OF_POINTS, nParcels_)
                .xmlAttr(fileAttr::NUMBER_OF_VERTS,  nParcels_)
                .closeTag();
        }
        else
        {
            format()
                .openTag(vtk::fileTag::PIECE)
                .xmlAttr(fileAttr::NUMBER_OF_POINTS, nParcels_)
                .closeTag();
        }
    }
}


void Foam::vtk::lagrangianWriter::writePoints()
{
    Cloud<passiveParticle> parcels(mesh_, cloudName_, false);
    nParcels_ = parcels.size();

    const uint64_t payLoad = (nParcels_ * 3 * sizeof(float));

    if (legacy_)
    {
        legacy::beginPoints(os_, nParcels_);
    }
    else
    {
        beginPiece(); // Tricky - hide in here

        format().tag(vtk::fileTag::POINTS)
            .openDataArray<float,3>(vtk::dataArrayAttr::POINTS)
            .closeTag();
    }

    format().writeSize(payLoad);

    forAllConstIters(parcels, iter)
    {
        const point pt(iter().position());

        vtk::write(format(), pt);
    }
    format().flush();

    if (!legacy_)
    {
        format()
            .endDataArray()
            .endTag(vtk::fileTag::POINTS);
    }
}


void Foam::vtk::lagrangianWriter::writeVertsLegacy()
{
    os_ << "VERTICES " << nParcels_ << ' ' << 2*nParcels_ << nl;

    // legacy has cells + connectivity together
    // count the number of vertices referenced

    for (label i=0; i < nParcels_; ++i)
    {
        format().write(label(1)); // Number of vertices for this cell (==1)
        format().write(i);
    }
    format().flush();
}


void Foam::vtk::lagrangianWriter::writeVerts()
{
    format().tag(vtk::fileTag::VERTS);

    // Same payload throughout
    const uint64_t payLoad = (nParcels_ * sizeof(label));

    //
    // 'connectivity'
    // = linear mapping onto points
    //
    {
        format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
            .closeTag();

        format().writeSize(payLoad);
        for (label i=0; i < nParcels_; ++i)
        {
            format().write(i);
        }
        format().flush();

        format().endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    // = linear mapping onto points (with 1 offset)
    //
    {
        format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
            .closeTag();

        format().writeSize(payLoad);
        for (label i=0; i < nParcels_; ++i)
        {
            format().write(i+1);
        }
        format().flush();

        format().endDataArray();
    }

    format().endTag(vtk::fileTag::VERTS);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::lagrangianWriter::lagrangianWriter
(
    const fvMesh& mesh,
    const word& cloudName,
    const fileName& baseName,
    const vtk::outputOptions outOpts,
    const bool dummyCloud
)
:
    mesh_(mesh),
    legacy_(outOpts.legacy()),
    useVerts_(false),
    format_(),
    cloudName_(cloudName),
    os_(),
    nParcels_(0)
{

    outputOptions opts(outOpts);
    opts.append(false);  // No append supported

    os_.open((baseName + (legacy_ ? ".vtk" : ".vtp")).c_str());
    format_ = opts.newFormatter(os_);

    const auto& title = mesh_.time().caseName();

    if (legacy_)
    {
        legacy::fileHeader(format(), title, vtk::fileTag::POLY_DATA);

        if (dummyCloud)
        {
            legacy::beginPoints(os_, nParcels_);
        }
        else
        {
            writePoints();
            if (useVerts_) writeVertsLegacy();
        }
    }
    else
    {
        // XML (inline)

        format()
            .xmlHeader()
            .xmlComment(title)
            .beginVTKFile(vtk::fileTag::POLY_DATA, "0.1");

        if (dummyCloud)
        {
            beginPiece();
        }
        else
        {
            writePoints();
            if (useVerts_) writeVerts();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::lagrangianWriter::~lagrangianWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::lagrangianWriter::beginParcelData(label nFields)
{
    const vtk::fileTag dataType =
    (
        useVerts_
      ? vtk::fileTag::CELL_DATA
      : vtk::fileTag::POINT_DATA
    );

    if (legacy_)
    {
        legacy::dataHeader(os_, dataType, nParcels_, nFields);
    }
    else
    {
        format().tag(dataType);
    }
}


void Foam::vtk::lagrangianWriter::endParcelData()
{
    const vtk::fileTag dataType =
    (
        useVerts_
      ? vtk::fileTag::CELL_DATA
      : vtk::fileTag::POINT_DATA
    );

    if (!legacy_)
    {
        format().endTag(dataType);
    }
}


void Foam::vtk::lagrangianWriter::writeFooter()
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
