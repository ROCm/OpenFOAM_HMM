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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkOutput::patchWriter::patchWriter
(
    const fvMesh& mesh,
    const fileName& baseName,
    const foamVtkOutput::outputOptions outOpts,
    const bool nearCellValue,
    const labelList& patchIDs
)
:
    mesh_(mesh),
    format_(),
    nearCellValue_(nearCellValue),
    patchIDs_(patchIDs),
    os_()
{
    outputOptions opts(outOpts);
    opts.legacy(true);  // Legacy only, no append

    os_.open((baseName + (opts.legacy() ? ".vtk" : ".vtp")).c_str());
    format_ = opts.newFormatter(os_);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (opts.legacy())
    {
        foamVtkOutput::legacy::fileHeader
        (
            format(),
            (
                patchIDs_.size() == 1
              ? patches[patchIDs_.first()].name()
              : "patches"
            )
        ) << "DATASET POLYDATA" << nl;
    }

    //------------------------------------------------------------------

    // Write topology
    nPoints_ = 0;
    nFaces_ = 0;
    label nFaceVerts = 0;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        nPoints_ += pp.nPoints();
        nFaces_ += pp.size();

        nFaceVerts += pp.size();
        forAll(pp, facei)
        {
            nFaceVerts += pp[facei].size();
        }
    }

    os_ << "POINTS " << nPoints_ << " float" << nl;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        foamVtkOutput::writeList(format(), pp.localPoints());
    }
    format().flush();

    os_ << "POLYGONS " << nFaces_ << ' ' << nFaceVerts << nl;

    label off = 0;
    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        forAll(pp, facei)
        {
            const face& f = pp.localFaces()[facei];

            format().write(f.size());
            forAll(f, fi)
            {
                format().write(off + f[fi]);
            }
        }

        off += pp.nPoints();
    }
    format().flush();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkOutput::patchWriter::~patchWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamVtkOutput::patchWriter::beginCellData(label nFields)
{
    foamVtkOutput::legacy::cellDataHeader(os(), nFaces_, nFields);
}


void Foam::foamVtkOutput::patchWriter::endCellData()
{}


void Foam::foamVtkOutput::patchWriter::beginPointData(label nFields)
{
    foamVtkOutput::legacy::pointDataHeader(os(), nPoints_, nFields);
}


void Foam::foamVtkOutput::patchWriter::endPointData()
{}


void Foam::foamVtkOutput::patchWriter::writePatchIDs()
{
    os_ << "patchID 1 " << nFaces_ << " float" << nl;

    forAll(patchIDs_, i)
    {
        const label patchId = patchIDs_[i];

        const polyPatch& pp = mesh_.boundaryMesh()[patchId];

        forAll(pp, facei)
        {
            format().write(patchId);
        }
    }
    format().flush();
}


// ************************************************************************* //
