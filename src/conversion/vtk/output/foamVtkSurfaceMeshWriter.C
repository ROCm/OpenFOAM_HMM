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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkOutput::surfaceMeshWriter::surfaceMeshWriter
(
    const bool binary,
    const indirectPrimitivePatch& pp,
    const word& name,
    const fileName& fName
)
:
    pp_(pp),
    format_(),
    os_(fName.c_str())
{
    format_ = foamVtkOutput::newFormatter
    (
        os_,
        (
            binary
          ? foamVtkOutput::LEGACY_BINARY
          : foamVtkOutput::LEGACY_ASCII
        )
    );

    foamVtkOutput::legacy::fileHeader(format(), name)
        << "DATASET POLYDATA" << nl;

    //------------------------------------------------------------------

    // Write topology
    label nFaceVerts = pp.size();

    forAll(pp, facei)
    {
        nFaceVerts += pp[facei].size();
    }

    os_ << "POINTS " << pp.nPoints() << " float" << nl;

    foamVtkOutput::writeList(format(), pp.localPoints());
    format().flush();

    os_ << "POLYGONS " << pp.size() << ' ' << nFaceVerts << nl;

    forAll(pp, facei)
    {
        const face& f = pp.localFaces()[facei];

        format().write(f.size());
        foamVtkOutput::writeList(format(), f);
    }
    format().flush();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamVtkOutput::surfaceMeshWriter::beginCellData(label nFields)
{
    foamVtkOutput::legacy::cellDataHeader(os(), pp_.size(), nFields);
}


void Foam::foamVtkOutput::surfaceMeshWriter::endCellData()
{}


// ************************************************************************* //
