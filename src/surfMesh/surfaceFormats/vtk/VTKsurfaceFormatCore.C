/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "VTKsurfaceFormatCore.H"
#include "clock.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::VTKsurfaceFormatCore::writeHeader
(
    vtk::formatter& format,
    const UList<point>& pts
)
{
    vtk::legacy::fileHeader
    (
        format,
        ("surface written " + clock::dateTime()),
        vtk::fileTag::POLY_DATA
    );

    vtk::legacy::beginPoints(format.os(), pts.size());

    vtk::writeList(format, pts);
    format.flush();
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeCellData
(
    vtk::formatter& format,
    const UList<surfZone>& zones
)
{
    // Zone ids as CellData

    // Number of faces covered by the zones
    label nFaces = 0;
    for (const auto& z : zones)
    {
        nFaces += z.size();
    }

    vtk::legacy::dataHeader
    (
        format.os(),
        vtk::fileTag::CELL_DATA,
        nFaces,
        1  // Only one field
    );

    vtk::legacy::intField
    (
        format.os(),
        "region",
        1, // nComponent
        nFaces
    );

    label zoneId = 0;
    for (const surfZone& zone : zones)
    {
        forAll(zone, i)
        {
            format.write(zoneId);
        }
        ++zoneId;
    }
    format.flush();
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeCellData
(
    vtk::formatter& format,
    const labelUList& zoneIds
)
{
    // Zone ids as CellData

    // Number of faces
    const label nFaces = zoneIds.size();

    vtk::legacy::dataHeader
    (
        format.os(),
        vtk::fileTag::CELL_DATA,
        nFaces,
        1  // Only one field
    );

    vtk::legacy::intField
    (
        format.os(),
        "region",
        1, // nComponent
        nFaces
    );

    vtk::writeList(format, zoneIds);
    format.flush();
}


// ************************************************************************* //
