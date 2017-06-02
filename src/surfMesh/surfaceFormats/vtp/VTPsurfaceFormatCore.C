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

#include "VTPsurfaceFormatCore.H"
#include "clock.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::VTPsurfaceFormatCore::writeHeader
(
    foamVtkOutput::formatter& format,
    const pointField& pts,
    const label nFaces
)
{
    // XML (inline)

    format
        .xmlHeader()
        .xmlComment("surface written " + clock::dateTime())
        .beginVTKFile(vtkFileTag::POLY_DATA, "0.1");

    // <Piece>
    format
        .openTag(vtkFileTag::PIECE)
        ( "NumberOfPoints", pts.size() )
        ( "NumberOfPolys",  nFaces )
        .closeTag();


    // Points

    const uint64_t payLoad = (pts.size()*3* sizeof(float));

    format.tag(vtkFileTag::POINTS)
        .openDataArray<float, 3>(vtkFileTag::POINTS)
        .closeTag();

    format.writeSize(payLoad);
    foamVtkOutput::writeList(format, pts);
    format.flush();

    format
        .endDataArray()
        .endTag(vtkFileTag::POINTS);
}


void Foam::fileFormats::VTPsurfaceFormatCore::writeFooter
(
    foamVtkOutput::formatter& format
)
{
    // Slight cheat. </Piece> too
    format.endTag(Foam::vtkFileTag::PIECE);

    format.endTag(vtkFileTag::POLY_DATA)
        .endVTKFile();
}



void Foam::fileFormats::VTPsurfaceFormatCore::writeCellData
(
    foamVtkOutput::formatter& format,
    const UList<surfZone>& zones
)
{
    // Zone ids as CellData

    // Number of faces covered by the zones
    uint64_t payLoad = 0;
    for (const auto& z : zones)
    {
        payLoad += z.size();
    }

    format.tag(vtkFileTag::CELL_DATA);
    format.openDataArray<label>("region")
        .closeTag();

    format.writeSize(payLoad * sizeof(label));

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
    format.endDataArray();

    format.endTag(vtkFileTag::CELL_DATA);
}


void Foam::fileFormats::VTPsurfaceFormatCore::writeCellData
(
    foamVtkOutput::formatter& format,
    const labelUList& zoneIds
)
{
    // Zone ids as CellData

    format.tag(vtkFileTag::CELL_DATA);
    format.openDataArray<label>("region")
        .closeTag();

    const uint64_t payLoad(zoneIds.size() * sizeof(label));

    format.writeSize(payLoad);
    foamVtkOutput::writeList(format, zoneIds);

    format.flush();
    format.endDataArray();

    format.endTag(vtkFileTag::CELL_DATA);

}


// ************************************************************************* //
