/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "foamVtkOutput.H"
#include "foamVtkAsciiFormatter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::foamVtkOutput::legacy::EXT = "vtk";


//! \cond fileScope
static inline std::ostream& legacyDataHeader
(
    std::ostream& os,
    const char* tag,
    const Foam::label nItems,
    const Foam::label nFields
)
{
    os  << tag << ' ' << nItems << '\n'
        << "FIELD attributes " << nFields << '\n';

    return os;
}
//! \endcond


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::foamVtkOutput::writeVtmFile
(
    std::ostream& os,
    const UList<fileName>& files
)
{
    const word& content = "vtkMultiBlockDataSet";

    foamVtkAsciiFormatter vtmFile(os);

    vtmFile
        .xmlHeader()
        .openTag("VTKFile")
        ( "type",        content )
        ( "version",     "1.0" )
        ( "byte_order",  foamVtkFormatter::byteOrder )
        ( "header_type", foamVtkFormatter::headerType )
        .closeTag();

    vtmFile.tag(content);

    forAll(files, i)
    {
        vtmFile
            .openTag("DataSet")
            ( "index", i )
            ( "file", files[i] )
            .closeTag(true);
    }

    vtmFile.endTag(content).endTag("VTKFile");

    return files.size();
}


std::ostream& Foam::foamVtkOutput::legacy::writeHeader
(
    std::ostream& os,
    const std::string& title,
    const bool binary
)
{
    os  << "# vtk DataFile Version 2.0" << nl
        << title << nl
        << (binary ? "BINARY" : "ASCII") << nl;

    return os;
}


std::ostream& Foam::foamVtkOutput::legacy::writeCellDataHeader
(
    std::ostream& os,
    const label nCells,
    const label nFields
)
{
    return legacyDataHeader(os, "CELL_DATA", nCells, nFields);
}


std::ostream& Foam::foamVtkOutput::legacy::writePointDataHeader
(
    std::ostream& os,
    const label nPoints,
    const label nFields
)
{
    return legacyDataHeader(os, "POINT_DATA", nPoints, nFields);
}


// ************************************************************************* //
