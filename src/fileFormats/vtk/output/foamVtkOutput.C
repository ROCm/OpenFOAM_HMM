/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "foamVtkFormatter.H"
#include "foamVtkAsciiFormatter.H"
#include "foamVtkBase64Formatter.H"
#include "foamVtkAppendBase64Formatter.H"
#include "foamVtkAppendRawFormatter.H"
#include "foamVtkLegacyAsciiFormatter.H"
#include "foamVtkLegacyRawFormatter.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * * //

const Foam::Enum<Foam::foamVtkOutput::legacy::textElemType>
Foam::foamVtkOutput::legacy::textNames
{
    { textElemType::POINTS,            "POINTS" },
    { textElemType::CELLS,             "CELLS" },
    { textElemType::POLYS,             "POLYGONS" },
    { textElemType::VERTS,             "VERTICES" },
    { textElemType::POLY_DATA,         "POLYDATA" },
    { textElemType::UNSTRUCTURED_GRID, "UNSTRUCTURED_GRID" },
};


const Foam::Enum<Foam::vtkFileTag>
Foam::foamVtkOutput::legacy::dataTypeNames
{
    { vtkFileTag::CELL_DATA,  "CELL_DATA" },
    { vtkFileTag::POINT_DATA, "POINT_DATA" }
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::foamVtkOutput::formatter>
Foam::foamVtkOutput::newFormatter
(
    std::ostream& os,
    const enum formatType fmtType,
    unsigned prec
)
{
    autoPtr<foamVtkOutput::formatter> fmt;

    switch (fmtType)
    {
        case formatType::INLINE_ASCII:
            fmt.set(new foamVtkOutput::asciiFormatter(os, prec));
            break;

        case formatType::INLINE_BASE64:
            fmt.set(new foamVtkOutput::base64Formatter(os));
            break;

        case formatType::APPEND_BASE64:
            fmt.set(new foamVtkOutput::appendBase64Formatter(os));
            break;

        case formatType::APPEND_BINARY:
            fmt.set(new foamVtkOutput::appendRawFormatter(os));
            break;

        case formatType::LEGACY_ASCII:
            fmt.set(new foamVtkOutput::legacyAsciiFormatter(os, prec));
            break;

        case formatType::LEGACY_BINARY:
            fmt.set(new foamVtkOutput::legacyRawFormatter(os));
            break;
    }

    return fmt;
}


Foam::label Foam::foamVtkOutput::writeVtmFile
(
    std::ostream& os,
    const UList<fileName>& files
)
{
    const word& content = "vtkMultiBlockDataSet";

    asciiFormatter vtmFile(os);

    vtmFile
        .xmlHeader()
        .beginVTKFile(content, "1.0");

    forAll(files, i)
    {
        vtmFile
            .openTag(vtkFileTag::DATA_SET)
            ( "index", i )
            ( "file", files[i] )
            .closeTag(true);
    }

    vtmFile.endTag(content).endVTKFile();

    return files.size();
}


std::ostream& Foam::foamVtkOutput::legacy::fileHeader
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


std::ostream& Foam::foamVtkOutput::legacy::fileHeader
(
    foamVtkOutput::formatter& fmt,
    const std::string& title
)
{
    return fileHeader(fmt.os(), title, isType<legacyRawFormatter>(fmt));
}


std::ostream& Foam::foamVtkOutput::legacy::dataHeader
(
    std::ostream& os,
    const vtkFileTag& dataTypeTag,
    const label nEntries,
    const label nFields
)
{
    os  << dataTypeNames[dataTypeTag] << ' ' << nEntries << nl
        << "FIELD attributes " << nFields << nl;

    return os;
}


std::ostream& Foam::foamVtkOutput::legacy::cellDataHeader
(
    std::ostream& os,
    const label nEntries,
    const label nFields
)
{
    return dataHeader(os, vtkFileTag::CELL_DATA, nEntries, nFields);
}


std::ostream& Foam::foamVtkOutput::legacy::pointDataHeader
(
    std::ostream& os,
    const label nEntries,
    const label nFields
)
{
    return dataHeader(os, vtkFileTag::POINT_DATA, nEntries, nFields);
}


// ************************************************************************* //
