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

const Foam::Enum
<
    Foam::vtk::fileTag
>
Foam::vtk::legacy::contentNames
{
    { vtk::fileTag::POLY_DATA, "POLYDATA" },
    { vtk::fileTag::UNSTRUCTURED_GRID, "UNSTRUCTURED_GRID" },
};


const Foam::Enum
<
    Foam::vtk::fileTag
>
Foam::vtk::legacy::dataTypeNames
{
    { vtk::fileTag::CELL_DATA, "CELL_DATA" },
    { vtk::fileTag::POINT_DATA, "POINT_DATA" }
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::vtk::formatter>
Foam::vtk::newFormatter
(
    std::ostream& os,
    const enum formatType fmtType,
    unsigned prec
)
{
    autoPtr<vtk::formatter> fmt;

    switch (fmtType)
    {
        case formatType::INLINE_ASCII:
            fmt.reset(new vtk::asciiFormatter(os, prec));
            break;

        case formatType::INLINE_BASE64:
            fmt.reset(new vtk::base64Formatter(os));
            break;

        case formatType::APPEND_BASE64:
            fmt.reset(new vtk::appendBase64Formatter(os));
            break;

        case formatType::APPEND_BINARY:
            fmt.reset(new vtk::appendRawFormatter(os));
            break;

        case formatType::LEGACY_ASCII:
            fmt.reset(new vtk::legacyAsciiFormatter(os, prec));
            break;

        case formatType::LEGACY_BINARY:
            fmt.reset(new vtk::legacyRawFormatter(os));
            break;
    }

    return fmt;
}


Foam::label Foam::vtk::writeVtmFile
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
            .openTag(vtk::fileTag::DATA_SET)
            .xmlAttr("index", i)
            .xmlAttr("file", files[i])
            .closeTag(true);
    }

    vtmFile.endTag(content).endVTKFile();

    return files.size();
}


std::ostream& Foam::vtk::legacy::fileHeader
(
    vtk::formatter& fmt,
    const std::string& title,
    const std::string& contentType
)
{
    std::ostream& os = fmt.os();

    fileHeader(os, title, isType<legacyRawFormatter>(fmt));
    if (!contentType.empty())
    {
        os << "DATASET " << contentType.c_str() << nl;
    }

    return os;
}


// ************************************************************************* //
