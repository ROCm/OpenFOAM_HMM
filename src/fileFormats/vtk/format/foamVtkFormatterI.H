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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::foamVtkOutput::formatter::formatter(std::ostream& os)
:
    os_(os),
    xmlTags_(),
    inTag_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline std::ostream& Foam::foamVtkOutput::formatter::os()
{
    return os_;
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::openTag(const vtkFileTag& tagEnum)
{
    return openTag(vtkFileTagNames[tagEnum]);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::endTag(const vtkFileTag& tagEnum)
{
    return endTag(vtkFileTagNames[tagEnum]);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::tag(const word& tagName)
{
    openTag(tagName);
    closeTag();

    return *this;
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::tag(const vtkFileTag& tagEnum)
{
    return tag(vtkFileTagNames[tagEnum]);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::beginVTKFile
(
    const vtkFileTag& contentType,
    const word& contentVersion,
    const bool leaveOpen
)
{
    return beginVTKFile
    (
        vtkFileTagNames[contentType],
        contentVersion,
        leaveOpen
    );
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::endDataArray()
{
    return endTag("DataArray");
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::endAppendedData()
{
    flush();     // flush any pending encoded content
    os_ << '\n'; // clear separation from content.
    return endTag("AppendedData");
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::endVTKFile()
{
    return endTag(vtkFileTag::VTK_FILE);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::xmlAttr
(
    const word& k,
    const int32_t v,
    const char quote
)
{
    return xmlAttribute(k, v, quote);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::xmlAttr
(
    const word& k,
    const int64_t v,
    const char quote
)
{
    return xmlAttribute(k, v, quote);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::xmlAttr
(
    const word& k,
    const uint64_t v,
    const char quote
)
{
    return xmlAttribute(k, v, quote);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::xmlAttr
(
    const word& k,
    const scalar v,
    const char quote
)
{
    return xmlAttribute(k, v, quote);
}


// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::operator()(const word& k, const std::string& v)
{
    return xmlAttr(k, v);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::operator()(const word& k, const int32_t v)
{
    return xmlAttr(k, v);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::operator()(const word& k, const int64_t v)
{
    return xmlAttr(k, v);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::operator()(const word& k, const uint64_t v)
{
    return xmlAttr(k, v);
}


inline Foam::foamVtkOutput::formatter&
Foam::foamVtkOutput::formatter::operator()(const word& k, const scalar v)
{
    return xmlAttr(k, v);
}


// ************************************************************************* //
