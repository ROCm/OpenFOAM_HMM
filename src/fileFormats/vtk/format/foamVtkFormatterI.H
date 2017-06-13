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

inline Foam::vtk::formatter::formatter(std::ostream& os)
:
    os_(os),
    xmlTags_(),
    inTag_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline std::ostream& Foam::vtk::formatter::os()
{
    return os_;
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::openTag(const vtk::fileTag& tagEnum)
{
    return openTag(vtk::fileTagNames[tagEnum]);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::endTag(const vtk::fileTag& tagEnum)
{
    return endTag(vtk::fileTagNames[tagEnum]);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::tag(const word& tagName)
{
    openTag(tagName);
    closeTag();

    return *this;
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::tag(const vtk::fileTag& tagEnum)
{
    return tag(vtk::fileTagNames[tagEnum]);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::beginVTKFile
(
    const vtk::fileTag& contentType,
    const word& contentVersion,
    const bool leaveOpen
)
{
    return beginVTKFile
    (
        vtk::fileTagNames[contentType],
        contentVersion,
        leaveOpen
    );
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::endDataArray()
{
    return endTag("DataArray");
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const int32_t v,
    const char quote
)
{
    return writeAttribute(k, v, quote);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const int64_t v,
    const char quote
)
{
    return writeAttribute(k, v, quote);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const uint64_t v,
    const char quote
)
{
    return writeAttribute(k, v, quote);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const scalar v,
    const char quote
)
{
    return writeAttribute(k, v, quote);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::xmlAttr
(
    const vtk::fileAttr& attrEnum,
    const int32_t v,
    const char quote
)
{
    return xmlAttr(vtk::fileAttrNames[attrEnum], v, quote);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::xmlAttr
(
    const vtk::fileAttr& attrEnum,
    const int64_t v,
    const char quote
)
{
    return xmlAttr(vtk::fileAttrNames[attrEnum], v, quote);
}


inline Foam::vtk::formatter&
Foam::vtk::formatter::xmlAttr
(
    const vtk::fileAttr& attrEnum,
    const uint64_t v,
    const char quote
)
{
    return xmlAttr(vtk::fileAttrNames[attrEnum], v, quote);
}


// ************************************************************************* //
