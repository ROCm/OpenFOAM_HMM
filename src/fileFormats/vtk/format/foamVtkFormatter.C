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

#include "foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::formatter::~formatter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

std::size_t Foam::vtk::formatter::encodedLength(std::size_t n) const
{
    return n;
}


void Foam::vtk::formatter::indent()
{
    label n = xmlTags_.size() * 2;
    while (n--)
    {
        os_ << ' ';
    }
}


Foam::vtk::formatter&
Foam::vtk::formatter::xmlHeader()
{
    if (inTag_)
    {
        WarningInFunction
            << "xml header, but already within a tag!"
            << endl;
    }

    os_ << "<?xml version='1.0'?>" << nl;

    return *this;
}


Foam::vtk::formatter&
Foam::vtk::formatter::xmlComment(const std::string& comment)
{
    if (inTag_)
    {
        WarningInFunction
            << "adding xml comment inside a tag??"
            << endl;
    }

    indent();
    os_ << "<!-- " << comment << " -->" << nl;

    return *this;
}


Foam::vtk::formatter&
Foam::vtk::formatter::openTag(const word& tagName)
{
    if (inTag_)
    {
        WarningInFunction
            << "open XML tag '" << tagName
            << "', but already within a tag!"
            << endl;
    }

    indent();
    os_ << '<' << tagName;

    xmlTags_.push(tagName);
    inTag_ = true;

    return *this;
}


Foam::vtk::formatter&
Foam::vtk::formatter::closeTag(const bool isEmpty)
{
    if (!inTag_)
    {
        WarningInFunction
            << "close XML tag, but not within a tag!"
            << endl;
    }

    if (isEmpty)
    {
        // eg, <tag ... />
        xmlTags_.pop();
        os_ << " /";
    }
    os_ << '>' << nl;

    inTag_ = false;

    return *this;
}


Foam::vtk::formatter&
Foam::vtk::formatter::endTag(const word& tagName)
{
    const word curr = xmlTags_.pop();
    indent();

    if (inTag_)
    {
        WarningInFunction
            << "adding XML endTag '" << curr
            << "' but already in another tag!"
            << endl;
    }

    // verify expected end tag
    if (!tagName.empty() && tagName != curr)
    {
        WarningInFunction
            << "expecting to end xml-tag '" << tagName
            << "' but found '" << curr << "' instead"
            << endl;
    }

    os_  << "</" << curr << '>' << nl;

    inTag_ = false;

    return *this;
}


Foam::vtk::formatter&
Foam::vtk::formatter::beginVTKFile
(
    const word& contentType,
    const word& contentVersion,
    const bool leaveOpen
)
{
    openTag(vtk::fileTag::VTK_FILE);
    xmlAttr("type",        contentType);
    xmlAttr("version",     contentVersion);
    xmlAttr("byte_order",  vtkPTraits<Foam::endian>::typeName);
    xmlAttr("header_type", vtkPTraits<headerType>::typeName);
    closeTag();

    openTag(contentType);
    if (!leaveOpen)
    {
        closeTag();
    }

    return *this;
}


Foam::vtk::formatter&
Foam::vtk::formatter::endVTKFile()
{
    return endTag(vtk::fileTag::VTK_FILE);
}


Foam::vtk::formatter&
Foam::vtk::formatter::beginAppendedData()
{
    openTag("AppendedData");
    xmlAttr("encoding", encoding());
    closeTag();
    os_ << '_';

    return *this;
}


Foam::vtk::formatter&
Foam::vtk::formatter::endAppendedData()
{
    flush();     // flush any pending encoded content
    os_ << nl;   // ensure clear separation from content.
    return endTag("AppendedData");
}


Foam::vtk::formatter&
Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const std::string& v,
    const char quote
)
{
    if (!inTag_)
    {
        WarningInFunction
            << "xml attribute '" << k << "' but not within a tag!"
            << endl;
    }

    os_ << ' ' << k << '=' << quote << v.c_str() << quote;

    return *this;
}


// ************************************************************************* //
