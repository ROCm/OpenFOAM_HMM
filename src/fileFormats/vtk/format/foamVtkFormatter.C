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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkFormatter::foamVtkFormatter(std::ostream& os)
:
    os_(os),
    xmlTags_(),
    inTag_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkFormatter::~foamVtkFormatter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

std::size_t Foam::foamVtkFormatter::encodedLength(std::size_t n) const
{
    return n;
}


void Foam::foamVtkFormatter::indent()
{
    label n = xmlTags_.size() * 2;
    while (n--)
    {
        os_ << ' ';
    }
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::xmlHeader()
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


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::comment(const std::string& text)
{
    if (inTag_)
    {
        WarningInFunction
            << "adding xml comment inside a tag??"
            << endl;
    }

    indent();
    os_ << "<!-- " << text << " -->" << nl;

    return *this;
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::openTag(const word& tag)
{
    if (inTag_)
    {
        WarningInFunction
            << "open XML tag '" << tag << "', but already within a tag!"
            << endl;
    }

    indent();
    os_ << '<' << tag;

    xmlTags_.push(tag);
    inTag_ = true;

    return *this;
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::closeTag(const bool isEmpty)
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


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::endTag(const word& tag)
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

    // verify inTag_
    if (!tag.empty() && tag != curr)
    {
        WarningInFunction
            << "expected to end xml-tag '" << tag
            << "' but found '" << curr << "' instead"
            << endl;
    }

    os_  << "</" << curr << '>' << nl;

    inTag_ = false;

    return *this;
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::tag(const word& tag)
{
    openTag(tag);
    closeTag();

    return *this;
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::beginVTKFile
(
    const word& contentType,
    const word& contentVersion,
    const bool leaveOpen
)
{
    openTag("VTKFile");
    xmlAttr("type",        contentType);
    xmlAttr("version",     contentVersion);
    xmlAttr("byte_order",  foamVtkPTraits<Foam::endian>::typeName);
    xmlAttr("header_type", foamVtkPTraits<headerType>::typeName);
    closeTag();

    openTag(contentType);
    if (!leaveOpen)
    {
        closeTag();
    }

    return *this;
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::beginAppendedData()
{
    openTag("AppendedData");
    xmlAttr("encoding", encoding());
    closeTag();
    os_ << '_';

    return *this;
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::xmlAttr
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


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::xmlAttr
(
    const word& k,
    const int32_t v,
    const char quote
)
{
    return xmlAttribute(k, v, quote);
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::xmlAttr
(
    const word& k,
    const int64_t v,
    const char quote
)
{
    return xmlAttribute(k, v, quote);
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::xmlAttr
(
    const word& k,
    const uint64_t v,
    const char quote
)
{
    return xmlAttribute(k, v, quote);
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::xmlAttr
(
    const word& k,
    const scalar v,
    const char quote
)
{
    return xmlAttribute(k, v, quote);
}


// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

Foam::foamVtkFormatter&
Foam::foamVtkFormatter::operator()(const word& k, const std::string& v)
{
    return xmlAttr(k, v);
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::operator()(const word& k, const int32_t v)
{
    return xmlAttr(k, v);
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::operator()(const word& k, const int64_t v)
{
    return xmlAttr(k, v);
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::operator()(const word& k, const uint64_t v)
{
    return xmlAttr(k, v);
}


Foam::foamVtkFormatter&
Foam::foamVtkFormatter::operator()(const word& k, const scalar v)
{
    return xmlAttr(k, v);
}


// ************************************************************************* //
