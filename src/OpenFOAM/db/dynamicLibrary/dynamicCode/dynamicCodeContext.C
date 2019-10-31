/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "OSHA1stream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::dynamicCodeContext::inplaceExpand
(
    string& code,
    const dictionary& dict
)
{
    stringOps::inplaceTrim(code);
    stringOps::inplaceExpand(code, dict);
}


unsigned Foam::dynamicCodeContext::addLineDirective
(
    string& code,
    label lineNum,
    const fileName& file
)
{
    ++lineNum;  // Change from 0-based to 1-based

    const auto len = code.length();

    if (lineNum > 0 && len && !file.empty())
    {
        code = "#line " + Foam::name(lineNum) + " \"" + file + "\"\n" + code;

        return (code.length() - len);
    }

    return 0;
}


unsigned Foam::dynamicCodeContext::addLineDirective
(
    string& code,
    label lineNum,
    const dictionary& dict
)
{
    return addLineDirective(code, lineNum, dict.name());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCodeContext::dynamicCodeContext()
:
    dict_(std::cref<dictionary>(dictionary::null))
{}


Foam::dynamicCodeContext::dynamicCodeContext(const dictionary& dict)
:
    dynamicCodeContext()
{
    setCodeContext(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::dynamicCodeContext::valid() const
{
    return &(dict_.get()) != &(dictionary::null);
}


void Foam::dynamicCodeContext::setCodeContext(const dictionary& dict)
{
    dict_ = std::cref<dictionary>(dict);
    sha1_.clear();

    // Expand dictionary entries.
    // Removing any leading/trailing whitespace is necessary for compilation
    // options, but is also convenient for includes and code body.

    const entry* eptr;

    options_.clear();
    sha1_.append("<codeOptions>");
    if ((eptr = dict.findEntry("codeOptions", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(options_);
        dynamicCodeContext::inplaceExpand(options_, dict);
        sha1_.append(options_);
        // No #line for options (Make/options)
    }

    libs_.clear();
    sha1_.append("<codeLibs>");
    if ((eptr = dict.findEntry("codeLibs", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(libs_);
        dynamicCodeContext::inplaceExpand(libs_, dict);
        sha1_.append(libs_);
        // No #line for libs (LIB_LIBS)
    }

    include_.clear();
    sha1_.append("<codeInclude>");
    if ((eptr = dict.findEntry("codeInclude", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(include_);
        dynamicCodeContext::inplaceExpand(include_, dict);
        sha1_.append(include_);
        addLineDirective(include_, eptr->startLineNumber(), dict);
    }

    code_.clear();
    sha1_.append("<code>");
    if ((eptr = dict.findEntry("code", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(code_);
        dynamicCodeContext::inplaceExpand(code_, dict);
        sha1_.append(code_);
        addLineDirective(code_, eptr->startLineNumber(), dict);
    }

    localCode_.clear();
    sha1_.append("<localCode>");
    if ((eptr = dict.findEntry("localCode", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(localCode_);
        dynamicCodeContext::inplaceExpand(localCode_, dict);
        sha1_.append(localCode_);
        addLineDirective(localCode_, eptr->startLineNumber(), dict);
    }
}


// ************************************************************************* //
