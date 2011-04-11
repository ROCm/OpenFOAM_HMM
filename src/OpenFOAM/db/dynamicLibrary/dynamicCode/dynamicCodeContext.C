/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "OSHA1stream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCodeContext::dynamicCodeContext(const dictionary& dict)
:
    dict_(dict),
    code_(),
    localCode_(),
    include_(),
    options_(),
    libs_()
{
    // expand dictionary entries

    {
        const entry& codeEntry = dict.lookupEntry("code", false, false);
        code_ = stringOps::trim(codeEntry.stream());
        stringOps::inplaceExpand(code_, dict);
        addLineDirective(code_, codeEntry.startLineNumber(), dict.name());
    }

    // note: removes any leading/trailing whitespace
    // - necessary for compilation options, convenient for includes
    // and body.

    // optional
    const entry* includePtr = dict.lookupEntryPtr
    (
        "codeInclude",
        false,
        false
    );
    if (includePtr)
    {
        include_ = stringOps::trim(includePtr->stream());
        stringOps::inplaceExpand(include_, dict);
        addLineDirective(include_, includePtr->startLineNumber(), dict.name());
    }

    // optional
    const entry* optionsPtr = dict.lookupEntryPtr
    (
        "codeOptions",
        false,
        false
    );
    if (optionsPtr)
    {
        options_ = stringOps::trim(optionsPtr->stream());
        stringOps::inplaceExpand(options_, dict);
    }

    // optional
    const entry* libsPtr = dict.lookupEntryPtr("codeLibs", false, false);
    if (libsPtr)
    {
        libs_ = stringOps::trim(libsPtr->stream());
        stringOps::inplaceExpand(libs_, dict);
    }

    // calculate SHA1 digest from include, options, localCode, code
    OSHA1stream os;
    os  << include_ << options_ << libs_ << localCode_ << code_;
    sha1_ = os.digest();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::dynamicCodeContext::addLineDirective
(
    string& code,
    const label lineNum,
    const fileName& name
)
{
    code = "#line " + Foam::name(lineNum) + " \"" + name + "\"\n" + code;
}


// ************************************************************************* //
