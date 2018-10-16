/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    // Expand dictionary entries

    // Note: removes any leading/trailing whitespace
    // - necessary for compilation options, convenient for includes
    // and body.

    const entry* codePtr = dict.findEntry("code", keyType::LITERAL);

    if (codePtr)
    {
        codePtr->readEntry(code_);
        stringOps::inplaceTrim(code_);
        stringOps::inplaceExpand(code_, dict);
    }

    const entry* includePtr = dict.findEntry("codeInclude", keyType::LITERAL);

    if (includePtr)
    {
        includePtr->readEntry(include_);
        stringOps::inplaceTrim(include_);
        stringOps::inplaceExpand(include_, dict);
    }

    const entry* optionsPtr = dict.findEntry("codeOptions", keyType::LITERAL);

    if (optionsPtr)
    {
        optionsPtr->readEntry(options_);
        stringOps::inplaceTrim(options_);
        stringOps::inplaceExpand(options_, dict);
    }

    const entry* libsPtr = dict.findEntry("codeLibs", keyType::LITERAL);

    if (libsPtr)
    {
        libsPtr->readEntry(libs_);
        stringOps::inplaceTrim(libs_);
        stringOps::inplaceExpand(libs_, dict);
    }

    const entry* localPtr = dict.findEntry("localCode", keyType::LITERAL);

    if (localPtr)
    {
        localPtr->readEntry(localCode_);
        stringOps::inplaceTrim(localCode_);
        stringOps::inplaceExpand(localCode_, dict);
    }

    // Calculate SHA1 digest from include, options, localCode, code
    OSHA1stream os;
    os  << include_ << options_ << libs_ << localCode_ << code_;
    sha1_ = os.digest();


    // Add line number after calculating sha1 since includes processorDDD
    // in path which differs between processors.

    if (codePtr)
    {
        addLineDirective(code_, codePtr->startLineNumber(), dict.name());
    }

    if (includePtr)
    {
        addLineDirective(include_, includePtr->startLineNumber(), dict.name());
    }

    // Do not add line directive to options_ (Make/options) and libs since
    // they are preprocessed as a single line at this point. Can be fixed.
    if (localPtr)
    {
        addLineDirective(localCode_, localPtr->startLineNumber(), dict.name());
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::dynamicCodeContext::addLineDirective
(
    string& code,
    const label lineNum,
    const fileName& name
)
{
    code = "#line " + Foam::name(lineNum + 1) + " \"" + name + "\"\n" + code;
}


// ************************************************************************* //
