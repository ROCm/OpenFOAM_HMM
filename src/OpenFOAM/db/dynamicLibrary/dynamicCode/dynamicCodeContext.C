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
    code_(stringOps::trim(dict["code"])),
    localCode_(),
    include_(),
    options_()
{
    // expand dictionary entries
    stringOps::inplaceExpand(code_, dict);

    // note: removes any leading/trailing whitespace
    // - necessary for compilation options, convenient for includes
    // and body.

    // optional
    if (dict.found("localCode"))
    {
        localCode_ = stringOps::trim(dict["localCode"]);
        stringOps::inplaceExpand(localCode_, dict);
    }

    // optional
    if (dict.found("codeInclude"))
    {
        include_ = stringOps::trim(dict["codeInclude"]);
        stringOps::inplaceExpand(include_, dict);
    }

    // optional
    if (dict.found("codeOptions"))
    {
        options_ = stringOps::trim(dict["codeOptions"]);
        stringOps::inplaceExpand(options_, dict);
    }

    // calculate SHA1 digest from include, options, localCode, code
    OSHA1stream os;
    os  << include_ << options_ << localCode_ << code_;
    sha1_ = os.digest();
}


// ************************************************************************* //
