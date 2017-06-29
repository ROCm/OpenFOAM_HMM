/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "TimePaths.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::TimePaths::detectProcessorCase()
{
    if (processorCase_)
    {
        return processorCase_;
    }

    // Look for "processor", but should really check for following digits too
    const std::string::size_type sep = globalCaseName_.rfind('/');
    const std::string::size_type pos = globalCaseName_.find
    (
        "processor",
        (sep == string::npos ? 0 : sep)
    );

    if (pos == 0)
    {
        globalCaseName_ = ".";
        processorCase_  = true;
    }
    else if (pos != string::npos && sep != string::npos && sep == pos-1)
    {
        globalCaseName_.resize(sep);
        processorCase_  = true;
    }

    return processorCase_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TimePaths::TimePaths
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    processorCase_(false),
    rootPath_(rootPath),
    distributed_(false),
    globalCaseName_(caseName),
    case_(caseName),
    system_(systemName),
    constant_(constantName)
{
    // Find out from case name whether a processor directory
    detectProcessorCase();
}


Foam::TimePaths::TimePaths
(
    const bool processorCase,
    const fileName& rootPath,
    const bool distributed,
    const fileName& globalCaseName,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    processorCase_(processorCase),
    rootPath_(rootPath),
    distributed_(distributed),
    globalCaseName_(globalCaseName),
    case_(caseName),
    system_(systemName),
    constant_(constantName)
{
    // For convenience: find out from case name whether it is a
    // processor directory and set processorCase flag so file searching
    // goes up one level.
    detectProcessorCase();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::TimePaths::caseSystem() const
{
    if (processorCase_)
    {
        return ".."/system();
    }
    else
    {
        return system();
    }
}


Foam::fileName Foam::TimePaths::caseConstant() const
{
    if (processorCase_)
    {
        return ".."/constant();
    }
    else
    {
        return constant();
    }
}



// ************************************************************************* //
