/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "inputModeEntry.H"
#include "dictionary.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::inputModeEntry::typeName
(
    Foam::functionEntries::inputModeEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include inputModeEntries
int Foam::functionEntries::inputModeEntry::debug(0);

namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeEntry,
        insert,
        primitiveEntryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeEntry,
        insert,
        dictionaryIstream
    );
}
}

// * * * * * * * * * * * * * * * * Private Data  * * * * * * * * * * * * * * //

Foam::label Foam::functionEntries::inputModeEntry::inputMode_ = imError;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionEntries::inputModeEntry::setMode(Istream& is)
{
    clear();

    word mode(is);
    if (mode == "merge")
    {
        inputMode_ = imMerge;
    }
    else if (mode == "overwrite")
    {
        inputMode_ = imOverwrite;
    }
    else if (mode == "error" || mode == "default")
    {
        inputMode_ = imError;
    }
    else
    {
        WarningIn("Foam::functionEntries::inputModeEntry::setMode(Istream&)")
            << "unsupported input mode " << mode
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::inputModeEntry::insert
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    setMode(is);
    return true;
}


bool Foam::functionEntries::inputModeEntry::insert
(
    dictionary& parentDict,
    Istream& is
)
{
    setMode(is);
    return true;
}


void Foam::functionEntries::inputModeEntry::clear()
{
    inputMode_ = imError;
}


bool Foam::functionEntries::inputModeEntry::merge()
{
    if (inputMode_ & imMerge)
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::functionEntries::inputModeEntry::overwrite()
{
    if (inputMode_ & imOverwrite)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
