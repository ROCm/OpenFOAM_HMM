/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Foam::functionEntries::inputModeEntry::inputMode
    Foam::functionEntries::inputModeEntry::mode_(MERGE);

namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeEntry,
        execute,
        dictionaryIstream
    );
}
}


const Foam::Enum
<
    Foam::functionEntries::inputModeEntry::inputMode
>
Foam::functionEntries::inputModeEntry::inputModeNames
{
    { inputMode::MERGE,  "merge" },
    { inputMode::OVERWRITE, "overwrite" },
    { inputMode::PROTECT, "protect" },
    { inputMode::WARN, "warn" },
    { inputMode::ERROR, "error" },
    // Aliases
    { inputMode::MERGE, "default" },
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::inputModeEntry::execute
(
    dictionary& unused,
    Istream& is
)
{
    const word modeName(is);

    // Bheaviour like Enum lookupOrFailsafe()
    if (inputModeNames.hasEnum(modeName))
    {
        mode_ = inputModeNames[modeName];
    }
    else
    {
        WarningInFunction
            << "Unsupported inputMode '" << modeName
            << "' ... defaulting to 'merge'"
            << endl;

        reset();
    }

    return true;
}


void Foam::functionEntries::inputModeEntry::reset()
{
    mode_ = MERGE;
}


void Foam::functionEntries::inputModeEntry::clear()
{
    reset();
}


bool Foam::functionEntries::inputModeEntry::merge()
{
    return mode_ == MERGE;
}


bool Foam::functionEntries::inputModeEntry::overwrite()
{
    return mode_ == OVERWRITE;
}


bool Foam::functionEntries::inputModeEntry::protect()
{
    return mode_ == PROTECT;
}

bool Foam::functionEntries::inputModeEntry::error()
{
    return mode_ == ERROR;
}


// ************************************************************************* //
