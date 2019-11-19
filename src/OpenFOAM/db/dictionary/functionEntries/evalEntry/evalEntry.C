/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "evalEntry.H"
#include "dictionary.H"
#include "stringOps.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        evalEntry,
        execute,
        dictionaryIstream,
        eval
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        evalEntry,
        execute,
        primitiveEntryIstream,
        eval
    );

} // End namespace functionEntry
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::functionEntries::evalEntry::evaluate
(
    const dictionary& parentDict,
    Istream& is
)
{
    #ifdef FULLDEBUG
    DetailInfo
        << "Using #eval at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << nl;
    #endif

    // String to evaluate
    string s;

    token tok(is);

    if (!tok.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get string to evaluate"
            << exit(FatalIOError);
        return 0;
    }

    if (tok.isString())
    {
        s = tok.stringToken();
    }
    else if (tok == token::BEGIN_BLOCK)
    {
        dynamic_cast<ISstream&>(is).getLine(s, token::END_BLOCK);
    }
    else
    {
        is.putBack(tok);

        FatalIOErrorInFunction(is)
            << "Invalid input for #eval" << nl
            << exit(FatalIOError);
    }

    #ifdef FULLDEBUG
    DetailInfo
        << "input: " << s << endl;
    #endif

    // Expanding with env-variables, with empty

    stringOps::inplaceRemoveComments(s);
    stringOps::inplaceExpand(s, parentDict, true, true);

    #ifdef FULLDEBUG
    DetailInfo
        << "expanded: " << s << endl;
    #endif

    return stringOps::toScalar(s);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::evalEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const scalar value = evaluate(parentDict, is);

    // The result as terminated token entry
    ITstream result
    (
        "eval",
        tokenList({token(value), token(token::END_STATEMENT)})
    );

    entry.read(parentDict, result);

    return true;
}


bool Foam::functionEntries::evalEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const scalar value = evaluate(parentDict, is);

    // The result as terminated token entry
    ITstream result
    (
        "eval",
        tokenList({token(value), token(token::END_STATEMENT)})
    );

    parentDict.read(result);

    return true;
}


// ************************************************************************* //
