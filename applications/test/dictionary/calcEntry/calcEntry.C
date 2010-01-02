/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "calcEntry.H"
#include "dictionary.H"
#include "addToMemberFunctionSelectionTable.H"

#include "ISstream.H"
#include "CocoParserErrors.H"
#include "calcEntryParser.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(calcEntry, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        calcEntry,
        execute,
        primitiveEntryIstream
    );

}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::calcEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    std::istream& iss = dynamicCast<ISstream>(is).stdStream();

    // define parser error handler
    CocoParserErrors<calcEntryInternal::Errors>
        myErrorHandler("calcEntryInternal::Parser");

    calcEntryInternal::Scanner scanner(iss);
    calcEntryInternal::Parser  parser(&scanner, &myErrorHandler);

    // Attach dictionary context
    parser.dict(parentDict);

    // Attach scalar functions
    // parser.functions(parentDict);

    parser.Parse();

    // make a small input list to contain the answer
    tokenList tokens(2);
    tokens[0] = parser.Result();
    tokens[1] = token::END_STATEMENT;

    entry.read(parentDict, ITstream("ParserResult", tokens)());

    return true;
}


// ************************************************************************* //
