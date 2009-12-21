/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
#include "IStringStream.H"
#include "OStringStream.H"
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
    Istream& istr
)
{
    static const int maxLen = 1024;
    static const int errLen = 80; // truncate error message for readability
    static char buf[maxLen];

    ISstream& is = dynamicCast<ISstream>(istr);

    // get the { ... } argument without the enclosing brace brackets
    //
    // THIS NEEDS REWORKING (LATER) ...

    char c = 0;

    if (!is.read(c).good() || c != token::BEGIN_BLOCK)
    {
        is.setBad();
        FatalIOErrorIn("functionEntries::calcEntry::execute()", is)
            << "Expected a '" << token::BEGIN_BLOCK
            << "', found '" << c << "'" << exit(FatalIOError);

        return false;
    }

    register int nChar = 0;
    buf[nChar++] = token::BEGIN_BLOCK;
    int listDepth = 1;       // already saw the first '{'

    while (is.get(c).good())
    {
        buf[nChar++] = c;
        if (nChar == maxLen)
        {
            buf[errLen] = '\0';

            FatalIOErrorIn("functionEntries::calcEntry::execute()", is)
                << "argument \"" << buf << "...\"\n"
                << "    is too long (max. " << maxLen << " characters)"
                << exit(FatalIOError);

            return false;
        }

        // handle nested blocks, even if we don't know what they'd
        // be useful for
        if (c == token::BEGIN_BLOCK)
        {
            ++listDepth;
        }
        else if (c == token::END_BLOCK)
        {
            if (--listDepth == 0)
            {
                // done reading - overwrite the final '}'
                // --nChar;
                break;
            }
        }
    }

    buf[nChar] = '\0';


    // emit some info
    Info<< "grabbed " << nChar << " characters:" << nl
        << "----------\n"
        << buf << nl
        << "----------\n"
        << nl;


    // define parser error handler
    CocoParserErrors<calcEntryInternal::Errors>
        myErrorHandler("calcEntry::Parser--");

    calcEntryInternal::Scanner scanner(buf, nChar);
    calcEntryInternal::Parser  parser(&scanner, &myErrorHandler);

    // Attach dictionary context
    parser.dict(parentDict);

    parser.Parse();

//    Info<<"got: " << parser.Result() << endl;

    tokenList tokens(2);
    tokens[0] = parser.Result();
    tokens[1] = token::END_STATEMENT;

//     Info<<"tokens[0] = " << tokens[0].info() <<nl;
//     Info<<"tokens[1] = " << tokens[1].info() <<nl;
//
//     {
//     const tokenList& toks = entry;
//
//     forAll(toks, tokI)
//     {
//         Info<< tokI <<":= " << toks[tokI].info() << endl;
//     }
//     }

    ITstream its("ParserResult", tokens);
    entry.read(parentDict, its);
//    Info<< "size: " << entry.size() << endl;
//    entry = newente;
//    entry.print(Info);

///     const tokenList& toks = entry;
///
///     forAll(toks, tokI)
///     {
///         Info<< tokI <<":= " << toks[tokI].info() << endl;
///     }
///
    return true;
}


// ************************************************************************* //
