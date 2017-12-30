/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "primitiveEntry.H"
#include "functionEntry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::primitiveEntry::acceptToken
(
    const token& tok,
    const dictionary& dict,
    Istream& is
)
{
    bool accept = tok.good();

    if (tok.isWord())
    {
        const word& key = tok.wordToken();

        accept =
        (
            disableFunctionEntries
         || key.size() == 1
         || (
                !(key[0] == '$' && expandVariable(key.substr(1), dict))
             && !(key[0] == '#' && expandFunction(key.substr(1), dict, is))
            )
        );
    }
    else if (tok.isVariable())
    {
        const string& key = tok.stringToken();

        accept =
        (
            disableFunctionEntries
         || key.size() <= 3
         || !(
                key[0] == '$'
             && key[1] == token::BEGIN_BLOCK
             && expandVariable(key.substr(1), dict)
            )
        );
    }

    return accept;
}


bool Foam::primitiveEntry::expandFunction
(
    const word& functionName,
    const dictionary& dict,
    Istream& is
)
{
    return functionEntry::execute(functionName, dict, *this, is);
}


bool Foam::primitiveEntry::read(const dictionary& dict, Istream& is)
{
    is.fatalCheck(FUNCTION_NAME);

    label depth = 0;
    token tok;

    while
    (
        !is.read(tok).bad() && tok.good()
     && !(tok == token::END_STATEMENT && depth == 0)
    )
    {
        if (tok.isPunctuation())
        {
            const char c = tok.pToken();
            if (c == token::BEGIN_BLOCK || c == token::BEGIN_LIST)
            {
                ++depth;
            }
            else if (c == token::END_BLOCK || c == token::END_LIST)
            {
                --depth;
            }
        }

        if (acceptToken(tok, dict, is))
        {
            newElmt(tokenIndex()++) = std::move(tok);
        }

        // With/without move: clear any old content and force to have a
        // known good token so that we can rely on it for the return value.

        tok = token::punctuationToken::NULL_TOKEN;
    }

    is.fatalCheck(FUNCTION_NAME);
    return tok.good();
}


void Foam::primitiveEntry::readEntry(const dictionary& dict, Istream& is)
{
    label keywordLineNumber = is.lineNumber();
    tokenIndex() = 0;

    if (read(dict, is))
    {
        setSize(tokenIndex());
        tokenIndex() = 0;
    }
    else
    {
        std::ostringstream os;
        os  << "ill defined primitiveEntry starting at keyword '"
            << keyword() << '\''
            << " on line " << keywordLineNumber
            << " and ending at line " << is.lineNumber();

        SafeFatalIOErrorInFunction
        (
            is,
            os.str()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const dictionary& dict,
    Istream& is
)
:
    entry(key),
    ITstream
    (
        is.name() + '.' + key,
        tokenList(10),
        is.format(),
        is.version()
    )
{
    readEntry(dict, is);
}


Foam::primitiveEntry::primitiveEntry(const keyType& key, Istream& is)
:
    entry(key),
    ITstream
    (
        is.name() + '.' + key,
        tokenList(10),
        is.format(),
        is.version()
    )
{
    readEntry(dictionary::null, is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::primitiveEntry::write(Ostream& os, const bool contentsOnly) const
{
    if (!contentsOnly)
    {
        os.writeKeyword(keyword());
    }

    bool addSpace = false;  // Separate from previous tokens with a space
    for (const token& tok : *this)
    {
        if (addSpace) os << token::SPACE;

        // Try to output token directly, with special handling in Ostreams.

        if (!os.write(tok))
        {
            os  << tok;   // Revert to normal '<<' output operator
        }

        addSpace = true;  // Separate from following tokens
    }

    if (!contentsOnly)
    {
        os  << token::END_STATEMENT << endl;
    }
}


void Foam::primitiveEntry::write(Ostream& os) const
{
    this->write(os, false);
}


// * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<primitiveEntry>& ip
)
{
    const primitiveEntry& e = ip.t_;

    e.print(os);

    const label nPrintTokens = 10;

    os  << "    primitiveEntry '" << e.keyword() << "' comprises ";

    for (label i=0; i<min(e.size(), nPrintTokens); i++)
    {
        os  << nl << "        " << e[i].info();
    }

    if (e.size() > nPrintTokens)
    {
        os  << " ...";
    }

    os  << endl;

    return os;
}


// ************************************************************************* //
