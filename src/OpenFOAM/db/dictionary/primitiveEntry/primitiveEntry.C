/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "dictionary.H"
#include "OSspecific.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveEntry::append(const UList<token>& varTokens)
{
    for (const token& tok : varTokens)
    {
        newElmt(tokenIndex()++) = tok;
    }
}


bool Foam::primitiveEntry::expandVariable
(
    const string& varName,
    const dictionary& dict
)
{
    if (varName.size() > 1 && varName[0] == token::BEGIN_BLOCK)
    {
        // Recursive substitution mode.
        // Content between {} is replaced with expansion.
        string expanded(varName.substr(1, varName.size()-2));

        // Substitute dictionary and environment variables.
        // Do not allow empty substitutions.
        stringOps::inplaceExpand(expanded, dict, true, false);

        return expandVariable(expanded, dict);
    }

    // Lookup variable name in the given dictionary WITHOUT pattern matching.
    // Having a pattern match means that in this example:
    // {
    //     internalField XXX;
    //     boundaryField { ".*" {YYY;} movingWall {value $internalField;}
    // }
    // The $internalField would be matched by the ".*" !!!

    // Recursive, non-patterns
    const entry* eptr = dict.lookupScopedEntryPtr(varName, true, false);
    if (!eptr)
    {
        // Not found - revert to environment variable
        const string str(getEnv(varName));

        if (str.empty())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Illegal dictionary entry or environment variable name "
                << varName << endl << "Valid dictionary entries are "
                << dict.toc() << exit(FatalIOError);

            return false;
        }

        append(tokenList(IStringStream('(' + str + ')')()));
    }
    else if (eptr->isDict())
    {
        // Found dictionary entry
        append(eptr->dict().tokens());
    }
    else
    {
        // Found primitive entry
        append(eptr->stream());
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primitiveEntry::primitiveEntry(const keyType& key, const ITstream& is)
:
    entry(key),
    ITstream(is)
{
    name() += '.' + keyword();
}


Foam::primitiveEntry::primitiveEntry(const keyType& key, const token& t)
:
    entry(key),
    ITstream(key, tokenList(1, t))
{}


Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const UList<token>& tokens
)
:
    entry(key),
    ITstream(key, tokens)
{}


Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const Xfer<List<token>>& tokens
)
:
    entry(key),
    ITstream(key, tokens)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::primitiveEntry::startLineNumber() const
{
    const tokenList& tokens = *this;

    if (tokens.empty())
    {
        return -1;
    }
    else
    {
        return tokens.first().lineNumber();
    }
}


Foam::label Foam::primitiveEntry::endLineNumber() const
{
    const tokenList& tokens = *this;

    if (tokens.empty())
    {
        return -1;
    }
    else
    {
        return tokens.last().lineNumber();
    }
}


Foam::ITstream& Foam::primitiveEntry::stream() const
{
    ITstream& is = const_cast<primitiveEntry&>(*this);
    is.rewind();
    return is;
}


const Foam::dictionary& Foam::primitiveEntry::dict() const
{
    FatalErrorInFunction
        << "Attempt to return primitive entry " << info()
        << " as a sub-dictionary"
        << abort(FatalError);

    return dictionary::null;
}


Foam::dictionary& Foam::primitiveEntry::dict()
{
    FatalErrorInFunction
        << "Attempt to return primitive entry " << info()
        << " as a sub-dictionary"
        << abort(FatalError);

    return const_cast<dictionary&>(dictionary::null);
}


// ************************************************************************* //
