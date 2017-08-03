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

#include "entry.H"
#include "primitiveEntry.H"
#include "dictionaryEntry.H"
#include "functionEntry.H"
#include "includeEntry.H"
#include "stringOps.H"
#include "dictionaryListEntry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::entry::getKeyword(keyType& keyword, token& keyToken, Istream& is)
{
    // Read the next valid token discarding spurious ';'s
    do
    {
        if
        (
            is.read(keyToken).bad()
         || is.eof()
         || !keyToken.good()
        )
        {
            return false;
        }
    }
    while (keyToken == token::END_STATEMENT);

    // If the token is a valid keyword set 'keyword' return true...
    if (keyToken.isWord())
    {
        keyword = keyToken.wordToken();
        return true;
    }
    else if (keyToken.isString())
    {
        // Enable wildcards
        keyword = keyToken.stringToken();
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::entry::getKeyword(keyType& keyword, Istream& is)
{
    token keyToken;
    const bool valid = getKeyword(keyword, keyToken, is);

    if (valid)
    {
        return true;
    }

    // Do some more checking
    if (keyToken == token::END_BLOCK || is.eof())
    {
        return false;
    }
    else
    {
        // Otherwise the token is invalid
        cerr<< "--> FOAM Warning :" << nl
            << "    From function "
            << FUNCTION_NAME << nl
            << "    in file " << __FILE__
            << " at line " << __LINE__ << nl
            << "    Reading " << is.name().c_str() << nl
            << "    found " << keyToken << nl
            << "    expected either " << token::END_BLOCK << " or EOF"
            << std::endl;
        return false;
    }
}


bool Foam::entry::New
(
    dictionary& parentDict,
    Istream& is,
    const entry::inputMode inMode
)
{
    // The inputMode for dealing with duplicate entries
    const entry::inputMode mode =
    (
        inMode == inputMode::GLOBAL
      ? globalInputMode
      : inMode
    );

    // If somehow the global itself is 'global' - this is a severe logic error.
    if (mode == inputMode::GLOBAL)
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot use 'GLOBAL' as an inputMode"
            << exit(FatalIOError);
    }

    is.fatalCheck(FUNCTION_NAME);

    keyType keyword;
    token keyToken;

    // Get the next keyword and if a valid keyword return true
    const bool valid = getKeyword(keyword, keyToken, is);

    if (!valid)
    {
        // Do some more checking
        if (keyToken == token::END_BLOCK || is.eof())
        {
            return false;
        }
        else if
        (
            keyToken.isLabel()
         || (keyToken.isPunctuation() && keyToken.pToken() == token::BEGIN_LIST)
        )
        {
            is.putBack(keyToken);
            return parentDict.add
            (
                new dictionaryListEntry(parentDict, is),
                false
            );
        }
        else
        {
            // Otherwise the token is invalid
            cerr<< "--> FOAM Warning :" << nl
                << "    From function "
                << FUNCTION_NAME << nl
                << "    in file " << __FILE__
                << " at line " << __LINE__ << nl
                << "    Reading " << is.name().c_str() << nl
                << "    found " << keyToken << nl
                << "    expected either " << token::END_BLOCK << " or EOF"
                << std::endl;
            return false;
        }
    }
    else if (keyword[0] == '#')
    {
        // Function entry

        if (disableFunctionEntries)
        {
            return parentDict.add
            (
                new functionEntry
                (
                    keyword,
                    parentDict,
                    is
                ),
                false
            );
        }
        else
        {
            const word functionName(keyword.substr(1), false);
            return functionEntry::execute(functionName, parentDict, is);
        }
    }
    else if
    (
        !disableFunctionEntries
     && keyword[0] == '$'
    )
    {
        // Substitution entry

        token nextToken(is);
        is.putBack(nextToken);

        if (keyword.size() > 2 && keyword[1] == token::BEGIN_BLOCK)
        {
            // Recursive substitution mode.
            // Content between {} is replaced with expansion.
            // Then let standard variable expansion deal with rest.
            string expanded = keyword.substr(2, keyword.size()-3);

            // Substitute dictionary and environment variables.
            // Do not allow empty substitutions.
            stringOps::inplaceExpand(expanded, parentDict, true, false);

            // Restore the '$' prefix.
            // Use replace since operator= is private
            keyword.std::string::replace(1, keyword.size()-1, expanded);
        }

        if (nextToken == token::BEGIN_BLOCK)
        {
            const word varName = keyword.substr(1);

            // Lookup the variable name in the given dictionary
            const entry* ePtr = parentDict.lookupScopedEntryPtr
            (
                varName,
                true,
                true
            );

            if (ePtr)
            {
                // Read as primitiveEntry
                const keyType newKeyword(ePtr->stream());

                return parentDict.add
                (
                    new dictionaryEntry(newKeyword, parentDict, is),
                    false
                );
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << "Attempt to use undefined variable " << varName
                    << " as keyword"
                    << exit(FatalIOError);
                return false;
            }
        }
        else
        {
            // Deal with duplicate entries (at least partially)
            const bool mergeEntry =
            (
                mode == inputMode::MERGE
             || mode == inputMode::OVERWRITE
            );

            parentDict.substituteScopedKeyword(keyword, mergeEntry);
        }

        return true;
    }
    else
    {
        // Normal entry

        token nextToken(is);
        is.putBack(nextToken);

        if (nextToken == token::END_LIST)
        {
            FatalIOErrorInFunction(is)
                << "Unexpected token encountered for "
                << keyword << " - " << nextToken.info()
                << exit(FatalIOError);
            return false;
        }

        // How to manage duplicate entries
        bool mergeEntry = false;

        // See (using exact match) if entry already present
        entry* existingPtr = parentDict.lookupEntryPtr
        (
            keyword,
            false,
            false
        );

        if (existingPtr)
        {
            if (mode == inputMode::MERGE)
            {
                mergeEntry = true;
            }
            else if (mode == inputMode::OVERWRITE)
            {
                // Clear existing dictionary so merge acts like overwrite
                if (existingPtr->isDict())
                {
                    existingPtr->dict().clear();
                }
                mergeEntry = true;
            }
            else if (mode == inputMode::PROTECT)
            {
                // Read and discard the entry.
                // Disable function/variable expansion to avoid side-effects
                const int oldFlag = entry::disableFunctionEntries;
                entry::disableFunctionEntries = 1;

                if (nextToken == token::BEGIN_BLOCK)
                {
                    dictionaryEntry dummy("dummy", parentDict, is);
                }
                else
                {
                    primitiveEntry  dummy("dummy", parentDict, is);
                }

                entry::disableFunctionEntries = oldFlag;
                return true;
            }
            else if (mode == inputMode::ERROR)
            {
                FatalIOErrorInFunction(is)
                    << "duplicate entry: " << keyword
                    << exit(FatalIOError);

                return false;
            }
        }


        if (nextToken == token::BEGIN_BLOCK)
        {
            return parentDict.add
            (
                new dictionaryEntry(keyword, parentDict, is),
                mergeEntry
            );
        }
        else
        {
            return parentDict.add
            (
                new primitiveEntry(keyword, parentDict, is),
                mergeEntry
            );
        }
    }
}


Foam::autoPtr<Foam::entry> Foam::entry::New(Istream& is)
{
    is.fatalCheck(FUNCTION_NAME);

    autoPtr<entry> ptr(nullptr);

    // Get the next keyword and if invalid return false
    keyType keyword;
    if (getKeyword(keyword, is))
    {
        // Keyword starts entry ...
        token nextToken(is);
        is.putBack(nextToken);

        if (nextToken == token::BEGIN_BLOCK)
        {
            // A sub-dictionary
            ptr.reset(new dictionaryEntry(keyword, dictionary::null, is));
        }
        else
        {
            ptr.reset(new primitiveEntry(keyword, is));
        }
    }

    return ptr;
}


// * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const entry& e)
{
    e.write(os);
    return os;
}


// ************************************************************************* //
