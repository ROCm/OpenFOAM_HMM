/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "dictionary.H"
#include "dictionaryEntry.H"
#include "stringOps.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    // file-scope
    //- Walk lists of patterns and regexps for an exact match
    //  or regular expression match
    template<class WcIterator, class ReIterator>
    static bool findInPatterns
    (
        const bool patternMatch,
        const word& keyword,
        WcIterator& wcIter,
        ReIterator& reIter
    )
    {
        while (wcIter.found())
        {
            if
            (
                patternMatch
              ? reIter()->match(keyword)
              : wcIter()->keyword() == keyword
            )
            {
                return true;
            }

            ++reIter;
            ++wcIter;
        }

        return false;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary::const_searcher Foam::dictionary::csearchDotScoped
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    std::string::size_type scopePos = keyword.find('.');

    if (scopePos == string::npos)
    {
        // Normal, non-scoped search
        return csearch(keyword, recursive, patternMatch);
    }

    if (scopePos == 0)
    {
        // Starting with a '.' -> go up for every further '.' found
        ++scopePos;

        const dictionary* dictPtr = this;
        for
        (
            string::const_iterator it = keyword.begin()+1;
            it != keyword.end() && *it == '.';
            ++scopePos, ++it
        )
        {
            // Go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = &dictPtr->parent_;
            }
            else
            {
                FatalIOErrorInFunction
                (
                    *this
                )   << "No parent of current dictionary when searching for "
                    << keyword.substr(1)
                    << exit(FatalIOError);

                return nullptr;
            }
        }

        return dictPtr->csearchDotScoped
        (
            keyword.substr(scopePos),
            false,
            patternMatch
        );
    }

    // The first word
    const_searcher finder = csearchDotScoped
    (
        keyword.substr(0, scopePos),
        false,
        patternMatch
    );

    // Fall back to finding key with '.' so e.g. if keyword is
    // a.b.c.d it would try
    // a.b, a.b.c, a.b.c.d

    if (!finder.found())
    {
        while (!finder.isDict())
        {
            scopePos = keyword.find('.', scopePos+1);

            // Local entry:
            finder = csearch(keyword.substr(0, scopePos), false, patternMatch);

            if (scopePos == string::npos)
            {
                // Parsed the whole word. Return entry or null.
                return finder;
            }
        }
    }

    if (finder.isDict())
    {
        return finder.dict().csearchDotScoped
        (
            keyword.substr(scopePos),
            false,
            patternMatch
        );
    }

    return finder;
}


Foam::dictionary::const_searcher Foam::dictionary::csearchSlashScoped
(
    const word& keyword,
    bool patternMatch
) const
{
    const dictionary* dictPtr = this;

    const auto slash = keyword.find('/');

    if (slash == string::npos)
    {
        // No slashes:
        // Can use normal (non-scoped) search at the current dictionary level
        return csearch(keyword, false, patternMatch);
    }
    else if (slash == 0)
    {
        // (isAbsolute)
        // Ascend to top-level
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }
    }

    // Split on '/'
    auto cmpts = stringOps::split<std::string>(keyword, '/');
    auto remaining = cmpts.size();

    for (const auto& cmpt : cmpts)
    {
        --remaining; // Decrement now so we can check (remaining == 0)

        if (cmpt == ".")
        {
            // "." - ignore
        }
        else if (cmpt == "..")
        {
            // ".." - go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = &dictPtr->parent_;
            }
            else
            {
                FatalIOErrorInFunction
                (
                    *dictPtr
                )   << "No parent of current dictionary when searching for "
                    << keyword << " at " << cmpt
                    << exit(FatalIOError);
                break;
            }
        }
        else
        {
            // Find entry
            const word key = word::validate(cmpt);

            auto finder = dictPtr->csearch(key, false, patternMatch);

            if (finder.found())
            {
                if (remaining)
                {
                    // Intermediate must be a dictionary
                    if (finder.isDict())
                    {
                        dictPtr = finder.dictPtr();
                    }
                    else
                    {
                        return const_searcher(dictPtr);
                    }
                }
                else
                {
                    // Last entry - done
                    return finder;
                }
            }
            else
            {
                break;
            }
        }
    }

    // Failed at this dictionary level
    return const_searcher(dictPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary::const_searcher Foam::dictionary::csearch
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    const_searcher finder(this);

    auto iter = hashedEntries_.cfind(keyword);

    if (iter.found())
    {
        finder.set(iter.object());
        return finder;
    }

    if (patternMatch && patterns_.size())
    {
        pattern_const_iterator wcLink = patterns_.begin();
        regexp_const_iterator  reLink = regexps_.begin();

        // Find in patterns using regular expressions only
        if (findInPatterns(patternMatch, keyword, wcLink, reLink))
        {
            finder.set(*wcLink);
            return finder;
        }
    }

    if (recursive && &parent_ != &dictionary::null)
    {
        return parent_.csearch
        (
            keyword,
            recursive,
            patternMatch
        );
    }

    return finder;
}


Foam::dictionary::const_searcher Foam::dictionary::search
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    return csearch(keyword, recursive, patternMatch);
}


Foam::dictionary::searcher Foam::dictionary::search
(
    const word& keyword,
    bool recursive,
    bool patternMatch
)
{
    const_searcher finder = csearch(keyword, recursive, patternMatch);

    return static_cast<const searcher&>(finder);
}


Foam::dictionary::const_searcher Foam::dictionary::csearchScoped
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    if (keyword.find('/') != string::npos)
    {
        return csearchSlashScoped(keyword, patternMatch);
    }

    if (keyword[0] == ':' || keyword[0] == '^')
    {
        // Ascend to top-level
        const dictionary* dictPtr = this;
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }

        return dictPtr->csearchDotScoped
        (
            keyword.substr(1),
            false,
            patternMatch
        );
    }

    return csearchDotScoped(keyword, recursive, patternMatch);
}


Foam::dictionary::const_searcher Foam::dictionary::searchScoped
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    return csearchScoped(keyword, recursive, patternMatch);
}


Foam::dictionary::searcher Foam::dictionary::searchScoped
(
    const word& keyword,
    bool recursive,
    bool patternMatch
)
{
    const_searcher finder = csearchScoped(keyword, recursive, patternMatch);

    return static_cast<const searcher&>(finder);
}


const Foam::dictionary* Foam::dictionary::cfindScopedDictPtr
(
    const fileName& dictPath
) const
{
    // Or warning
    if (dictPath.empty())
    {
        return nullptr;
    }

    const dictionary* dictPtr = this;
    if (fileName::isAbsolute(dictPath))
    {
        // Ascend to top-level
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }
    }

    fileName path = dictPath.clean();
    const wordList cmpts = path.components();

    for (const word& cmpt : cmpts)
    {
        if (cmpt == ".")
        {
            // "." - ignore
        }
        else if (cmpt == "..")
        {
            // ".." - go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = &dictPtr->parent_;
            }
            else
            {
                FatalIOErrorInFunction
                (
                    *dictPtr
                )   << "No parent for dictionary while searching "
                    << path
                    << exit(FatalIOError);

                return nullptr;
            }
        }
        else
        {
            // Non-recursive, no patternMatch
            // -> can do direct lookup, without csearch(cmpt, false, false);

            auto iter = dictPtr->hashedEntries_.cfind(cmpt);

            if (iter.found())
            {
                const entry *eptr = iter.object();

                if (eptr->isDict())
                {
                    dictPtr = eptr->dictPtr();
                }
                else
                {
                    FatalIOErrorInFunction
                    (
                        *dictPtr
                    )   << "Found entry '" << cmpt
                        << "' but not a dictionary, while searching scoped"
                        << nl
                        << "    " << path
                        << exit(FatalIOError);

                    return nullptr;
                }
            }
            else
            {
                return nullptr;
            }
        }
    }

    return dictPtr;
}


const Foam::dictionary* Foam::dictionary::findScopedDictPtr
(
    const fileName& dictPath
) const
{
    return cfindScopedDictPtr(dictPath);
}


Foam::dictionary* Foam::dictionary::findScopedDictPtr
(
    const fileName& dictPath
)
{
    const dictionary* ptr = cfindScopedDictPtr(dictPath);
    return const_cast<dictionary*>(ptr);
}


Foam::dictionary* Foam::dictionary::makeScopedDictPtr(const fileName& dictPath)
{
    // Or warning
    if (dictPath.empty())
    {
        return nullptr;
    }

    dictionary* dictPtr = this;
    if (fileName::isAbsolute(dictPath))
    {
        // Ascend to top-level
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = const_cast<dictionary*>(&dictPtr->parent_);
        }
    }

    // Work on a copy, without any assumptions
    std::string path = dictPath;
    fileName::clean(path);

    // Split on '/'
    auto cmpts = stringOps::split(path, '/');

    for (const auto& cmpt : cmpts)
    {
        if (cmpt == ".")
        {
            // "." - ignore
        }
        else if (cmpt == "..")
        {
            // ".." - go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = const_cast<dictionary*>(&dictPtr->parent_);
            }
            else
            {
                FatalIOErrorInFunction
                (
                    *dictPtr
                )   << "No parent for dictionary while searching "
                    << path
                    << exit(FatalIOError);

                return nullptr;
            }
        }
        else
        {
            // Non-recursive, no patternMatch
            // -> can do direct lookup, without csearch(cmptName, false, false);
            const word cmptName(cmpt.str(), false);

            auto iter = dictPtr->hashedEntries_.find(cmptName);

            if (iter.found())
            {
                entry *eptr = iter.object();

                if (eptr->isDict())
                {
                    dictPtr = eptr->dictPtr();
                }
                else
                {
                    FatalIOErrorInFunction
                    (
                        *dictPtr
                    )   << "Cannot create sub-dictionary entry '" << cmptName
                        << "' - a non-dictionary entry is in the way"
                        << nl << "Encountered in scope" << nl
                        << "    " << path
                        << exit(FatalIOError);

                    return nullptr;
                }
            }
            else
            {
                dictionaryEntry *eptr =
                    new dictionaryEntry(cmptName, *dictPtr, dictionary());

                // Add *without* merging, since we just checked that the entry
                // doesn't exist and to ensure that the pointer remains valid.

                if (dictPtr->add(eptr, false))  // NO merge
                {
                    dictPtr = eptr;
                }
                else
                {
                    // Note: a failed add() deletes the eptr passed
                    return nullptr;
                }
            }
        }
    }

    return dictPtr;
}


bool Foam::dictionary::remove(const word& keyword)
{
    auto iter = hashedEntries_.find(keyword);

    if (iter.found())
    {
        // Delete from patterns
        pattern_iterator wcLink = patterns_.begin();
        regexp_iterator  reLink = regexps_.begin();

        // Find in pattern using exact match only
        if (findInPatterns(false, keyword, wcLink, reLink))
        {
            patterns_.remove(wcLink);
            regexps_.remove(reLink);
        }

        parent_type::remove(iter());
        delete iter();
        hashedEntries_.erase(iter);

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::dictionary::changeKeyword
(
    const keyType& oldKeyword,
    const keyType& newKeyword,
    bool overwrite
)
{
    // No change
    if (oldKeyword == newKeyword)
    {
        return false;
    }

    // Check that oldKeyword exists and can be changed
    auto iter = hashedEntries_.find(oldKeyword);

    if (!iter.found())
    {
        return false;
    }

    if (iter()->keyword().isPattern())
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "Old keyword "<< oldKeyword
            << " is a pattern."
            << "Pattern replacement not yet implemented."
            << exit(FatalIOError);
    }


    auto iter2 = hashedEntries_.find(newKeyword);

    // newKeyword already exists
    if (iter2.found())
    {
        if (overwrite)
        {
            if (iter2()->keyword().isPattern())
            {
                // Delete from patterns
                pattern_iterator wcLink = patterns_.begin();
                regexp_iterator  reLink = regexps_.begin();

                // Find in patterns using exact match only
                if (findInPatterns(false, iter2()->keyword(), wcLink, reLink))
                {
                    patterns_.remove(wcLink);
                    regexps_.remove(reLink);
                }
            }

            parent_type::replace(iter2(), iter());
            delete iter2();
            hashedEntries_.erase(iter2);
        }
        else
        {
            IOWarningInFunction
            (
                *this
            )   << "cannot rename keyword "<< oldKeyword
                << " to existing keyword " << newKeyword
                << " in dictionary " << name() << endl;
            return false;
        }
    }

    // Change name and HashTable, but leave DL-List untouched
    iter()->keyword() = newKeyword;
    iter()->name() = name() + '.' + newKeyword;
    hashedEntries_.erase(oldKeyword);
    hashedEntries_.insert(newKeyword, iter());

    if (newKeyword.isPattern())
    {
        patterns_.insert(iter());
        regexps_.insert
        (
            autoPtr<regExp>(new regExp(newKeyword))
        );
    }

    return true;
}


// ************************************************************************* //
