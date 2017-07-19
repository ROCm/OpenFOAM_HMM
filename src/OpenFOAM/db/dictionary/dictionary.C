/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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
#include "primitiveEntry.H"
#include "dictionaryEntry.H"
#include "regExp.H"
#include "OSHA1stream.H"
#include "DynamicList.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(dictionary, 0);
}

const Foam::dictionary Foam::dictionary::null;

bool Foam::dictionary::writeOptionalEntries
(
    Foam::debug::infoSwitch("writeOptionalEntries", 0)
);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

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

const Foam::entry* Foam::dictionary::lookupScopedSubEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    string::size_type dotPos = keyword.find('.');

    if (dotPos == string::npos)
    {
        // Non-scoped lookup
        return lookupEntryPtr(keyword, recursive, patternMatch);
    }
    else if (dotPos == 0)
    {
        // Starting with a '.' -> go up for every further '.' found
        ++dotPos;

        const dictionary* dictPtr = this;
        for
        (
            string::const_iterator it = keyword.begin()+1;
            it != keyword.end() && *it == '.';
            ++dotPos, ++it
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

        return dictPtr->lookupScopedSubEntryPtr
        (
            keyword.substr(dotPos),
            false,
            patternMatch
        );
    }
    else
    {
        // The first word
        const entry* entPtr = lookupScopedSubEntryPtr
        (
            keyword.substr(0, dotPos),
            false,
            patternMatch
        );

        if (!entPtr)
        {
            // Fall back to finding key with '.' so e.g. if keyword is
            // a.b.c.d it would try
            // a.b, a.b.c, a.b.c.d

            while (true)
            {
                dotPos = keyword.find('.', dotPos+1);

                entPtr = lookupEntryPtr
                (
                    keyword.substr(0, dotPos),
                    false,
                    patternMatch
                );

                if (dotPos == string::npos)
                {
                    // Parsed the whole word. Return entry or null.
                    return entPtr;
                }

                if (entPtr && entPtr->isDict())
                {
                    return entPtr->dict().lookupScopedSubEntryPtr
                    (
                        keyword.substr(dotPos),
                        false,
                        patternMatch
                    );
                }
            }
        }

        if (entPtr->isDict())
        {
            return entPtr->dict().lookupScopedSubEntryPtr
            (
                keyword.substr(dotPos),
                false,
                patternMatch
            );
        }
    }

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionary::dictionary()
:
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary(const fileName& name)
:
    dictionaryName(name),
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    const dictionary& dict
)
:
    dictionaryName(dict.name()),
    parent_type(dict, *this),
    parent_(parentDict)
{
    forAllIter(parent_type, *this, iter)
    {
        hashedEntries_.insert(iter().keyword(), &iter());

        if (iter().keyword().isPattern())
        {
            patterns_.insert(&iter());
            regexps_.insert
            (
                autoPtr<regExp>(new regExp(iter().keyword()))
            );
        }
    }
}


Foam::dictionary::dictionary
(
    const dictionary& dict
)
:
    dictionaryName(dict.name()),
    parent_type(dict, *this),
    parent_(dictionary::null)
{
    forAllIter(parent_type, *this, iter)
    {
        hashedEntries_.insert(iter().keyword(), &iter());

        if (iter().keyword().isPattern())
        {
            patterns_.insert(&iter());
            regexps_.insert
            (
                autoPtr<regExp>(new regExp(iter().keyword()))
            );
        }
    }
}


Foam::dictionary::dictionary
(
    const dictionary* dictPtr
)
:
    parent_(dictionary::null)
{
    if (dictPtr)
    {
        operator=(*dictPtr);
    }
}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    const Xfer<dictionary>& dict
)
:
    parent_(parentDict)
{
    transfer(dict());
    name() = parentDict.name() + '.' + name();
}


Foam::dictionary::dictionary
(
    const Xfer<dictionary>& dict
)
:
    parent_(dictionary::null)
{
    transfer(dict());
}


Foam::autoPtr<Foam::dictionary> Foam::dictionary::clone() const
{
    return autoPtr<dictionary>(new dictionary(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dictionary::~dictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::dictionary::topDict() const
{
    const dictionary& p = parent();

    if (&p != this && !p.name().empty())
    {
        return p.topDict();
    }
    else
    {
        return *this;
    }
}


Foam::label Foam::dictionary::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::label Foam::dictionary::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::SHA1Digest Foam::dictionary::digest() const
{
    OSHA1stream os;

    // Process entries
    forAllConstIter(parent_type, *this, iter)
    {
        os << *iter;
    }

    return os.digest();
}


Foam::tokenList Foam::dictionary::tokens() const
{
    // Serialize dictionary into a string
    OStringStream os;
    write(os, false);
    IStringStream is(os.str());

    // Parse string as tokens
    DynamicList<token> tokens;
    token t;
    while (is.read(t))
    {
        tokens.append(t);
    }

    return tokenList(tokens.xfer());
}


bool Foam::dictionary::found
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    if (hashedEntries_.found(keyword))
    {
        return true;
    }

    if (patternMatch && patterns_.size())
    {
        pattern_const_iterator wcLink = patterns_.begin();
        regexp_const_iterator  reLink = regexps_.begin();

        // Find in patterns using regular expressions only
        if (findInPatterns(patternMatch, keyword, wcLink, reLink))
        {
            return true;
        }
    }

    if (recursive && &parent_ != &dictionary::null)
    {
        return parent_.found(keyword, recursive, patternMatch);
    }

    return false;
}


const Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    auto iter = hashedEntries_.cfind(keyword);

    if (iter.found())
    {
        return iter();
    }

    if (patternMatch && patterns_.size())
    {
        pattern_const_iterator wcLink = patterns_.begin();
        regexp_const_iterator  reLink = regexps_.begin();

        // Find in patterns using regular expressions only
        if (findInPatterns(patternMatch, keyword, wcLink, reLink))
        {
            return wcLink();
        }
    }

    if (recursive && &parent_ != &dictionary::null)
    {
        return parent_.lookupEntryPtr(keyword, recursive, patternMatch);
    }

    return nullptr;
}


Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
)
{
    auto iter = hashedEntries_.find(keyword);

    if (iter.found())
    {
        return iter();
    }

    if (patternMatch && patterns_.size())
    {
        pattern_iterator wcLink = patterns_.begin();
        regexp_iterator  reLink = regexps_.begin();

        // Find in patterns using regular expressions only
        if (findInPatterns(patternMatch, keyword, wcLink, reLink))
        {
            return wcLink();
        }
    }

    if (recursive && &parent_ != &dictionary::null)
    {
        return const_cast<dictionary&>(parent_).lookupEntryPtr
        (
            keyword,
            recursive,
            patternMatch
        );
    }

    return nullptr;
}


const Foam::entry& Foam::dictionary::lookupEntry
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return *entryPtr;
}


Foam::ITstream& Foam::dictionary::lookup
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    return lookupEntry(keyword, recursive, patternMatch).stream();
}


const Foam::entry* Foam::dictionary::lookupScopedEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    if ((keyword[0] == ':' || keyword[0] == '^'))
    {
        // Go up to top level
        const dictionary* dictPtr = this;
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }

        return dictPtr->lookupScopedSubEntryPtr
        (
            keyword.substr(1),
            false,
            patternMatch
        );
    }

    return lookupScopedSubEntryPtr
    (
        keyword,
        recursive,
        patternMatch
    );
}


bool Foam::dictionary::substituteScopedKeyword
(
    const word& keyword,
    bool mergeEntry
)
{
    // Drop first character of keyword to get the var-name, already validated.
    const word varName(keyword.substr(1), false);

    // Lookup the variable name in the given dictionary
    const entry* ePtr = lookupScopedEntryPtr(varName, true, true);

    // If defined insert its entries into this dictionary
    if (ePtr != nullptr)
    {
        const dictionary& addDict = ePtr->dict();

        forAllConstIter(parent_type, addDict, iter)
        {
            add(iter(), mergeEntry);
        }

        return true;
    }

    return false;
}


bool Foam::dictionary::isDict(const word& keyword) const
{
    // Find non-recursive with patterns
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return entryPtr->isDict();
    }
    else
    {
        return false;
    }
}


const Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return &entryPtr->dict();
    }
    else
    {
        return nullptr;
    }
}


Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword)
{
    entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return &entryPtr->dict();
    }
    else
    {
        return nullptr;
    }
}


const Foam::dictionary& Foam::dictionary::subDict(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }
    return entryPtr->dict();
}


Foam::dictionary& Foam::dictionary::subDict(const word& keyword)
{
    entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }
    return entryPtr->dict();
}


Foam::dictionary Foam::dictionary::subOrEmptyDict
(
    const word& keyword,
    const bool mustRead
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        if (mustRead)
        {
            FatalIOErrorInFunction
            (
                *this
            )   << "keyword " << keyword << " is undefined in dictionary "
                << name()
                << exit(FatalIOError);
            return entryPtr->dict();
        }
        else
        {
            return dictionary(*this, dictionary(name() + '.' + keyword));
        }
    }
    else
    {
        return entryPtr->dict();
    }
}


const Foam::dictionary& Foam::dictionary::optionalSubDict
(
    const word& keyword
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return entryPtr->dict();
    }
    else
    {
        return *this;
    }
}


Foam::wordList Foam::dictionary::toc() const
{
    wordList keys(size());

    label nKeys = 0;
    forAllConstIter(parent_type, *this, iter)
    {
        keys[nKeys++] = iter().keyword();
    }

    return keys;
}


Foam::wordList Foam::dictionary::sortedToc() const
{
    return hashedEntries_.sortedToc();
}


Foam::List<Foam::keyType> Foam::dictionary::keys(bool patterns) const
{
    List<keyType> keys(size());

    label nKeys = 0;
    forAllConstIter(parent_type, *this, iter)
    {
        if (iter().keyword().isPattern() ? patterns : !patterns)
        {
            keys[nKeys++] = iter().keyword();
        }
    }
    keys.setSize(nKeys);

    return keys;
}


bool Foam::dictionary::add(entry* entryPtr, bool mergeEntry)
{
    auto iter = hashedEntries_.find(entryPtr->keyword());

    if (mergeEntry && iter.found())
    {
        // Merge dictionary with dictionary
        if (iter()->isDict() && entryPtr->isDict())
        {
            iter()->dict().merge(entryPtr->dict());
            delete entryPtr;

            return true;
        }


        // Replace existing dictionary with entry or vice versa
        parent_type::replace(iter(), entryPtr);
        delete iter();
        hashedEntries_.erase(iter);

        if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
        {
            entryPtr->name() = name() + '.' + entryPtr->keyword();

            if (entryPtr->keyword().isPattern())
            {
                patterns_.insert(entryPtr);
                regexps_.insert
                (
                    autoPtr<regExp>(new regExp(entryPtr->keyword()))
                );
            }

            return true;
        }
        else
        {
            IOWarningInFunction((*this))
                << "problem replacing entry "<< entryPtr->keyword()
                << " in dictionary " << name() << endl;

            parent_type::remove(entryPtr);
            delete entryPtr;
            return false;
        }
    }


    if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
    {
        entryPtr->name() = name() + '.' + entryPtr->keyword();
        parent_type::append(entryPtr);

        if (entryPtr->keyword().isPattern())
        {
            patterns_.insert(entryPtr);
            regexps_.insert
            (
                autoPtr<regExp>(new regExp(entryPtr->keyword()))
            );
        }

        return true;
    }
    else
    {
        IOWarningInFunction((*this))
            << "attempt to add entry "<< entryPtr->keyword()
            << " which already exists in dictionary " << name()
            << endl;

        delete entryPtr;
        return false;
    }
}


void Foam::dictionary::add(const entry& e, bool mergeEntry)
{
    add(e.clone(*this).ptr(), mergeEntry);
}


void Foam::dictionary::add(const keyType& k, const word& w, bool overwrite)
{
    add(new primitiveEntry(k, token(w)), overwrite);
}


void Foam::dictionary::add
(
    const keyType& k,
    const Foam::string& s,
    bool overwrite
)
{
    add(new primitiveEntry(k, token(s)), overwrite);
}


void Foam::dictionary::add(const keyType& k, const label l, bool overwrite)
{
    add(new primitiveEntry(k, token(l)), overwrite);
}


void Foam::dictionary::add(const keyType& k, const scalar s, bool overwrite)
{
    add(new primitiveEntry(k, token(s)), overwrite);
}


void Foam::dictionary::add
(
    const keyType& k,
    const dictionary& d,
    bool mergeEntry
)
{
    add(new dictionaryEntry(k, *this, d), mergeEntry);
}


void Foam::dictionary::set(entry* entryPtr)
{
    // Find non-recursive with patterns
    entry* existingPtr = lookupEntryPtr(entryPtr->keyword(), false, true);

    // Clear dictionary so merge acts like overwrite
    if (existingPtr && existingPtr->isDict())
    {
        existingPtr->dict().clear();
    }
    add(entryPtr, true);
}


void Foam::dictionary::set(const entry& e)
{
    set(e.clone(*this).ptr());
}


void Foam::dictionary::set(const keyType& k, const dictionary& d)
{
    set(new dictionaryEntry(k, *this, d));
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
    bool forceOverwrite
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
        if (forceOverwrite)
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


bool Foam::dictionary::merge(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalIOErrorInFunction(*this)
            << "attempted merge to self for dictionary " << name()
            << abort(FatalIOError);
    }

    bool changed = false;

    forAllConstIter(parent_type, dict, iter)
    {
        auto fnd = hashedEntries_.find(iter().keyword());

        if (fnd.found())
        {
            // Recursively merge sub-dictionaries
            // TODO: merge without copying
            if (fnd()->isDict() && iter().isDict())
            {
                if (fnd()->dict().merge(iter().dict()))
                {
                    changed = true;
                }
            }
            else
            {
                add(iter().clone(*this).ptr(), true);
                changed = true;
            }
        }
        else
        {
            // Not found - just add
            add(iter().clone(*this).ptr());
            changed = true;
        }
    }

    return changed;
}


void Foam::dictionary::clear()
{
    parent_type::clear();
    hashedEntries_.clear();
    patterns_.clear();
    regexps_.clear();
}


void Foam::dictionary::transfer(dictionary& dict)
{
    // Changing parents probably doesn't make much sense,
    // but what about the names?
    name() = dict.name();

    parent_type::transfer(dict);
    hashedEntries_.transfer(dict.hashedEntries_);
    patterns_.transfer(dict.patterns_);
    regexps_.transfer(dict.regexps_);
}


Foam::Xfer<Foam::dictionary> Foam::dictionary::xfer()
{
    return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::ITstream& Foam::dictionary::operator[](const word& keyword) const
{
    return lookup(keyword);
}


void Foam::dictionary::operator=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    name() = rhs.name();
    clear();

    // Create clones of the entries in the given dictionary
    // resetting the parentDict to this dictionary

    forAllConstIter(parent_type, rhs, iter)
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator+=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted addition assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(parent_type, rhs, iter)
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator|=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(parent_type, rhs, iter)
    {
        if (!found(iter().keyword()))
        {
            add(iter().clone(*this).ptr());
        }
    }
}


void Foam::dictionary::operator<<=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(parent_type, rhs, iter)
    {
        set(iter().clone(*this).ptr());
    }
}


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

Foam::dictionary Foam::operator+
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary sum(dict1);
    sum += dict2;
    return sum;
}


Foam::dictionary Foam::operator|
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary sum(dict1);
    sum |= dict2;
    return sum;
}


// ************************************************************************* //
