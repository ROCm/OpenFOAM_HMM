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

    return *this;
}


Foam::label Foam::dictionary::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }

    return -1;
}


Foam::label Foam::dictionary::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }

    return -1;
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
    // Serialize dictionary entries into a string
    OStringStream os;

    // Process entries
    forAllConstIter(parent_type, *this, iter)
    {
        os << *iter;
    }

    // String re-parsed as a list of tokens
    return ITstream::parse(os.str());
}


bool Foam::dictionary::found
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    return csearch(keyword, recursive, patternMatch).found();
}


const Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    return csearch(keyword, recursive, patternMatch).ptr();
}


Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
)
{
    return search(keyword, recursive, patternMatch).ptr();
}


const Foam::entry& Foam::dictionary::lookupEntry
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    const const_searcher finder(csearch(keyword, recursive, patternMatch));

    if (!finder.found())
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return finder.ref();
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
    return csearchScoped(keyword, recursive, patternMatch).ptr();
}


bool Foam::dictionary::substituteKeyword(const word& keyword, bool mergeEntry)
{
    if (keyword.size() < 2)
    {
        return false;
    }

    // Drop leading '$' to get the var-name, already validated as word.
    const word varName(keyword.substr(1), false);

    // Lookup the variable name in the given dictionary
    const const_searcher finder(csearch(varName, true, true));

    // If defined insert its entries into this dictionary
    if (finder.found())
    {
        const dictionary& addDict = finder.dict();

        forAllConstIters(addDict, iter)
        {
            add(iter(), mergeEntry);
        }

        return true;
    }

    return false;
}


bool Foam::dictionary::substituteScopedKeyword
(
    const word& keyword,
    bool mergeEntry
)
{
    if (keyword.size() < 2)
    {
        return false;
    }

    // Drop leading '$' to get the var-name, already validated as word.
    const word varName(keyword.substr(1), false);

    // Lookup the variable name in the given dictionary
    const const_searcher finder(csearchScoped(varName, true, true));

    // If defined insert its entries into this dictionary
    if (finder.found())
    {
        const dictionary& addDict = finder.dict();

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
    return csearch(keyword, false, true).isDict();
}


const Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword) const
{
    // Find non-recursive with patterns
    return csearch(keyword, false, true).dictPtr();
}


Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword)
{
    // Find non-recursive with patterns
    return search(keyword, false, true).dictPtr();
}


const Foam::dictionary& Foam::dictionary::subDict(const word& keyword) const
{
    // Find non-recursive with patterns
    const const_searcher finder(csearch(keyword, false, true));

    if (!finder.found())
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return finder.dict();
}


Foam::dictionary& Foam::dictionary::subDict(const word& keyword)
{
    // Find non-recursive with patterns
    searcher finder = search(keyword, false, true);

    if (!finder.found())
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return finder.dict();
}


Foam::dictionary Foam::dictionary::subOrEmptyDict
(
    const word& keyword,
    const bool mustRead
) const
{
    // Find non-recursive with patterns
    const const_searcher finder(csearch(keyword, false, true));

    if (finder.isDict())
    {
        // Found and a sub-dictionary
        return finder.dict();
    }

    if (mustRead)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword
            << " is not a sub-dictionary in dictionary "
            << name()
            << exit(FatalIOError);
    }

    if (finder.found())
    {
        IOWarningInFunction((*this))
            << "keyword " << keyword
            << " found but not a sub-dictionary in dictionary "
            << name() << endl;
    }

    return dictionary(*this, dictionary(name() + '.' + keyword));
}


const Foam::dictionary& Foam::dictionary::optionalSubDict
(
    const word& keyword
) const
{
    const const_searcher finder(csearch(keyword, false, true));

    if (finder.isDict())
    {
        // Found and a sub-dictionary
        return finder.dict();
    }

    if (finder.found())
    {
        IOWarningInFunction((*this))
            << "keyword " << keyword
            << " found but not a sub-dictionary in dictionary "
            << name() << endl;
    }

    return *this;
}


Foam::wordList Foam::dictionary::toc() const
{
    wordList keys(size());

    label n = 0;
    forAllConstIters(*this, iter)
    {
        keys[n++] = iter().keyword();
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

    label n = 0;
    forAllConstIters(*this, iter)
    {
        if (iter().keyword().isPattern() ? patterns : !patterns)
        {
            keys[n++] = iter().keyword();
        }
    }
    keys.setSize(n);

    return keys;
}


Foam::entry* Foam::dictionary::add(entry* entryPtr, bool mergeEntry)
{
    if (!entryPtr)
    {
        return nullptr;
    }

    auto iter = hashedEntries_.find(entryPtr->keyword());

    if (mergeEntry && iter.found())
    {
        // Merge dictionary with dictionary
        if (iter()->isDict() && entryPtr->isDict())
        {
            iter()->dict().merge(entryPtr->dict());

            delete entryPtr;
            return iter();   // pointer to existing dictionary
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

            return entryPtr;  // now an entry in the dictionary
        }


        IOWarningInFunction((*this))
            << "problem replacing entry "<< entryPtr->keyword()
            << " in dictionary " << name() << endl;

        parent_type::remove(entryPtr);

        delete entryPtr;
        return nullptr;
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

        return entryPtr;  // now an entry in the dictionary
    }


    IOWarningInFunction((*this))
        << "attempt to add entry " << entryPtr->keyword()
        << " which already exists in dictionary " << name()
        << endl;

    delete entryPtr;
    return nullptr;
}


Foam::entry* Foam::dictionary::add(const entry& e, bool mergeEntry)
{
    return add(e.clone(*this).ptr(), mergeEntry);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const word& v,
    bool overwrite
)
{
    return add(new primitiveEntry(k, token(v)), overwrite);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const Foam::string& v,
    bool overwrite
)
{
    return add(new primitiveEntry(k, token(v)), overwrite);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const label v,
    bool overwrite
)
{
    return add(new primitiveEntry(k, token(v)), overwrite);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const scalar v,
    bool overwrite
)
{
    return add(new primitiveEntry(k, token(v)), overwrite);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const dictionary& v,
    bool mergeEntry
)
{
    return add(new dictionaryEntry(k, *this, v), mergeEntry);
}


Foam::entry* Foam::dictionary::set(entry* entryPtr)
{
    if (!entryPtr)
    {
        return nullptr;
    }

    // Find non-recursive with patterns
    searcher finder(search(entryPtr->keyword(), false, true));

    // Clear dictionary so merge acts like overwrite
    if (finder.isDict())
    {
        finder.dict().clear();
    }

    return add(entryPtr, true);
}


Foam::entry* Foam::dictionary::set(const entry& e)
{
    return set(e.clone(*this).ptr());
}


Foam::entry* Foam::dictionary::set(const keyType& k, const dictionary& v)
{
    return set(new dictionaryEntry(k, *this, v));
}


bool Foam::dictionary::merge(const dictionary& dict)
{
    if (this == &dict)
    {
        FatalIOErrorInFunction(*this)
            << "attempted merge to self for dictionary " << name()
            << abort(FatalIOError);
    }

    bool changed = false;

    forAllConstIters(dict, iter)
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

    forAllConstIters(rhs, iter)
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator+=(const dictionary& rhs)
{
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted addition assignment to self for dictionary "
            << name()
            << abort(FatalIOError);
    }

    forAllConstIters(rhs, iter)
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator|=(const dictionary& rhs)
{
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary "
            << name()
            << abort(FatalIOError);
    }

    forAllConstIters(rhs, iter)
    {
        if (!found(iter().keyword()))
        {
            add(iter().clone(*this).ptr());
        }
    }
}


void Foam::dictionary::operator<<=(const dictionary& rhs)
{
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary "
            << name()
            << abort(FatalIOError);
    }

    forAllConstIters(rhs, iter)
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
