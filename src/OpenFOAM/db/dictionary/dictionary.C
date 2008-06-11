/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "dictionary.H"
#include "primitiveEntry.H"
#include "dictionaryEntry.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

defineTypeNameAndDebug(Foam::dictionary, 0);

const Foam::dictionary Foam::dictionary::null;

#define DICTIONARY_INPLACE_MERGE

// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

bool Foam::dictionary::add(entry* ePtr, bool mergeEntry)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(ePtr->keyword());

    if (mergeEntry && iter != hashedEntries_.end())
    {
        // merge dictionary with dictionary
        if (iter()->isDict() && ePtr->isDict())
        {
            iter()->dict().merge(ePtr->dict());
            delete ePtr;

            return true;
        }
        else
        {
            // replace existing dictionary with entry or vice versa
#ifdef DICTIONARY_INPLACE_MERGE
            IDLList<entry>::replace(iter(), ePtr);
            delete iter();
            hashedEntries_.erase(iter);

            if (hashedEntries_.insert(ePtr->keyword(), ePtr))
            {
                ePtr->name() = name_ + "::" + ePtr->keyword();
                return true;
            }
            else
            {
                IOWarningIn("dictionary::add(entry* ePtr)", (*this))
                    << "problem replacing entry "<< ePtr->keyword()
                    << " in dictionary " << name()
                    << endl;

                IDLList<entry>::remove(ePtr);
                delete ePtr;
                return false;
            }
#else
            remove(ePtr->keyword());
#endif
        }
    }

    if (hashedEntries_.insert(ePtr->keyword(), ePtr))
    {
        ePtr->name() = name_ + "::" + ePtr->keyword();
        IDLList<entry>::append(ePtr);

        return true;
    }
    else
    {
        IOWarningIn("dictionary::add(entry* ePtr)", (*this))
            << "attempt to add entry "<< ePtr->keyword()
            << " which already exists in dictionary " << name()
            << endl;

        delete ePtr;
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionary::dictionary()
:
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    const dictionary& dict
)
:
    IDLList<entry>(dict, *this),
    name_(dict.name()),
    parent_(parentDict)
{
    for
    (
        IDLList<entry>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        hashedEntries_.insert(iter().keyword(), &iter());
    }
}


Foam::dictionary::dictionary
(
    const dictionary& dict
)
:
    IDLList<entry>(dict, *this),
    name_(dict.name()),
    parent_(dictionary::null)
{
    for
    (
        IDLList<entry>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        hashedEntries_.insert(iter().keyword(), &iter());
    }
}


Foam::autoPtr<Foam::dictionary> Foam::dictionary::clone() const
{
    return autoPtr<dictionary>(new dictionary(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dictionary::~dictionary()
{
    // cerr<< "~dictionary() " << name() << " " << long(this) << std::endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


bool Foam::dictionary::found(const word& keyword, bool recursive) const
{
    if (hashedEntries_.found(keyword))
    {
        return true;
    }
    else if (recursive && &parent_ != &dictionary::null)
    {
        return parent_.found(keyword, recursive);
    }
    else
    {
        return false;
    }
}


const Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive
) const
{
    HashTable<entry*>::const_iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        if (recursive && &parent_ != &dictionary::null)
        {
            return parent_.lookupEntryPtr(keyword, recursive);
        }
        else
        {
            return NULL;
        }
    }

    return iter();
}


Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive
)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        if (recursive && &parent_ != &dictionary::null)
        {
            return const_cast<dictionary&>(parent_).lookupEntryPtr
            (
                keyword,
                recursive
            );
        }
        else
        {
            return NULL;
        }
    }

    return iter();
}


const Foam::entry& Foam::dictionary::lookupEntry
(
    const word& keyword,
    bool recursive
) const
{
    const entry* ePtr = lookupEntryPtr(keyword, recursive);

    if (ePtr == NULL)
    {
        // If keyword not found print error message ...
        FatalIOErrorIn
        (
            "dictionary::lookupEntry(const word& keyword) const",
            *this
        )   << " keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return *ePtr;
}


Foam::ITstream& Foam::dictionary::lookup
(
    const word& keyword,
    bool recursive
) const
{
    return lookupEntry(keyword, recursive).stream();
}


bool Foam::dictionary::isDict(const word& keyword) const
{
    if (const entry* entryPtr = lookupEntryPtr(keyword))
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
    if (const entry* entryPtr = lookupEntryPtr(keyword))
    {
        return &entryPtr->dict();
    }
    else
    {
        return NULL;
    }
}


const Foam::dictionary& Foam::dictionary::subDict(const word& keyword) const
{
    if (const entry* entryPtr = lookupEntryPtr(keyword))
    {
        return entryPtr->dict();
    }
    else
    {
        // If keyword not found print error message ...
        FatalIOErrorIn
        (
            "dictionary::subDict(const word& keyword) const",
            *this
        )   << " keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);

        return entryPtr->dict();
    }
}


Foam::dictionary& Foam::dictionary::subDict(const word& keyword)
{
    if (entry* entryPtr = lookupEntryPtr(keyword))
    {
        return entryPtr->dict();
    }
    else
    {
        // If keyword not found print error message ...
        FatalIOErrorIn
        (
            "dictionary::subDict(const word& keyword) const",
            *this
        )   << " keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);

        return entryPtr->dict();
    }
}


Foam::wordList Foam::dictionary::toc() const
{
    wordList keywords(size());

    label i = 0;
    for
    (
        IDLList<entry>::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        keywords[i++] = iter().keyword();
    }

    return keywords;
}


void Foam::dictionary::add(const entry& e)
{
    add(e.clone(*this).ptr());
}

void Foam::dictionary::add(const word& keyword, const token& t)
{
    add(new primitiveEntry(keyword, t));
}

void Foam::dictionary::add(const word& keyword, const word& w)
{
    add(new primitiveEntry(keyword, token(w)));
}

void Foam::dictionary::add(const word& keyword, const Foam::string& s)
{
    add(new primitiveEntry(keyword, token(s)));
}

void Foam::dictionary::add(const word& keyword, const label l)
{
    add(new primitiveEntry(keyword, token(l)));
}

void Foam::dictionary::add(const word& keyword, const scalar s)
{
    add(new primitiveEntry(keyword, token(s)));
}

void Foam::dictionary::add(const word& keyword, const ITstream& tokens)
{
    add(new primitiveEntry(keyword, tokens));
}

void Foam::dictionary::add(const word& keyword, const tokenList& tokens)
{
    add(new primitiveEntry(keyword, tokens));
}

void Foam::dictionary::add(const word& keyword, const dictionary& dict)
{
    add(new dictionaryEntry(keyword, *this, dict));
}


bool Foam::dictionary::remove(const word& Keyword)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(Keyword);

    if (iter != hashedEntries_.end())
    {
        IDLList<entry>::remove(iter());
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
    const word& oldKeyword,
    const word& newKeyword,
    bool forceOverwrite
)
{
    // no change
    if (oldKeyword == newKeyword)
    {
        return false;
    }

    HashTable<entry*>::iterator iter = hashedEntries_.find(oldKeyword);

    // oldKeyword not found - do nothing
    if (iter == hashedEntries_.end())
    {
        return false;
    }

    HashTable<entry*>::iterator iter2 = hashedEntries_.find(newKeyword);

    // newKeyword already exists
    if (iter2 != hashedEntries_.end())
    {
        if (forceOverwrite)
        {
            IDLList<entry>::replace(iter2(), iter());
            delete iter2();
            hashedEntries_.erase(iter2);
        }
        else
        {
            // could issue warning if desired
            return false;
        }
    }

    // change name and HashTable, but leave DL-List untouched
    iter()->keyword() = newKeyword;
    iter()->name() = name_ + "::" + newKeyword;
    hashedEntries_.erase(oldKeyword);
    hashedEntries_.insert(newKeyword, iter());

    return true;
}


bool Foam::dictionary::merge(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::merge(const dictionary&)")
            << "attempted merge to self for dictionary " << name()
            << abort(FatalError);
    }

    bool changed = false;

    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        const word& keyword = iter().keyword();

        HashTable<entry*>::iterator iter2 = hashedEntries_.find(keyword);

        if (iter2 != hashedEntries_.end())
        {
            // Recursively merge sub-dictionaries
            if (iter2()->isDict() && iter().isDict())
            {
                // without copying and without remove/add?
                // this certainly looks ugly and doesn't necessarily
                // retain the original sort order (perhaps nice to have)
                if
                (
                    iter2()->dict().merge(iter().dict())
                )
                {
                    changed = true;
                }
            }
            else
            {
#ifdef DICTIONARY_INPLACE_MERGE
                add(iter().clone(*this).ptr(), true);
#else
                remove(keyword);
                add(iter().clone(*this)());
#endif
                changed = true;
            }
        }
        else
        {
            // not found - just add
            add(iter().clone(*this)());
            changed = true;
        }
    }

    return changed;
}


void Foam::dictionary::clear()
{
    IDLList<entry>::clear();
    hashedEntries_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::ITstream& Foam::dictionary::operator[](const word& keyword) const
{
    return lookup(keyword);
}


void Foam::dictionary::operator=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator=(const dictionary&)")
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    // Clear the current entries
    IDLList<entry>::clear();

    // Create clones of the entries in the given dictionary
    // resetting the parentDict to this dictionary
    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        IDLList<entry>::append(iter().clone(*this).ptr());
    }

    name_ = dict.name();

    hashedEntries_.clear();

    for
    (
        IDLList<entry>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        hashedEntries_.insert(iter().keyword(), &iter());
    }
}


void Foam::dictionary::operator+=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator+=(const dictionary&)")
            << "attempted addition assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        add(iter().clone(*this)());
    }
}


void Foam::dictionary::operator|=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator|=(const dictionary&)")
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        if (!found(iter().keyword()))
        {
            add(iter().clone(*this)());
        }
    }
}


void Foam::dictionary::operator<<=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator<<=(const dictionary&)")
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    for
    (
        IDLList<entry>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        remove(iter().keyword());
        add(iter().clone(*this)());
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
