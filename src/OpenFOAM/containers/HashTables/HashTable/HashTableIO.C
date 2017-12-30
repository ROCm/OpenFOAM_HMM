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

#include "HashTable.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(Istream& is, const label size)
:
    HashTableCore(),
    size_(0),
    capacity_(HashTableCore::canonicalSize(size)),
    table_(nullptr)
{
    if (capacity_)
    {
        table_ = new node_type*[capacity_];

        for (label i=0; i < capacity_; ++i)
        {
            table_[i] = nullptr;
        }
    }

    operator>>(is, *this);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Ostream& Foam::HashTable<T, Key, Hash>::printInfo(Ostream& os) const
{
    label used = 0;
    label maxChain = 0;
    unsigned avgChain = 0;

    for (label i=0; i < capacity_; ++i)
    {
        label count = 0;
        for (node_type* ep = table_[i]; ep; ep = ep->next_)
        {
            ++count;
        }

        if (count)
        {
            ++used;
            avgChain += count;

            if (maxChain < count)
            {
                maxChain = count;
            }
        }
    }

    os  << "HashTable<T,Key,Hash>"
        << " elements:" << size() << " slots:" << used << "/" << capacity_
        << " chaining(avg/max):" << (used ? (float(avgChain)/used) : 0)
        << "/" << maxChain << endl;

    return os;
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::HashTable<T, Key, Hash>::writeKeys
(
    Ostream& os,
    const label shortListLen
) const
{
    // Similar to UList::writeList version except the following:
    // - the keys can never be uniform
    // - never write in binary

    label i = this->size();

    if (i <= 1 || !shortListLen || (i <= shortListLen))
    {
        // Write size and start delimiter
        os << i << token::BEGIN_LIST;

        i = 0;
        for (const_iterator iter = this->cbegin(); iter != this->cend(); ++iter)
        {
            if (i++) os << token::SPACE;
            os << iter.key();
        }

        os << token::END_LIST;  // End delimiter
    }
    else
    {
        // Write size and start delimiter
        os << nl << i << nl << token::BEGIN_LIST << nl;

        for (const_iterator iter = this->cbegin(); iter != this->cend(); ++iter)
        {
            os << iter.key() << nl;
        }

        os << token::END_LIST << nl;  // End delimiter
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    HashTable<T, Key, Hash>& L
)
{
    is.fatalCheck(FUNCTION_NAME);

    // Anull existing table
    L.clear();

    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck
    (
        "operator>>(Istream&, HashTable&) : "
        "reading first token"
    );

    if (firstToken.isLabel())
    {
        const label s = firstToken.labelToken();

        // Read beginning of contents
        const char delimiter = is.readBeginList("HashTable");

        if (s)
        {
            if (2*s > L.capacity_)
            {
                L.resize(2*s);
            }

            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<s; ++i)
                {
                    Key key;
                    is >> key;
                    L.insert(key, pTraits<T>(is));

                    is.fatalCheck
                    (
                        "operator>>(Istream&, HashTable&) : "
                        "reading entry"
                    );
                }
            }
            else
            {
                FatalIOErrorInFunction
                (
                    is
                )   << "incorrect first token, '(', found " << firstToken.info()
                    << exit(FatalIOError);
            }
        }

        // Read end of contents
        is.readEndList("HashTable");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction
            (
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);

            Key key;
            is >> key;
            L.insert(key, pTraits<T>(is));

            is.fatalCheck
            (
                "operator>>(Istream&, HashTable&) : "
                "reading entry"
            );

            is >> lastToken;
        }
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    is.fatalCheck(FUNCTION_NAME);

    return is;
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const HashTable<T, Key, Hash>& tbl
)
{
    const label sz = tbl.size();

    if (sz)
    {
        // Size and start list delimiter
        os << nl << sz << nl << token::BEGIN_LIST << nl;

        // Contents
        for (auto iter = tbl.cbegin(); iter != tbl.cend(); ++iter)
        {
            os << iter.key() << token::SPACE << iter.object() << nl;
        }

        os << token::END_LIST;    // End list delimiter
    }
    else
    {
        // Empty hash table
        os << sz << token::BEGIN_LIST << token::END_LIST;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
