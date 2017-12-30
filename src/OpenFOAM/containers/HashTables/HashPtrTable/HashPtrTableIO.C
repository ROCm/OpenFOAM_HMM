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

#include "HashPtrTable.H"
#include "Istream.H"
#include "Ostream.H"
#include "INew.H"
#include "dictionary.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class INew>
void Foam::HashPtrTable<T, Key, Hash>::read(Istream& is, const INew& inewt)
{
    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck
    (
        "HashPtrTable::read(Istream&, const INew&) : "
        "reading first token"
    );

    if (firstToken.isLabel())
    {
        const label s = firstToken.labelToken();

        // Read beginning of contents
        const char delimiter = is.readBeginList("HashPtrTable");

        if (s)
        {
            if (2*s > this->capacity())
            {
                this->resize(2*s);
            }

            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<s; ++i)
                {
                    Key key;
                    is >> key;
                    this->insert(key, inewt(key, is).ptr());

                    is.fatalCheck
                    (
                        "HashPtrTable::read(Istream&, const INew&) : "
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
        is.readEndList("HashPtrTable");
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
            this->insert(key, inewt(key, is).ptr());

            is.fatalCheck
            (
                "HashPtrTable::read(Istream&, const INew&) : "
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
}


template<class T, class Key, class Hash>
template<class INew>
void Foam::HashPtrTable<T, Key, Hash>::read
(
    const dictionary& dict,
    const INew& inewt
)
{
    forAllConstIter(dictionary, dict, iter)
    {
        const word& k = iter().keyword();

        this->insert(k, inewt(dict.subDict(k)).ptr());
    }
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::write(Ostream& os) const
{
    for (const_iterator iter = this->cbegin(); iter != this->cend(); ++iter)
    {
        const T* ptr = iter.object();
        if (ptr)
        {
            ptr->write(os);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class INew>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(Istream& is, const INew& inewt)
{
    this->read(is, inewt);
}


template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(Istream& is)
{
    this->read(is, INew<T>());
}


template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(const dictionary& dict)
{
    this->read(dict, INew<T>());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Istream& Foam::operator>>(Istream& is, HashPtrTable<T, Key, Hash>& tbl)
{
    tbl.clear();
    tbl.read(is, INew<T>());

    return is;
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const HashPtrTable<T, Key, Hash>& tbl
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
            const T* ptr = iter.object();

            os << iter.key();
            if (ptr)
            {
                os << token::SPACE << *ptr;
            }
            os << nl;
        }
        os << token::END_LIST; // End list delimiter
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
