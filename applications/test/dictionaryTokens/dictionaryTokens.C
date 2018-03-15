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

#include "dictionaryTokens.H"
#include "IOstream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::token Foam::dictionaryTokens::keywordToken(const entry& e)
{
    const keyType& k = e.keyword();

    if (k.empty())
    {
        return token::undefinedToken;
    }
    if (k.isPattern())
    {
        return token(static_cast<string>(k)); // quoted
    }
    else
    {
        return token(static_cast<word>(k)); // unquoted
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::dictionaryTokens::setIterator() const
{
    primIter_.clear();
    dictIter_.clear();

    if (entryIter_ != dict_.cend())
    {
        const entry& base = *entryIter_;

        if (isA<primitiveEntry>(base))
        {
            primIter_.reset
            (
                new dictionaryTokens::primitive_iterator
                (
                    dynamicCast<const primitiveEntry&>(base)
                )
            );

            return true;
        }
        else if (isA<dictionaryListEntry>(base))
        {
            // Must check for isA<dictionaryListEntry> before checking
            // for isA<dictionaryEntry> !

            dictIter_.reset
            (
                new dictionaryTokens::dictionary_iterator
                (
                    dynamicCast<const dictionaryListEntry&>(base)
                )
            );

            return true;
        }
        else if (isA<dictionaryEntry>(base))
        {
            dictIter_.reset
            (
                new dictionaryTokens::dictionary_iterator
                (
                    dynamicCast<const dictionaryEntry&>(base)
                )
            );

            return true;
        }

    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryTokens::dictionaryTokens(const dictionary& dict)
:
    dict_(dict),
    entryIter_(dict_.cbegin()),
    primIter_(nullptr),
    dictIter_(nullptr)
{
    rewind();
}


Foam::dictionaryTokens::primitive_iterator::primitive_iterator
(
    const primitiveEntry& e
)
:
    tokensPtr_(&(static_cast<const tokenList&>(e))),
    key_(dictionaryTokens::keywordToken(e)),
    end_(token::punctuationToken::END_STATEMENT),
    pos_((key_.good() ? -1 : 0))
{}


Foam::dictionaryTokens::dictionary_iterator::dictionary_iterator
(
    const dictionaryEntry& e
)
:
    key_(dictionaryTokens::keywordToken(e)),
    lbrace_(token::punctuationToken::BEGIN_BLOCK),
    rbrace_(token::punctuationToken::END_BLOCK),
    state_(key_.good() ? states::KEY : states::OPEN),
    dictTokens_(e.dict())
{}


Foam::dictionaryTokens::dictionary_iterator::dictionary_iterator
(
    const dictionaryListEntry& e
)
:
    key_(e.size()),
    lbrace_(token::punctuationToken::BEGIN_LIST),
    rbrace_(token::punctuationToken::END_LIST),
    state_(states::KEY),
    dictTokens_(e.dict())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dictionaryTokens::good() const
{
    return
    (
        entryIter_ != dict_.cend()
     &&
        (
            (primIter_.valid() && primIter_().good())
         || (dictIter_.valid() && dictIter_().good())
        )
    );
}


bool Foam::dictionaryTokens::primitive_iterator::good() const
{
    return (tokensPtr_ && pos_ <= tokensPtr_->size());
}


bool Foam::dictionaryTokens::dictionary_iterator::good() const
{
    return (state_ != states::END);
}


const Foam::token& Foam::dictionaryTokens::operator*() const
{
    if (good())
    {
        if (primIter_.valid()) return *(primIter_());
        if (dictIter_.valid()) return *(dictIter_());
    }

    return token::undefinedToken;
}


const Foam::token&
Foam::dictionaryTokens::primitive_iterator::operator*() const
{
    if (good())
    {
        if (pos_ == -1)
        {
            return key_;
        }
        else if (pos_ >= tokensPtr_->size())
        {
            return end_; // The trailing ';'
        }

        return tokensPtr_->operator[](pos_);
    }

    return token::undefinedToken;
}


const Foam::token&
Foam::dictionaryTokens::dictionary_iterator::operator*() const
{
    if (good())
    {
        if (state_ == states::KEY)
        {
            return key_;  // keyword
        }
        if (state_ == states::OPEN)
        {
            return lbrace_;  // Opening '{'
        }
        if (state_ == states::CONTENT)
        {
            return *(dictTokens_);
        }
        if (state_ == states::CLOSE)
        {
            return rbrace_;  // Closing '}'
        }
    }

    return token::undefinedToken;
}


bool Foam::dictionaryTokens::operator++()
{
    bool ok = good();

    if (ok)
    {
        if (primIter_.valid()) ok = ++(primIter_());
        if (dictIter_.valid()) ok = ++(dictIter_());

        if (!ok)
        {
            ++entryIter_; // Next entry
            setIterator();
        }
    }

    return ok;
}


bool Foam::dictionaryTokens::primitive_iterator::operator++()
{
    // Advance good iterators.
    //
    // Going beyond trailing ';' makes it into an end iterator

    if (tokensPtr_ && (++pos_ > tokensPtr_->size()))
    {
        tokensPtr_ = nullptr;
        return false;
    }

    return this->good();
}


bool Foam::dictionaryTokens::dictionary_iterator::operator++()
{
    if
    (
        state_ == states::KEY
     || state_ == states::OPEN
     || state_ == states::CLOSE
    )
    {
        ++state_;
    }
    else if (state_ == states::CONTENT && !(++dictTokens_))
    {
        ++state_;
    }

    return good();
}


void Foam::dictionaryTokens::rewind()
{
    entryIter_ = dict_.cbegin();
    setIterator();
}


// ************************************************************************* //
