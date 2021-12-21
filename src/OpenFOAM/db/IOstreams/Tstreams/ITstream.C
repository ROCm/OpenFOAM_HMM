/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "error.H"
#include "ITstream.H"
#include "StringStream.H"
#include "UIListStream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Failsafe read-access.
// Return the token at location, or undefinedToken.
inline static const token& peekTokenAt
(
    const UList<token>& list,
    const label i
)
{
    return
    (
        i >= 0 && i < list.size()
      ? list[i]
      : token::undefinedToken
    );
}


// Convert input sequence into a list of tokens.
// Return the number of tokens in the resulting list.
static label parseStream(ISstream& is, tokenList& tokens)
{
    label nTok = 0;

    tokens.clear();
    tokens.resize(64, token());

    token tok;
    while (!is.read(tok).bad() && tok.good())
    {
        tokens.newElmt(nTok++) = std::move(tok);
    }

    tokens.resize(nTok);

    return nTok;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tokenList Foam::ITstream::parse
(
    const UList<char>& input,
    IOstreamOption streamOpt
)
{
    UIListStream is(input, streamOpt);

    tokenList tokens;
    parseStream(is, tokens);
    return tokens;
}


Foam::tokenList Foam::ITstream::parse
(
    const std::string& input,
    IOstreamOption streamOpt
)
{
    UIListStream is(input.data(), input.length(), streamOpt);

    tokenList tokens;
    parseStream(is, tokens);
    return tokens;
}


Foam::tokenList Foam::ITstream::parse
(
    const char* input,
    IOstreamOption streamOpt
)
{
    UIListStream is(input, strlen(input), streamOpt);

    tokenList tokens;
    parseStream(is, tokens);
    return tokens;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ITstream::reserveCapacity
(
    const label nElem,
    const bool lazy
)
{
    if (lazy)
    {
        // Reserve - leave excess capacity for further appends

        label n = tokenList::size();

        if (nElem > n)
        {
            if (!n) n = 1;  // Avoid dead-lock when starting from zero-sized

            do
            {
                n *= 2;
            }
            while (nElem >= n);

            tokenList::resize(n);
        }
    }
    else
    {
        // Strict capacity
        tokenList::resize(nElem);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ITstream::ITstream(const ITstream& is)
:
    Istream(static_cast<IOstreamOption>(is)),
    tokenList(is),
    name_(is.name_),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream(ITstream&& is)
:
    Istream(static_cast<IOstreamOption>(is)),
    tokenList(std::move(static_cast<tokenList&>(is))),
    name_(std::move(is.name_)),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream
(
    IOstreamOption streamOpt,
    const string& name
)
:
    Istream(streamOpt.format(), streamOpt.version()),
    tokenList(),
    name_(name),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream
(
    const Foam::zero,
    const string& name,
    IOstreamOption streamOpt
)
:
    ITstream(streamOpt, name)
{}


Foam::ITstream::ITstream
(
    const string& name,
    const UList<token>& tokens,
    IOstreamOption streamOpt
)
:
    Istream(streamOpt.format(), streamOpt.version()),
    tokenList(tokens),
    name_(name),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream
(
    const string& name,
    List<token>&& tokens,
    IOstreamOption streamOpt
)
:
    Istream(streamOpt.format(), streamOpt.version()),
    tokenList(std::move(tokens)),
    name_(name),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream
(
    const UList<char>& input,
    IOstreamOption streamOpt,
    const string& name
)
:
    ITstream(streamOpt, name)
{
    UIListStream is(input, streamOpt);

    parseStream(is, static_cast<tokenList&>(*this));
    ITstream::rewind();
}


Foam::ITstream::ITstream
(
    const std::string& input,
    IOstreamOption streamOpt,
    const string& name
)
:
    ITstream(streamOpt, name)
{
    UIListStream is(input.data(), input.length(), streamOpt);

    parseStream(is, static_cast<tokenList&>(*this));
    ITstream::rewind();
}


Foam::ITstream::ITstream
(
    const char* input,
    IOstreamOption streamOpt,
    const string& name
)
:
    ITstream(streamOpt, name)
{
    UIListStream is(input, strlen(input), streamOpt);

    parseStream(is, static_cast<tokenList&>(*this));
    ITstream::rewind();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ITstream::print(Ostream& os) const
{
    os  << "ITstream : " << name_.c_str() << ", line ";

    const tokenList& toks = *this;

    if (toks.empty())
    {
        os  << lineNumber();
    }
    else
    {
        os  << toks.first().lineNumber();

        if (toks.first().lineNumber() < toks.last().lineNumber())
        {
            os  << '-' << toks.last().lineNumber();
        }
    }
    os  << ", ";

    IOstream::print(os);
}


std::string Foam::ITstream::toString() const
{
    // NOTE: may wish to have special handling if there is a single token
    // and it is already a string or word

    OStringStream buf;
    unsigned i = 0;
    for (const token& tok : *this)
    {
        if (i++)
        {
            buf << ' ';
        }
        buf << tok;
    }

    return buf.str();
}


const Foam::token& Foam::ITstream::peekFirst() const
{
    return peekTokenAt(*this, 0);
}


const Foam::token& Foam::ITstream::peekLast() const
{
    return peekTokenAt(*this, tokenList::size()-1);
}


const Foam::token& Foam::ITstream::peek() const
{
    // Use putback token if it exists
    if (Istream::hasPutback())
    {
        return Istream::peekBack();
    }

    return peekTokenAt(*this, tokenIndex_);
}


void Foam::ITstream::seek(label pos)
{
    lineNumber_ = 0;
    const tokenList& toks = *this;
    const label nToks = toks.size();

    if (!pos)
    {
        // Seek begin (rewind)
        tokenIndex_ = 0;

        if (nToks)
        {
            lineNumber_ = toks.first().lineNumber();
        }

        setOpened();
        setGood();
    }
    else if (pos < 0 || pos >= nToks)
    {
        // Seek end or seek is out of range
        tokenIndex_ = nToks;

        if (nToks)
        {
            lineNumber_ = toks.last().lineNumber();
        }

        setEof();
    }
    else
    {
        // Seek middle (from the beginning)
        tokenIndex_ = pos;

        if (nToks)
        {
            lineNumber_ = toks[tokenIndex_].lineNumber();
        }

        setOpened();
        setGood();
    }
}


void Foam::ITstream::skip(label n)
{
    const tokenList& toks = *this;
    const label nToks = toks.size();

    if (n < 0)
    {
        // Move backwards
        while (n++ && tokenIndex_)
        {
            --tokenIndex_;
        }

        if (tokenIndex_ < nToks)
        {
            lineNumber_ = toks[tokenIndex_].lineNumber();
            setOpened();
            setGood();
        }
    }
    else if (n > 0)
    {
        // Move forward
        while (n-- && tokenIndex_ < nToks)
        {
            ++tokenIndex_;
        }

        if (tokenIndex_ < nToks)
        {
            lineNumber_ = toks[tokenIndex_].lineNumber();
            setOpened();
            setGood();
        }
        else
        {
            setEof();
        }
    }
}


Foam::Istream& Foam::ITstream::read(token& tok)
{
    // Use putback token if it exists
    if (Istream::getBack(tok))
    {
        lineNumber_ = tok.lineNumber();
        return *this;
    }

    tokenList& toks = *this;
    const label nToks = toks.size();

    if (tokenIndex_ < nToks)
    {
        tok = toks[tokenIndex_++];
        lineNumber_ = tok.lineNumber();

        if (tokenIndex_ == nToks)
        {
            setEof();
        }
    }
    else
    {
        if (eof())
        {
            FatalIOErrorInFunction(*this)
                << "attempt to read beyond EOF"
                << exit(FatalIOError);
            setBad();
        }
        else
        {
            setEof();
        }

        tok.reset();

        if (nToks)
        {
            tok.lineNumber(toks.last().lineNumber());
        }
        else
        {
            tok.lineNumber(this->lineNumber());
        }
    }

    return *this;
}


Foam::Istream& Foam::ITstream::read(char&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(word&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(string&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(label&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(floatScalar&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(doubleScalar&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::readRaw(char*, std::streamsize)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(char*, std::streamsize)
{
    NotImplemented;
    return *this;
}


void Foam::ITstream::rewind()
{
    seek(0);
}


void Foam::ITstream::append(const token& t, const bool lazy)
{
    reserveCapacity(tokenIndex_ + 1, lazy);
    tokenList& toks = *this;

    toks[tokenIndex_] = t;  // copy append
    ++tokenIndex_;
}


void Foam::ITstream::append(token&& t, const bool lazy)
{
    reserveCapacity(tokenIndex_ + 1, lazy);
    tokenList& toks = *this;

    toks[tokenIndex_] = std::move(t);  // move append
    ++tokenIndex_;
}


void Foam::ITstream::append(const UList<token>& newTokens, const bool lazy)
{
    reserveCapacity(tokenIndex_ + newTokens.size(), lazy);
    tokenList& toks = *this;

    for (const token& t : newTokens)
    {
        toks[tokenIndex_] = t;  // copy append
        ++tokenIndex_;
    }
}


void Foam::ITstream::append(List<token>&& newTokens, const bool lazy)
{
    reserveCapacity(tokenIndex_ + newTokens.size(), lazy);
    tokenList& toks = *this;

    for (token& t : newTokens)
    {
        toks[tokenIndex_] = std::move(t);  // move append
        ++tokenIndex_;
    }

    newTokens.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::ITstream::operator=(const ITstream& is)
{
    // Self-assignment is a no-op
    if (this != &is)
    {
        Istream::operator=(is);
        tokenList::operator=(is);
        name_ = is.name_;
        rewind();
    }
}


void Foam::ITstream::operator=(const UList<token>& toks)
{
    tokenList::operator=(toks);
    rewind();
}


void Foam::ITstream::operator=(List<token>&& toks)
{
    tokenList::operator=(std::move(toks));
    rewind();
}


// ************************************************************************* //
