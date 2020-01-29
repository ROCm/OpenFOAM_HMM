/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::ITstream::parseStream(ISstream& is, tokenList& tokens)
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


Foam::tokenList Foam::ITstream::parse
(
    const UList<char>& input,
    streamFormat format
)
{
    UIListStream is(input, format, IOstream::currentVersion);

    tokenList tokens;
    parseStream(is, tokens);
    return tokens;
}


Foam::tokenList Foam::ITstream::parse
(
    const std::string& input,
    streamFormat format
)
{
    UIListStream is
    (
        input.data(),
        input.size(),
        format,
        IOstream::currentVersion
    );

    tokenList tokens;
    parseStream(is, tokens);
    return tokens;
}


Foam::tokenList Foam::ITstream::parse
(
    const char* input,
    streamFormat format
)
{
    UIListStream is(input, strlen(input), format, IOstream::currentVersion);

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

Foam::ITstream::ITstream
(
    const string& name,
    const UList<char>& input,
    streamFormat format,
    versionNumber version
)
:
    Istream(format, version),
    tokenList(),
    name_(name),
    tokenIndex_(0)
{
    UIListStream is(input, format, version);

    parseStream(is, static_cast<tokenList&>(*this));
    ITstream::rewind();
}


Foam::ITstream::ITstream
(
    const string& name,
    const std::string& input,
    streamFormat format,
    versionNumber version
)
:
    Istream(format, version),
    tokenList(),
    name_(name),
    tokenIndex_(0)
{
    UIListStream is(input.data(), input.size(), format, version);

    parseStream(is, static_cast<tokenList&>(*this));
    ITstream::rewind();
}


Foam::ITstream::ITstream
(
    const string& name,
    const char* input,
    streamFormat format,
    versionNumber version
)
:
    Istream(format, version),
    tokenList(),
    name_(name),
    tokenIndex_(0)
{
    UIListStream is(input, strlen(input), format, version);

    parseStream(is, static_cast<tokenList&>(*this));
    ITstream::rewind();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ITstream::print(Ostream& os) const
{
    os  << "ITstream : " << name_.c_str() << ", line ";

    if (size())
    {
        os  << tokenList::first().lineNumber();

        if (tokenList::first().lineNumber() < tokenList::last().lineNumber())
        {
            os  << '-' << tokenList::last().lineNumber();
        }
    }
    else
    {
        os  << lineNumber();
    }

    os  << ", ";

    IOstream::print(os);
}


std::string Foam::ITstream::toString() const
{
    OStringStream buf;

    const tokenList& tokens = *this;

    label len = tokens.size();

    // NOTE: may wish to have special handling if there is a single token
    // and it is already a string or word

    for (const token& tok : tokens)
    {
        buf << tok;

        if (--len)
        {
            buf << ' ';
        }
    }

    return buf.str();
}


Foam::Istream& Foam::ITstream::read(token& tok)
{
    // Return the put back token if it exists
    if (Istream::getBack(tok))
    {
        lineNumber_ = tok.lineNumber();
        return *this;
    }

    if (tokenIndex_ < size())
    {
        tok = operator[](tokenIndex_++);
        lineNumber_ = tok.lineNumber();

        if (tokenIndex_ == size())
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

        if (size())
        {
            tok.lineNumber() = tokenList::last().lineNumber();
        }
        else
        {
            tok.lineNumber() = lineNumber();
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


void Foam::ITstream::seek(label pos)
{
    lineNumber_ = 0;
    tokenList& toks = *this;

    if (!pos)
    {
        // Seek begin (rewind)
        tokenIndex_ = 0;

        if (!toks.empty())
        {
            lineNumber_ = toks.first().lineNumber();
        }

        setOpened();
        setGood();
    }
    else if (pos < 0 || pos >= toks.size())
    {
        // Seek end or seek is out of range
        tokenIndex_ = toks.size();

        if (!toks.empty())
        {
            lineNumber_ = toks.last().lineNumber();
        }

        setEof();
    }
    else
    {
        // Seek middle (from the beginning)
        tokenIndex_ = pos;

        if (!toks.empty())
        {
            lineNumber_ = toks[tokenIndex_].lineNumber();
        }

        setOpened();
        setGood();
    }
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


void Foam::ITstream::append(const tokenList& newTokens, const bool lazy)
{
    reserveCapacity(tokenIndex_ + newTokens.size(), lazy);
    tokenList& toks = *this;

    for (const token& t : newTokens)
    {
        toks[tokenIndex_] = t;  // copy append
        ++tokenIndex_;
    }
}


void Foam::ITstream::append(tokenList&& newTokens, const bool lazy)
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
    Istream::operator=(is);
    tokenList::operator=(is);
    name_ = is.name_;

    rewind();
}


void Foam::ITstream::operator=(const tokenList& toks)
{
    tokenList::operator=(toks);

    rewind();
}


void Foam::ITstream::operator=(tokenList&& toks)
{
    tokenList::operator=(std::move(toks));

    rewind();
}


// ************************************************************************* //
