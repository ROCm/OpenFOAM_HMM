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

#include "error.H"
#include "ITstream.H"
#include "UIListStream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::ITstream::parseStream(ISstream& is, tokenList& tokens)
{
    label nTok = 0;

    tokens.clear();
    tokens.setSize(64, token::undefinedToken);

    token tok;
    while (!is.read(tok).bad() && tok.good())
    {
        tokens.newElmt(nTok++) = std::move(tok);
    }

    tokens.setSize(nTok);

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
    os  << "ITstream : " << name_.c_str();

    if (size())
    {
        if (begin()->lineNumber() == rbegin()->lineNumber())
        {
            os  << ", line " << begin()->lineNumber() << ", ";
        }
        else
        {
            os  << ", lines " << begin()->lineNumber()
                << '-' << rbegin()->lineNumber() << ", ";
        }
    }
    else
    {
        os  << ", line " << lineNumber() << ", ";
    }

    IOstream::print(os);
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
            FatalIOErrorInFunction
            (
                *this
            )   << "attempt to read beyond EOF"
                << exit(FatalIOError);

            setBad();
        }
        else
        {
            setEof();
        }

        tok = token::undefinedToken;

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


Foam::Istream& Foam::ITstream::read(char*, std::streamsize)
{
    NotImplemented;
    return *this;
}


void Foam::ITstream::rewind()
{
    tokenIndex_ = 0;
    lineNumber_ = 0;

    if (size())
    {
        lineNumber_ = tokenList::first().lineNumber();
    }

    setOpened();
    setGood();
}


// ************************************************************************* //
