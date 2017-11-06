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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ITstream::toTokenList(ISstream& is)
{
    tokenIndex_ = 0;

    token tok;

    while (!is.read(tok).bad() && tok.good())
    {
        newElmt(tokenIndex()++) = std::move(tok);
    }

    tokenList::setSize(tokenIndex());

    setOpened();
    ITstream::rewind();
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
    tokenList(16, token::undefinedToken),
    name_(name),
    tokenIndex_(0)
{
    UIListStream is(input, format, version);

    toTokenList(is);
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
    tokenList(16, token::undefinedToken),
    name_(name),
    tokenIndex_(0)
{
    UIListStream is(input.data(), input.size(), format, version);

    toTokenList(is);
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
    tokenList(16, token::undefinedToken),
    name_(name),
    tokenIndex_(0)
{
    const size_t len = strlen(input);
    UIListStream is(input, len, format, version);

    toTokenList(is);
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


Foam::Istream& Foam::ITstream::read(token& t)
{
    // Return the put back token if it exists
    if (Istream::getBack(t))
    {
        lineNumber_ = t.lineNumber();
        return *this;
    }

    if (tokenIndex_ < size())
    {
        t = operator[](tokenIndex_++);
        lineNumber_ = t.lineNumber();

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

        t = token::undefinedToken;

        if (size())
        {
            t.lineNumber() = tokenList::last().lineNumber();
        }
        else
        {
            t.lineNumber() = lineNumber();
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

    if (size())
    {
        lineNumber_ = tokenList::first().lineNumber();
    }

    setGood();
}


// ************************************************************************* //
