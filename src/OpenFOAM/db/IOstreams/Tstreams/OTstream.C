/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
#include "OTstream.H"
#include <cctype>

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::OTstream::write(const token& tok)
{
    if (tok.good())
    {
        append(tok);
        return true;
    }

    return false;
}


Foam::Ostream& Foam::OTstream::write(const char c)
{
    if (!std::isspace(c) && std::isprint(c))
    {
        // Should generally work, but need to verify corner cases
        append(token(token::punctuationToken(c)));
    }

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const char* str)
{
    const word nonWhiteChars(string::validate<word>(str));

    if (nonWhiteChars.size() == 1)
    {
        // Like punctuation
        write(nonWhiteChars[0]);
    }
    else if (nonWhiteChars.size())
    {
        // As a word
        write(nonWhiteChars);
    }

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const word& str)
{
    append(token(str)); // tokenType::WORD

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const string& str)
{
    append(token(str)); // tokenType::STRING

    return *this;
}


Foam::Ostream& Foam::OTstream::writeQuoted
(
    const std::string& str,
    const bool quoted
)
{
    if (quoted)
    {
        append(token(string(str))); // tokenType::STRING
    }
    else if (!str.empty())
    {
        append(token(word(str, false))); // tokenType::WORD
    }

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const int32_t val)
{
    append(token(label(val))); // tokenType::LABEL

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const int64_t val)
{
    append(token(label(val))); // tokenType::LABEL

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const floatScalar val)
{
    append(token(val)); // tokenType::FLOAT

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const doubleScalar val)
{
    append(token(val)); // tokenType::DOUBLE

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const char* data, std::streamsize count)
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    NotImplemented;

    return *this;
}


Foam::Ostream& Foam::OTstream::writeRaw
(
    const char* data,
    std::streamsize count
)
{
    // No check for format() == BINARY since this is either done in the
    // beginRawWrite() method, or the caller knows what they are doing.

    NotImplemented;

    return *this;
}


bool Foam::OTstream::beginRawWrite(std::streamsize count)
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    NotImplemented;

    return true;
}


void Foam::OTstream::print(Ostream& os) const
{
    os  << "OTstream : " << name().c_str() << ", " << size() << " tokens, ";
    IOstream::print(os);
}


// ************************************************************************* //
