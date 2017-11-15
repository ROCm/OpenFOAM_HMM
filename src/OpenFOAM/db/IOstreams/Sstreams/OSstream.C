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

#include "error.H"
#include "token.H"
#include "OSstream.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OSstream::write(const token& tok)
{
    // Direct token handling only for some types

    switch (tok.type())
    {
        case token::tokenType::FLAG :
        {
            // silently consume the flag
            return true;
        }

        case token::tokenType::VERBATIMSTRING :
        {
            write(char(token::HASH));
            write(char(token::BEGIN_BLOCK));
            writeQuoted(tok.stringToken(), false);
            write(char(token::HASH));
            write(char(token::END_BLOCK));

            return true;
        }

        case token::tokenType::VARIABLE :
        {
            writeQuoted(tok.stringToken(), false);

            return true;
        }

        default:
            break;
    }

    return false;
}


Foam::Ostream& Foam::OSstream::write(const char c)
{
    os_ << c;
    if (c == token::NL)
    {
        ++lineNumber_;
    }
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const char* str)
{
    lineNumber_ += stringOps::count(str, token::NL);
    os_ << str;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const word& str)
{
    os_ << str;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::writeQuoted
(
    const std::string& str,
    const bool quoted
)
{
    if (!quoted)
    {
        // Output unquoted, only advance line number on newline
        lineNumber_ += stringOps::count(str, token::NL);
        os_ << str;

        setState(os_.rdstate());
        return *this;
    }


    // Output with surrounding quotes and backslash escaping
    os_ << token::BEGIN_STRING;

    unsigned backslash = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        const char c = *iter;

        if (c == '\\')
        {
            ++backslash;
            continue; // only output after escaped character is known
        }
        else if (c == token::NL)
        {
            ++lineNumber_;
            ++backslash;    // backslash escape for newline
        }
        else if (c == token::END_STRING)
        {
            ++backslash;    // backslash escape for quote
        }

        // output all pending backslashes
        while (backslash)
        {
            os_ << '\\';
            --backslash;
        }

        os_ << c;
    }

    // silently drop any trailing backslashes
    // they would otherwise appear like an escaped end-quote
    os_ << token::END_STRING;

    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const string& str)
{
    return writeQuoted(str, true);
}


Foam::Ostream& Foam::OSstream::write(const int32_t val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const int64_t val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const floatScalar val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const doubleScalar val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write
(
    const char* data,
    const std::streamsize count
)
{
    beginRaw(count);
    writeRaw(data, count);
    endRaw();

    return *this;
}


Foam::Ostream& Foam::OSstream::beginRaw
(
    const std::streamsize count
)
{
    if (format() != BINARY)
    {
        FatalIOErrorInFunction(*this)
            << "stream format not binary"
            << abort(FatalIOError);
    }

    os_ << token::BEGIN_LIST;

    setState(os_.rdstate());

    return *this;
}


Foam::Ostream& Foam::OSstream::writeRaw
(
    const char* data,
    std::streamsize count
)
{
    // No check for format() == BINARY since this is either done in the
    // beginRaw() method, or the caller knows what they are doing.

    os_.write(data, count);
    setState(os_.rdstate());

    return *this;
}


Foam::Ostream& Foam::OSstream::endRaw()
{
    os_ << token::END_LIST;
    setState(os_.rdstate());

    return *this;
}


void Foam::OSstream::indent()
{
    for (unsigned short i = 0; i < indentLevel_*indentSize_; ++i)
    {
        os_ << ' ';
    }
}


void Foam::OSstream::flush()
{
    os_.flush();
}


void Foam::OSstream::endl()
{
    write('\n');
    os_.flush();
}


std::ios_base::fmtflags Foam::OSstream::flags() const
{
    return os_.flags();
}


std::ios_base::fmtflags Foam::OSstream::flags(const ios_base::fmtflags f)
{
    return os_.flags(f);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int Foam::OSstream::width() const
{
    return os_.width();
}


int Foam::OSstream::width(const int w)
{
    return os_.width(w);
}


int Foam::OSstream::precision() const
{
    return os_.precision();
}


int Foam::OSstream::precision(const int p)
{
    return os_.precision(p);
}


// ************************************************************************* //
