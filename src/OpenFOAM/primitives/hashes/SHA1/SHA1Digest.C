/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "SHA1Digest.H"
#include "IOstreams.H"
#include <cstring>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::SHA1Digest Foam::SHA1Digest::null;

static const char hexChars[] = "0123456789abcdef";

// The char '0' == 0
static constexpr int offsetZero = int('0');

// The char 'A' (or 'a') == 10
static constexpr int offsetUpper = int('A') - 10;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read hexadecimal value, ignoring leading or intermediate '_'
static unsigned char readHexDigit(Istream& is)
{
    // Silently ignore leading or intermediate '_'
    char c = 0;
    do
    {
        is.read(c);
    }
    while (c == '_');

    if (isdigit(c))
    {
        return int(c) - offsetZero;
    }
    else if (!isxdigit(c))
    {
        FatalIOErrorInFunction(is)
            << "Illegal hex digit: '" << c << "'"
            << exit(FatalIOError);
    }

    return toupper(c) - offsetUpper;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SHA1Digest::SHA1Digest()
{
    clear();
}


Foam::SHA1Digest::SHA1Digest(Istream& is)
{
    clear();
    read(is);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SHA1Digest::clear()
{
    dig_.fill(0);  // Same as memset(dig_.data(), 0, dig_.size());
}


bool Foam::SHA1Digest::empty() const
{
    for (const auto& byteVal : dig_)
    {
        if (byteVal)
        {
            return false;
        }
    }

    return true;
}


std::string Foam::SHA1Digest::str(const bool prefixed) const
{
    std::string buf;
    unsigned nChar = 0;

    if (prefixed)
    {
        buf.resize(1 + 2*dig_.size());
        buf[nChar++] = '_';
    }
    else
    {
        buf.resize(2*dig_.size());
    }

    for (const auto& byteVal : dig_)
    {
        buf[nChar++] = hexChars[((byteVal >> 4) & 0xF)];  // Upper nibble
        buf[nChar++] = hexChars[(byteVal & 0xF)];         // Lower nibble
    }

    return buf;
}


Foam::Istream& Foam::SHA1Digest::read(Istream& is)
{
    for (auto& byteVal : dig_)
    {
        const unsigned char upp = readHexDigit(is);
        const unsigned char low = readHexDigit(is);

        byteVal = (upp << 4) + low;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::SHA1Digest::write(Ostream& os, const bool prefixed) const
{
    if (prefixed)
    {
        os.write('_');
    }

    for (const auto& byteVal : dig_)
    {
        os.write(hexChars[((byteVal >> 4) & 0xF)]);  // Upper nibble
        os.write(hexChars[(byteVal & 0xF)]);         // Lower nibble
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::SHA1Digest::operator==(const SHA1Digest& rhs) const
{
    return (dig_ == rhs.dig_);
}


bool Foam::SHA1Digest::operator==(const std::string& hexdigits) const
{
    // Null or empty string is not an error - interpret as '0000..'
    if (hexdigits.empty())
    {
        return empty();
    }

    // Skip possible '_' prefix
    unsigned nChar = 0;
    if (hexdigits[0] == '_')
    {
        ++nChar;
    }

    // Incorrect length - can never match
    if (hexdigits.size() != nChar + 2*dig_.size())
    {
        return false;
    }

    for (const auto& byteVal : dig_)
    {
        const char upp = hexChars[((byteVal >> 4) & 0xF)];  // Upper nibble
        const char low = hexChars[(byteVal & 0xF)];         // Lower nibble

        if (upp != hexdigits[nChar++]) return false;
        if (low != hexdigits[nChar++]) return false;
    }

    return true;
}


bool Foam::SHA1Digest::operator==(const char* hexdigits) const
{
    // Null or empty string is not an error - interpret as '0000..'
    if (!hexdigits || !*hexdigits)
    {
        return empty();
    }

    // Skip possible '_' prefix
    unsigned nChar = 0;
    if (hexdigits[0] == '_')
    {
        ++nChar;
    }

    // Incorrect length - can never match
    if (strlen(hexdigits) != nChar + 2*dig_.size())
    {
        return false;
    }

    for (const auto& byteVal : dig_)
    {
        const char upp = hexChars[((byteVal >> 4) & 0xF)];
        const char low = hexChars[(byteVal & 0xF)];

        if (upp != hexdigits[nChar++]) return false;
        if (low != hexdigits[nChar++]) return false;
    }

    return true;
}


bool Foam::SHA1Digest::operator!=(const SHA1Digest& rhs) const
{
    return !operator==(rhs);
}


bool Foam::SHA1Digest::operator!=(const std::string& rhs) const
{
    return !operator==(rhs);
}


bool Foam::SHA1Digest::operator!=(const char* rhs) const
{
    return !operator==(rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, SHA1Digest& dig)
{
    return dig.read(is);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const SHA1Digest& dig)
{
    return dig.write(os);
}


// ************************************************************************* //
