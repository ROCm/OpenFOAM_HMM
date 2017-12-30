/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "word.H"
#include "debug.H"
#include <cctype>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::word::typeName = "word";
int Foam::word::debug(Foam::debug::debugSwitch(word::typeName, 0));
const Foam::word Foam::word::null;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::word::validate(const std::string& s, const bool prefix)
{
    word out;
    out.resize(s.size() + (prefix ? 1 : 0));

    std::string::size_type count = 0;

    // As per validate, but optionally detect if the first character
    // is a digit, which we'd like to avoid having since this will
    // cause parse issues when read back later.
    for (auto iter = s.cbegin(); iter != s.cend(); ++iter)
    {
        const char c = *iter;

        if (word::valid(c))
        {
            if (!count && prefix && isdigit(c))
            {
                // First valid character was a digit - prefix with '_'
                out[count++] = '_';
            }

            out[count++] = c;
        }
    }

    out.resize(count);

    return out;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::word::lessExt() const
{
    const size_type i = find_ext();

    if (i == npos)
    {
        return *this;
    }

    return substr(0, i);
}


Foam::word Foam::word::ext() const
{
    return string::ext();
}


Foam::word& Foam::word::ext(const word& ending)
{
    string::ext(ending);
    return *this;
}


bool Foam::word::hasExt(const word& ending) const
{
    return string::hasExt(ending);
}


bool Foam::word::hasExt(const wordRe& ending) const
{
    return string::hasExt(ending);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::word Foam::operator&(const word& a, const word& b)
{
    if (a.size())
    {
        if (b.size())
        {
            // Two non-empty words: can concatenate and perform camel case
            word camelCase(a + b);
            camelCase[a.size()] = char(toupper(b[0]));

            return camelCase;
        }

        return a;
    }

    // Or, if the first string is empty (or both are empty)

    return b;
}


// ************************************************************************* //
