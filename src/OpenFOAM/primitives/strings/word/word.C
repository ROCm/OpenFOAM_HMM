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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::word::typeName = "word";
int Foam::word::debug(Foam::debug::debugSwitch(word::typeName, 0));
const Foam::word Foam::word::null;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::word::validated(const std::string& s)
{
    std::string::size_type count = 0;
    bool prefix = false;

    // Count number of valid characters and detect if the first character
    // happens to be a digit, which we'd like to avoid having since this
    // will cause parse issues when read back later.
    for (std::string::const_iterator it = s.cbegin(); it != s.cend(); ++it)
    {
        const char c = *it;

        if (word::valid(c))
        {
            if (!count && isdigit(c))
            {
                // First valid character was a digit - prefix with '_'
                prefix = true;
                ++count;
            }

            ++count;
        }
    }

    if (count == s.size() && !prefix)
    {
        return word(s, false);  // Already checked, can just return as word
    }

    word out;
    out.resize(count);
    count = 0;

    // Copy valid content.
    if (prefix)
    {
        out[count++] = '_';
    }

    for (std::string::const_iterator it = s.cbegin(); it != s.cend(); ++it)
    {
        const char c = *it;

        if (word::valid(c))
        {
            out[count++] = c;
        }
    }

    out.resize(count);

    return out;
}


// ************************************************************************* //
