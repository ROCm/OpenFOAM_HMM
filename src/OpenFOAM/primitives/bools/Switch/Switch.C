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

#include "Switch.H"
#include "error.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::Switch::names[Foam::Switch::INVALID+1] =
{
    "false", "true",
    "no",    "yes",
    "off",   "on",
    "none",  "true",  // Is there a reasonable counterpart to "none"?
    "invalid"
};


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

Foam::Switch::switchType Foam::Switch::asEnum
(
    const std::string& str,
    const bool allowInvalid
)
{
    const std::string::size_type len = str.size();
    switch (len)
    {
        case 1: // (f|n|t|y) - single-character forms
        {
            switch (str[0])
            {
                case 'f': return switchType::FALSE;
                case 'n': return switchType::NO;
                case 't': return switchType::TRUE;
                case 'y': return switchType::YES;
            }
            break;
        }
        case 2: // (no|on)
        {
            if (str == names[switchType::NO]) return switchType::NO;
            if (str == names[switchType::ON]) return switchType::ON;
            break;
        }
        case 3: // (off|yes)
        {
            if (str == names[switchType::OFF]) return switchType::OFF;
            if (str == names[switchType::YES]) return switchType::YES;
            break;
        }
        case 4: // (none|true)
        {
            if (str == names[switchType::NONE]) return switchType::NONE;
            if (str == names[switchType::TRUE]) return switchType::TRUE;
            break;
        }
        case 5: // (false)
        {
            if (str == names[switchType::FALSE]) return switchType::FALSE;
            break;
        }
    }

    if (!allowInvalid)
    {
        FatalErrorInFunction
            << "Unknown switch word " << str << nl
            << abort(FatalError);
    }

    return switchType::INVALID;
}


Foam::Switch Foam::Switch::lookupOrAddToDict
(
    const word& name,
    dictionary& dict,
    const Switch& defaultValue
)
{
    return dict.lookupOrAddDefault<Switch>(name, defaultValue);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::Switch::valid() const
{
    return switch_ <= switchType::NONE;
}


const char* Foam::Switch::asText() const
{
    return names[switch_];
}


bool Foam::Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent<Switch>(name, *this);
}


// ************************************************************************* //
