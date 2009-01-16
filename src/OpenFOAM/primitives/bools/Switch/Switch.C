/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Switch.H"
#include "error.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// NB: values chosen such that bitwise '&' 0x1 yields the bool value

const char* Foam::Switch::names[Foam::Switch::INVALID+1] =
{
    "false", "true",
    "off",   "on",
    "no",    "yes",
    "n",     "y",
    "invalid"
};


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //


Foam::Switch::switchType Foam::Switch::asEnum
(
    const word& val,
    const bool ignoreError
)
{
    for (int sw = 0; sw < INVALID; sw++)
    {
        if (val == names[sw])
        {
            if (sw == NO_1)
            {
                return NO;
            }
            else if (sw == YES_1)
            {
                return YES;
            }
            else
            {
                return switchType(sw);
            }
        }
    }

    if (!ignoreError)
    {
        FatalErrorIn("Switch::asEnum(const word&)")
            << "unknown switch word " << val
            << abort(FatalError);
    }

    return INVALID;
}


Foam::Switch::switchType Foam::Switch::asEnum(const bool val)
{
    return val ? ON : OFF;
}


bool Foam::Switch::asBool
(
    const word& val,
    const bool ignoreError
)
{
    switchType sw = asEnum(val, true);

    // catch error here
    if (sw == INVALID && !ignoreError)
    {
        FatalErrorIn("Switch::asBool(const word&)")
            << "unknown switch word " << val
            << abort(FatalError);
    }

    return asBool(sw);
}


bool Foam::Switch::asBool(const switchType& val)
{
    return (val & 0x1);
}


Foam::word Foam::Switch::asWord(const bool val)
{
    return word((val ? names[ON] : names[OFF]), false);
}


Foam::word Foam::Switch::asWord(const switchType& val)
{
    return word(names[val], false);
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

bool Foam::Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent(name, bool_);
}


// ************************************************************************* //
