/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "Switch.H"
#include "scalar.H"
#include "error.H"
#include "dictionary.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace
{
static_assert
(
    Foam::Switch::INVALID+1 == 9,
    "Switch::switchType does not have 9 entries"
);

//- The names corresponding to the Switch::switchType enumeration.
//  Includes extra entries for "invalid".
static const char* names[9] =
{
    "false", "true",
    "no",    "yes",
    "off",   "on",
    "none",  "any",
    "invalid"
};

} // End anonymous namespace

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

const char* Foam::Switch::name(const bool b) noexcept
{
    return names[(b ? 1 : 0)];
}


Foam::Switch::switchType Foam::Switch::parse
(
    const std::string& str,
    bool allowBad
)
{
    switch (str.size())
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
        case 3: // (off|yes|any)
        {
            if (str == names[switchType::OFF]) return switchType::OFF;
            if (str == names[switchType::YES]) return switchType::YES;
            if (str == names[switchType::ANY]) return switchType::ANY;
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

    if (!allowBad)
    {
        FatalErrorInFunction
            << "Unknown switch word " << str << nl
            << abort(FatalError);
    }

    return switchType::INVALID;
}


Foam::Switch Foam::Switch::getOrAddToDict
(
    const word& name,
    dictionary& dict,
    const Switch deflt
)
{
    return dict.getOrAdd<Switch>(name, deflt);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Switch::Switch(const float val, const float tol)
:
    switch_((mag(val) > tol) ? switchType::TRUE : switchType::FALSE)
{}


Foam::Switch::Switch(const double val, const double tol)
:
    switch_((mag(val) > tol) ? switchType::TRUE : switchType::FALSE)
{}


Foam::Switch::Switch
(
    const word& key,
    const dictionary& dict
)
{
    const word str(dict.get<word>(key, keyType::LITERAL));

    (*this) = parse(str, true);

    if (!valid())
    {
        FatalIOErrorInFunction(dict)
            << "Expected 'true/false', 'on/off' ... found " << str << nl
            << exit(FatalIOError);
    }
}


Foam::Switch::Switch
(
    const word& key,
    const dictionary& dict,
    const Switch deflt
)
:
    Switch(deflt)
{
    const entry* eptr = dict.findEntry(key, keyType::LITERAL);

    if (eptr)
    {
        const word str(eptr->get<word>());

        (*this) = parse(str, true);

        if (!valid())
        {
            // Found entry, but was bad input

            FatalIOErrorInFunction(dict)
                << "Expected 'true/false', 'on/off' ... found " << str << nl
                << exit(FatalIOError);
        }
    }
}


Foam::Switch::Switch(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::Switch::valid() const noexcept
{
    return switch_ != switchType::INVALID;
}


Foam::Switch::switchType Foam::Switch::type() const noexcept
{
    return switchType(switch_);
}


const char* Foam::Switch::c_str() const noexcept
{
    return names[(switch_ & 0x0F)];
}


std::string Foam::Switch::str() const
{
    return names[(switch_ & 0x0F)];
}


bool Foam::Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent<Switch>(name, *this);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, Switch& sw)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get bool"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        sw = bool(t.labelToken());
    }
    else if (t.isWord())
    {
        // Permit invalid value, but catch immediately for better messages
        sw = Switch(t.wordToken(), true);

        if (!sw.valid())
        {
            FatalIOErrorInFunction(is)
                << "Expected 'true/false', 'on/off' ... found "
                << t.wordToken()
                << exit(FatalIOError);
            is.setBad();
            return is;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected bool, found "
            << t.info()
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const Switch& sw)
{
    os << sw.c_str();
    return os;
}


// ************************************************************************* //
