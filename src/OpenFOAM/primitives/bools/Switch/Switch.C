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

Foam::Switch Foam::Switch::lookupOrAddToDict
(
    const word& name,
    dictionary& dict,
    const Switch& defaultValue
)
{
    return dict.lookupOrAddDefault<Switch>(name, defaultValue);
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// NOTE: possible alternative implementation
// make a direct bool, handle assignments and use switchTypes instead of word
// for the word representation ...
//
//   //- Possible word representions
//   enum switchTypes
//   {
//       OFF = 0, ON = 1,
//       FALSE = 2, TRUE = 3,
//       NO = 4, YES = 5
//   };


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

bool Foam::Switch::asBool(const word& val)
{
    if (val == "on" || val == "true" || val == "yes" || val == "y")
    {
        return true;
    }
    else if (val == "off" || val == "false" || val == "no" || val == "n")
    {
        return false;
    }
    else
    {
        FatalErrorIn("Switch::asBool(const word&)")
            << "unknown switch word " << val
            << abort(FatalError);
    }

    return false;
}


Foam::word Foam::Switch::asWord(const bool val)
{
    return word((val ? "on" : "off"), false);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent(name, bool_);
}


// ************************************************************************* //
