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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::Switch::wordValue(const bool val) const
{
    return word(val ? "on" : "off");
}


bool Foam::Switch::boolValue(const word& val) const
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
        FatalErrorIn("Switch::boolValue(const word&) const")
            << "unknown switch word " << val
            << abort(FatalError);
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent(name, logicalSwitch_);
}


// ************************************************************************* //
