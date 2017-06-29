/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "bool.H"
#include "Switch.H"
#include "error.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, bool& b)
{
    if (is.good())
    {
        b = Switch(is);
    }

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const bool b)
{
    // we could also write as text string without any difficulty
    // os  << (b ? "true" : "false");
    os.write(label(b));
    os.check(FUNCTION_NAME);
    return os;
}


bool Foam::readBool(Istream& is)
{
    bool b;
    is  >> b;

    return b;
}


// ************************************************************************* //
