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

#include "NASCore.H"
#include "StringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::NASCore::NASCore()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::fileFormats::NASCore::parseNASCoord(const string& s)
{
    scalar value = 0;

    const size_t expSign = s.find_last_of("+-");

    if (expSign != std::string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar exponent = 0;

        // Parse as per strtod/strtof - allowing trailing space or [Ee]
        readScalar(s.substr(0, expSign).c_str(), value);  // mantissa
        readScalar(s.substr(expSign+1).c_str(),  exponent);

        if (s[expSign] == '-')
        {
            exponent = -exponent;
        }

        value *= ::pow(10, exponent);
    }
    else
    {
        readScalar(s.c_str(), value);
    }

    return value;
}


// ************************************************************************* //
