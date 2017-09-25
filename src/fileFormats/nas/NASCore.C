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

#include "NASCore.H"
#include "parsing.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::scalar Foam::fileFormats::NASCore::readNasScalar(const string& str)
{
    const auto signPos = str.find_last_of("+-");

    if
    (
        signPos == std::string::npos
     || signPos == 0
     || str[signPos-1] == 'E' || str[signPos-1] == 'e'
     || isspace(str[signPos-1])
    )
    {
        // A normal number format
        return readScalar(str);
    }


    // Nastran compact number format.
    // Eg, "1234-2" instead of "1234E-2"

    scalar value = 0;
    int exponent = 0; // Any integer

    if
    (
        readScalar(str.substr(0, signPos), value)   // Mantissa
     && readInt(str.substr(signPos),  exponent)     // Exponent (with sign)
    )
    {
        // Note: this does not catch underflow/overflow
        // (especially when scalar is a float)
        value *= ::pow(10, exponent);
    }
    else
    {
        FatalIOErrorInFunction("unknown")
            << parsing::errorNames[parsing::errorType::GENERAL] << str
            << exit(FatalIOError);

        value = 0;
    }

    return value;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::NASCore::NASCore()
{}


// ************************************************************************* //
