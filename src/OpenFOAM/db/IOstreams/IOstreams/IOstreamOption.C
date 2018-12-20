/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "IOstreamOption.H"
#include "error.H"
#include "Enum.H"
#include "Switch.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::IOstreamOption::versionNumber
    Foam::IOstreamOption::originalVersion(0,5);

const Foam::IOstreamOption::versionNumber
    Foam::IOstreamOption::currentVersion(2,0);


const Foam::Enum
<
    Foam::IOstreamOption::streamFormat
>
Foam::IOstreamOption::formatNames
({
    { streamFormat::ASCII, "ascii" },
    { streamFormat::BINARY, "binary" },
});


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

Foam::IOstreamOption::streamFormat
Foam::IOstreamOption::formatEnum(const word& formatName)
{
    // Handle bad input graciously
    if (formatNames.found(formatName))
    {
        return formatNames[formatName];
    }

    WarningInFunction
        << "Unknown format specifier '" << formatName
        << "', using 'ascii'" << endl;

    return streamFormat::ASCII;
}


Foam::IOstreamOption::compressionType
Foam::IOstreamOption::compressionEnum(const word& compName)
{
    // Handle bad input graciously

    const Switch sw(compName, true);
    if (sw.valid())
    {
        return
        (
            sw
          ? compressionType::COMPRESSED
          : compressionType::UNCOMPRESSED
        );
    }

    WarningInFunction
        << "Unknown compression specifier '" << compName
        << "', assuming no compression" << endl;

    return compressionType::UNCOMPRESSED;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const IOstreamOption::streamFormat& sf
)
{
    os << IOstreamOption::formatNames[sf];
    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const IOstreamOption::versionNumber& vn
)
{
    // Emit as char sequence instead of as individual characters
    // in case this is needed for sending in parallel.
    os  << vn.str().c_str();
    return os;
}


// ************************************************************************* //
