/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
#include "dictionary.H"
#include "Enum.H"
#include "Switch.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::IOstreamOption::versionNumber Foam::IOstreamOption::currentVersion;

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
Foam::IOstreamOption::formatEnum
(
    const word& formatName,
    const streamFormat deflt
)
{
    // Handle bad input graciously. A no-op for an empty string

    if (!formatName.empty())
    {
        if (formatNames.found(formatName))
        {
            return formatNames[formatName];
        }

        // Fall-through to warning

        WarningInFunction
            << "Unknown format specifier '" << formatName
            << "', using '" << formatNames[deflt] << "'\n";
    }

    return deflt;
}


Foam::IOstreamOption::streamFormat
Foam::IOstreamOption::formatEnum
(
    const word& key,
    const dictionary& dict,
    const streamFormat deflt
)
{
    return formatNames.getOrDefault(key, dict, deflt, true); // failsafe=true
}


Foam::IOstreamOption::compressionType
Foam::IOstreamOption::compressionEnum
(
    const word& compName,
    const compressionType deflt
)
{
    // Handle bad input graciously. A no-op for an empty string

    if (!compName.empty())
    {
        const Switch sw = Switch::find(compName);

        if (sw.good())
        {
            return
            (
                sw
              ? compressionType::COMPRESSED
              : compressionType::UNCOMPRESSED
            );
        }

        // Fall-through to warning

        WarningInFunction
            << "Unknown compression specifier '" << compName
            << "', using compression "
            << (deflt ? "on" : "off" ) << nl;
    }

    return deflt;
}


Foam::IOstreamOption::compressionType
Foam::IOstreamOption::compressionEnum
(
    const word& key,
    const dictionary& dict,
    const compressionType deflt
)
{
    return
    (
        Switch(key, dict, Switch(bool(deflt)), true) // failsafe=true
      ? compressionType::COMPRESSED
      : compressionType::UNCOMPRESSED
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOstreamOption::versionNumber::versionNumber(const std::string& verNum)
:
    versionNumber(readFloat(verNum))
{}


Foam::IOstreamOption::versionNumber::versionNumber(const token& tok)
:
    versionNumber()
{
    if (tok.isStringType())
    {
        (*this) = versionNumber(tok.stringToken());
    }
    else if (tok.isNumber())
    {
        // Accept integer or floating-point
        // Eg, '2.0' becomes '2' after foamDictionary -expand
        (*this) = versionNumber(float(tok.number()));
    }
    else
    {
        WarningInFunction
            << "Wrong token for version - expected word/number, found "
            << tok.info() << nl;
    }
}


Foam::IOstreamOption::versionNumber::versionNumber
(
    const word& key,
    const dictionary& dict
)
:
    versionNumber()
{
    token tok;

    if (dict.readIfPresent<token>(key, tok, keyType::LITERAL))
    {
        (*this) = versionNumber(tok);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const IOstreamOption::streamFormat& fmt
)
{
    os << IOstreamOption::formatNames[fmt];
    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const IOstreamOption::versionNumber& ver
)
{
    // Emit unquoted char sequence (eg, word)
    // for correct behaviour when sending in parallel

    os.writeQuoted(ver.str(), false);
    return os;
}


// ************************************************************************* //
