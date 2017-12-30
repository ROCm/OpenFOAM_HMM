/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
#include "IOmanip.H"
#include "Ostream.H"
#include "parsing.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fileFormats::NASCore::fieldFormat
>
Foam::fileFormats::NASCore::fieldFormatNames
{
    { fieldFormat::SHORT, "short" },
    { fieldFormat::LONG,  "long" },
    { fieldFormat::FREE,  "free" },
};


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::scalar Foam::fileFormats::NASCore::readNasScalar(const std::string& str)
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
     && readInt(str.substr(signPos), exponent)      // Exponent (with sign)
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


std::string Foam::fileFormats::NASCore::nextNasField
(
    const std::string& str,
    std::string::size_type& pos,
    std::string::size_type len
)
{
    const auto beg = pos;
    const auto end = str.find(',', pos);

    if (end == std::string::npos)
    {
        pos = beg + len;    // Continue after field width
    }
    else
    {
        len = (end - beg);  // Efffective width
        pos = end + 1;      // Continue after comma
    }

    return str.substr(beg, len);
}


void Foam::fileFormats::NASCore::setPrecision
(
    Ostream& os,
    const fieldFormat format
)
{
    os.setf(ios_base::scientific);

    // Capitalise the E marker
    os.setf(ios_base::uppercase);

    const label offset = 7;

    label prec = 16 - offset;
    switch (format)
    {
        case fieldFormat::SHORT :
        {
            prec = 8 - offset;
            break;
        }

        case fieldFormat::LONG :
        case fieldFormat::FREE :
        {
            prec = 16 - offset;
            break;
        }
    }

    os.precision(prec);
}


Foam::Ostream& Foam::fileFormats::NASCore::writeKeyword
(
    Ostream& os,
    const word& keyword,
    const fieldFormat format
)
{
    os.setf(ios_base::left);

    switch (format)
    {
        case fieldFormat::SHORT :
        {
            os  << setw(8) << keyword;
            break;
        }

        case fieldFormat::LONG :
        {
            os  << setw(8) << word(keyword + '*');
            break;
        }

        case fieldFormat::FREE :
        {
            os  << keyword;
            break;
        }
    }

    os.unsetf(ios_base::left);

    return os;
}


// ************************************************************************* //
