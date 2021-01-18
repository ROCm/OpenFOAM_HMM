/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "IOobject.H"
#include "objectRegistry.H"
#include "foamVersion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// A banner corresponding to this:
//
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  VERSION                               |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

Foam::Ostream&
Foam::IOobject::writeBanner(Ostream& os, const bool noSyntaxHint)
{
    // The version padded with spaces to fit after "Version:  "
    // - initialized with zero-length string to detect if it has been populated
    static char paddedVersion[39] = "";

    if (!*paddedVersion)
    {
        // Populate: like strncpy but without trailing '\0'

        std::size_t len = foamVersion::version.length();
        if (len > 38)
        {
            len = 38;
        }

        std::memset(paddedVersion, ' ', 38);
        std::memcpy(paddedVersion, foamVersion::version.c_str(), len);
        paddedVersion[38] = '\0';
    }

    os  <<
        "/*--------------------------------";

    if (noSyntaxHint)
    {
        // Without syntax hint
        os  << "---------";
    }
    else
    {
        // With syntax hint
        os  << "*- C++ -*";
    }

    os  <<
        "----------------------------------*\\\n"
        "| =========                 |"
        "                                                 |\n"
        "| \\\\      /  F ield         |"
        " OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     |"
        " Version:  " << paddedVersion << "|\n"
        "|   \\\\  /    A nd           |"
        " Website:  www.openfoam.com                      |\n"
        "|    \\\\/     M anipulation  |"
        "                                                 |\n"
        "\\*-----------------------------------------"
        "----------------------------------*/\n";

    return os;
}


Foam::Ostream& Foam::IOobject::writeDivider(Ostream& os)
{
    os  <<
        "// * * * * * * * * * * * * * * * * * "
        "* * * * * * * * * * * * * * * * * * * * //\n";

    return os;
}


Foam::Ostream& Foam::IOobject::writeEndDivider(Ostream& os)
{
    os  << "\n\n"
        "// *****************************************"
        "******************************** //\n";

    return os;
}


bool Foam::IOobject::writeHeader
(
    Ostream& os,
    const word& objectType,
    const bool noArchAscii
) const
{
    if (!os.good())
    {
        InfoInFunction
            << "No stream open for write" << nl
            << os.info() << endl;

        return false;
    }

    IOobject::writeBanner(os)
        << "FoamFile" << nl
        << '{' << nl
        << "    version     " << os.version() << ';' << nl
        << "    format      " << os.format() << ';' << nl;

    if (os.format() == IOstream::BINARY || !noArchAscii)
    {
        // Arch information (BINARY: always, ASCII: can disable)
        os  << "    arch        " << foamVersion::buildArch << ';' << nl;
    }
    if (!note().empty())
    {
        os  << "    note        " << note() << ';' << nl;
    }

    os  << "    class       ";
    if (objectType.empty())
    {
        // Empty type not allowed - use 'dictionary' fallback
        os  << "dictionary";
    }
    else
    {
        os  << objectType;
    }
    os  << ';' << nl;

    os  << "    location    " << instance()/db().dbDir()/local() << ';' << nl
        << "    object      " << name() << ';' << nl
        << '}' << nl;

    writeDivider(os) << nl;

    return true;
}


bool Foam::IOobject::writeHeader(Ostream& os) const
{
    return writeHeader(os, type());
}


// ************************************************************************* //
