/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Description
    Writes the header description of the File to the stream
    associated with the File.

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
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

Foam::Ostream& Foam::IOobject::writeBanner(Ostream& os, bool noHint)
{
    // The version padded with spaces to fit after "Version:  "
    // - initialized with zero-length string to detect if it has been populated
    static char paddedVersion[39] = "";

    if (!*paddedVersion)
    {
        // Populate: like strncpy but without trailing '\0'
        const char *p = Foam::FOAMversion;

        memset(paddedVersion, ' ', 38);
        for (int i = 0; *p && i < 38; ++i)
        {
            paddedVersion[i] = *p++;
        }
        paddedVersion[38] = '\0';
    }

    os  <<
        "/*--------------------------------";

    if (noHint)
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
        " Web:      www.OpenFOAM.com                      |\n"
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


bool Foam::IOobject::writeHeader(Ostream& os, const word& type) const
{
    if (!os.good())
    {
        InfoInFunction
            << "No stream open for write" << nl
            << os.info() << endl;

        return false;
    }

    writeBanner(os)
        << "FoamFile\n{\n"
        << "    version     " << os.version() << ";\n"
        << "    format      " << os.format() << ";\n"
        << "    class       " << type << ";\n";

    if (os.format() == IOstream::BINARY)
    {
        os  << "    arch        " << Foam::FOAMbuildArch << ";\n";
    }

    if (!note().empty())
    {
        os  << "    note        " << note() << ";\n";
    }

    os  << "    location    " << instance()/db().dbDir()/local() << ";\n"
        << "    object      " << name() << ";\n"
        << "}" << nl;

    writeDivider(os) << nl;

    return true;
}


bool Foam::IOobject::writeHeader(Ostream& os) const
{
    return writeHeader(os, type());
}


// ************************************************************************* //
