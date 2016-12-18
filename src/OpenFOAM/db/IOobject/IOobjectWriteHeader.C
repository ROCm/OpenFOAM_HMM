/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
#include "endian.H"
#include "label.H"
#include "scalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// file-scope
// Hint about machine endian, OpenFOAM label and scalar sizes
static const std::string archHint =
(
#ifdef WM_LITTLE_ENDIAN
    "LSB"
#elif defined (WM_BIG_ENDIAN)
    "MSB"
#else
    "???"
#endif
    ";label="  + std::to_string(8*sizeof(Foam::label))
  + ";scalar=" + std::to_string(8*sizeof(Foam::scalar))
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        os  << "    arch        " << archHint << ";\n";
    }

    if (!note().empty())
    {
        os  << "    note        " << note() << ";\n";
    }

    os  << "    location    " << instance()/db().dbDir()/local() << ";\n"
        << "    object      " << name() << ";\n"
        << "}" << nl;

    writeDivider(os) << endl;

    return true;
}


bool Foam::IOobject::writeHeader(Ostream& os) const
{
    return writeHeader(os, type());
}


// ************************************************************************* //
