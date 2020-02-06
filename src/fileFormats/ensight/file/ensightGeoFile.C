/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "ensightGeoFile.H"
#include "foamVersion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightGeoFile::initialize()
{
    writeBinaryHeader();

    // Description line 1
    write("Ensight Geometry File");
    newline();

    // Description line 2
    write(string("Written by OpenFOAM " + std::to_string(foamVersion::api)));
    newline();

    write("node id assign");
    newline();

    write("element id assign");
    newline();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightGeoFile::ensightGeoFile
(
    const fileName& pathname,
    IOstream::streamFormat format
)
:
    ensightFile(pathname, format)
{
    initialize();
}


Foam::ensightGeoFile::ensightGeoFile
(
    const fileName& path,
    const fileName& name,
    IOstream::streamFormat format
)
:
    ensightFile(path, name, format)
{
    initialize();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::ensightGeoFile::writeKeyword(const keyType& key)
{
    writeString(key);
    newline();

    return *this;
}


//
// Convenience Output Methods
//

void Foam::ensightGeoFile::beginPart
(
    const label index,
    const string& description
)
{
    beginPart(index);
    writeString(description);
    newline();
}


void Foam::ensightGeoFile::beginCoordinates(const label npoints)
{
    writeString(ensightFile::coordinates);
    newline();

    write(npoints);
    newline();
}


// ************************************************************************* //
