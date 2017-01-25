/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "STLReader.H"
#include "Map.H"
#include "IFstream.H"

#undef DEBUG_STLBINARY

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileFormats::STLReader::readBINARY
(
    const fileName& filename
)
{
    sorted_ = true;
    format_ = UNKNOWN;

    label nTris = 0;
    autoPtr<istream> streamPtr = readBinaryHeader(filename, nTris);

    if (!streamPtr.valid())
    {
        FatalErrorInFunction
            << "Error reading file " << filename
            << " or file " << filename + ".gz"
            << exit(FatalError);
    }

    istream& is = streamPtr();

#ifdef DEBUG_STLBINARY
    Info<< "# " << nTris << " facets" << endl;
    label prevZone = -1;
#endif

    points_.setSize(3*nTris);
    zoneIds_.setSize(nTris);

    Map<label> lookup;
    DynamicList<label> dynSizes;

    label ptI = 0;
    label zoneI = -1;
    forAll(zoneIds_, facei)
    {
        // Read STL triangle
        STLtriangle stlTri(is);

        // transcribe the vertices of the STL triangle -> points
        points_[ptI++] = stlTri.a();
        points_[ptI++] = stlTri.b();
        points_[ptI++] = stlTri.c();

        // interpret STL attribute as a zone
        const label origId = stlTri.attrib();

        Map<label>::const_iterator fnd = lookup.find(origId);
        if (fnd != lookup.end())
        {
            if (zoneI != fnd())
            {
                // group appeared out of order
                sorted_ = false;
            }
            zoneI = fnd();
        }
        else
        {
            zoneI = dynSizes.size();
            lookup.insert(origId, zoneI);
            dynSizes.append(0);
        }

        zoneIds_[facei] = zoneI;
        dynSizes[zoneI]++;

#ifdef DEBUG_STLBINARY
        if (prevZone != zoneI)
        {
            if (prevZone != -1)
            {
                Info<< "endsolid zone" << prevZone << nl;
            }
            prevZone = zoneI;

            Info<< "solid zone" << prevZone << nl;
        }

        stlTri.print(Info);
#endif
    }

#ifdef DEBUG_STLBINARY
    if (prevZone != -1)
    {
        Info<< "endsolid zone" << prevZone << nl;
    }
#endif

    names_.clear();
    sizes_.transfer(dynSizes);

    format_ = BINARY;
    return true;
}


bool Foam::fileFormats::STLReader::readFile
(
    const fileName& filename,
    const STLFormat& format
)
{
    if (format == UNKNOWN ? detectBinaryHeader(filename) : format == BINARY)
    {
        return readBINARY(filename);
    }
    else
    {
        return readASCII(filename);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::STLReader::STLReader
(
    const fileName& filename
)
:
    sorted_(true),
    points_(),
    zoneIds_(),
    names_(),
    sizes_(),
    format_(STLCore::UNKNOWN)
{
    // Auto-detect ASCII/BINARY format
    readFile(filename, STLCore::UNKNOWN);
}


Foam::fileFormats::STLReader::STLReader
(
    const fileName& filename,
    const STLFormat& format
)
:
    sorted_(true),
    points_(),
    zoneIds_(),
    names_(),
    sizes_(),
    format_(STLCore::UNKNOWN)
{
    // Manually specified ASCII/BINARY format
    readFile(filename, format);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileFormats::STLReader::~STLReader()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fileFormats::STLReader::clear()
{
    sorted_ = true;
    points_.clear();
    zoneIds_.clear();
    names_.clear();
    sizes_.clear();
    format_ = UNKNOWN;
}


// ************************************************************************* //
