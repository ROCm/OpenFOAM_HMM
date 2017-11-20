/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "STARCDCore.H"
#include "ListOps.H"
#include "clock.H"
#include "PackedBoolList.H"
#include "DynamicList.H"
#include "StringStream.H"
#include "OSspecific.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fileFormats::STARCDCore::fileHeader
>
Foam::fileFormats::STARCDCore::fileHeaders_
{
    { fileHeader::HEADER_CEL, "PROSTAR_CELL" },
    { fileHeader::HEADER_VRT, "PROSTAR_VERTEX" },
    { fileHeader::HEADER_BND, "PROSTAR_BOUNDARY" }
};


const Foam::Enum
<
    Foam::fileFormats::STARCDCore::fileExt
>
Foam::fileFormats::STARCDCore::fileExtensions_
{
    { fileExt::CEL_FILE, "cel" },
    { fileExt::VRT_FILE, "vrt" },
    { fileExt::BND_FILE, "bnd" },
    { fileExt::INP_FILE, "inp" }
};


const char* const Foam::fileFormats::STARCDCore::defaultBoundaryName =
    "Default_Boundary_Region";

const char* const Foam::fileFormats::STARCDCore::defaultSolidBoundaryName =
    "Default_Boundary_Solid";


const Foam::Map<Foam::FixedList<int, 6>>
Foam::fileFormats::STARCDCore::foamToStarFaceAddr =
{
    { starcdHex,   { 4, 5, 2, 3, 0, 1 } },
    { starcdPrism, { 0, 1, 4, 5, 2, -1 } },
    { starcdTet,   { 5, 4, 2, 0, -1, -1 } },
    { starcdPyr,   { 0, 4, 3, 5, 2, -1 } }
};


const Foam::Map<Foam::FixedList<int, 6>>
Foam::fileFormats::STARCDCore::starToFoamFaceAddr =
{
    { starcdHex,   { 4, 5, 2, 3, 0, 1 } },
    { starcdPrism, { 0, 1, 4, -1, 2, 3 } },
    { starcdTet,   { 3, -1, 2, -1, 1, 0 } },
    { starcdPyr,   { 0, -1, 4, 2, 1, 3 } }
};


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

bool Foam::fileFormats::STARCDCore::readHeader
(
    IFstream& is,
    const enum fileHeader header
)
{
    if (!is.good())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    word  magic;
    label majorVersion;

    string line;

    is.getLine(line);
    IStringStream(line)() >> magic;

    is.getLine(line);
    IStringStream(line)() >> majorVersion;

    // add other checks ...
    if (magic != fileHeaders_[header])
    {
        Info<< "header mismatch " << fileHeaders_[header]
            << "  " << is.name()
            << endl;

        return false;
    }

    return true;
}


void Foam::fileFormats::STARCDCore::writeHeader
(
    Ostream& os,
    const enum fileHeader header
)
{
    os  << fileHeaders_[header] << nl
        << 4000
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileFormats::STARCDCore::starFileName
(
    const fileName& base,
    const enum fileExt ext
)
{
    return base + '.' + fileExtensions_[ext];
}


void Foam::fileFormats::STARCDCore::removeFiles(const fileName& base)
{
    Foam::rm(starFileName(base, VRT_FILE));
    Foam::rm(starFileName(base, CEL_FILE));
    Foam::rm(starFileName(base, BND_FILE));
    Foam::rm(starFileName(base, INP_FILE));
}


Foam::label Foam::fileFormats::STARCDCore::readPoints
(
    IFstream& is,
    List<point>& points,
    List<label>& ids
)
{
    label maxId = 0;

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    readHeader(is, HEADER_VRT);

    // reuse memory if possible
    DynamicList<point> dynPoints(points.xfer());
    DynamicList<label> dynPointId(ids.xfer());    // STAR-CD index of points

    dynPoints.clear();
    dynPointId.clear();

    {
        label lineLabel;
        scalar x, y, z;

        while ((is >> lineLabel).good())
        {
            maxId = max(maxId, lineLabel);
            is >> x >> y >> z;

            dynPoints.append(point(x, y, z));
            dynPointId.append(lineLabel);
        }
    }

    points.transfer(dynPoints);
    ids.transfer(dynPointId);

    return maxId;
}


void Foam::fileFormats::STARCDCore::writePoints
(
    Ostream& os,
    const UList<point>& points,
    const scalar scaleFactor
)
{
    writeHeader(os, HEADER_VRT);

    // Set the precision of the points data to 10
    os.precision(10);

    // force decimal point for Fortran input
    os.setf(std::ios::showpoint);

    forAll(points, ptI)
    {
        // convert [m] -> [mm] etc
        os
            << ptI + 1 << ' '
            << scaleFactor * points[ptI].x() << ' '
            << scaleFactor * points[ptI].y() << ' '
            << scaleFactor * points[ptI].z() << '\n';
    }
    os.flush();
}


// ************************************************************************* //
