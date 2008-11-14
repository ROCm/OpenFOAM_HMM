/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "STARCDsurfaceFormatCore.H"
#include "clock.H"
#include "OSspecific.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//! @cond localscope
const int starcdShellShape = 3;
const int starcdShellType  = 4;
//! @endcond localscope


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileFormats::STARCDsurfaceFormatCore::readHeader
(
    IFstream& is,
    const word& signature
)
{
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::STARCDsurfaceFormatCore::readHeader(...)"
        )
            << "cannot read " << signature  << "  " << is.name()
            << abort(FatalError);
    }

    word header;
    label majorVersion;

    string line;

    is.getLine(line);
    IStringStream(line)() >> header;

    is.getLine(line);
    IStringStream(line)() >> majorVersion;

    // add other checks ...
    if (header != signature)
    {
        Info<< "header mismatch " << signature << "  " << is.name()
            << endl;
    }

    return true;
}


void Foam::fileFormats::STARCDsurfaceFormatCore::writeHeader
(
    Ostream& os,
    const char* filetype
)
{
    os  << "PROSTAR_" << filetype << nl
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


void Foam::fileFormats::STARCDsurfaceFormatCore::writePoints
(
    Ostream& os,
    const pointField& pointLst
)
{
    writeHeader(os, "VERTEX");

    // Set the precision of the points data to 10
    os.precision(10);

    // force decimal point for Fortran input
    os.setf(std::ios::showpoint);

    forAll(pointLst, ptI)
    {
        os
            << ptI + 1 << " "
            << pointLst[ptI].x() << " "
            << pointLst[ptI].y() << " "
            << pointLst[ptI].z() << nl;
    }
    os.flush();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
