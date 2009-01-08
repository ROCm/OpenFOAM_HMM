/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "SMESHsurfaceFormatCore.H"
#include "clock.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::SMESHsurfaceFormatCore::writeHeader
(
    Ostream& os,
    const pointField& pointLst,
    const label nFaces
)
{
    // Write header
    os << "# tetgen .smesh file written " << clock::dateTime().c_str() << nl;
    os << "# <points count=\"" << pointLst.size() << "\">" << endl;
    os << pointLst.size() << " 3" << nl;    // 3: dimensions

    // Write vertex coords
    forAll(pointLst, ptI)
    {
        os  << ptI
            << ' ' << pointLst[ptI].x()
            << ' ' << pointLst[ptI].y()
            << ' ' << pointLst[ptI].z() << nl;
    }
    os  << "# </points>" << nl
        << nl
        << "# <faces count=\"" << nFaces << "\">" << endl;

    os  << nFaces << " 1" << endl;   // one attribute: region number
}


void Foam::fileFormats::SMESHsurfaceFormatCore::writeTail(Ostream& os)
{
    os  << "# </faces>" << nl
        << nl
        << "# no holes or regions:" << nl
        << '0' << nl        // holes
        << '0' << endl;     // regions
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
