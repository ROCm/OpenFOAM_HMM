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

#include "SMESHsurfaceFormat.H"
#include "clock.H"
#include "IStringStream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::SMESHsurfaceFormat<Face>::writeHead
(
    Ostream& os,
    const pointField& pointLst,
    const List<Face>& faceLst
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
        << "# <faces count=\"" << faceLst.size() << "\">" << endl;

    os  << faceLst.size() << " 1" << endl;   // one attribute: region number
}


template<class Face>
void Foam::fileFormats::SMESHsurfaceFormat<Face>::writeTail(Ostream& os)
{
    os  << "# </faces>" << nl
        << nl
        << "# no holes or regions:" << nl
        << '0' << nl        // holes
        << '0' << endl;     // regions
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::SMESHsurfaceFormat<Face>::SMESHsurfaceFormat()
:
    ParentType()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::SMESHsurfaceFormat<Face>::write
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();

    writeHead(os, surf.points(), faceLst);

    labelList faceMap;
    List<surfGroup> patchLst = surf.sortedRegions(faceMap);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const Face& f = faceLst[faceMap[faceIndex++]];

            os  << f.size();
            forAll(f, fp)
            {
                os << ' ' << f[fp];
            }
            os << ' ' << patchI << endl;
        }
    }

    writeTail(os);
}


template<class Face>
void Foam::fileFormats::SMESHsurfaceFormat<Face>::write
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    writeHead(os, surf.points(), faceLst);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const face& f = faceLst[faceIndex++];

            os  << f.size();
            forAll(f, fp)
            {
                os << ' ' << f[fp];
            }
            os << ' ' << patchI << endl;
        }
    }

    writeTail(os);
}

// ************************************************************************* //
