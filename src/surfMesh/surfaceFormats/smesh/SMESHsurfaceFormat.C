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

#include "SMESHsurfaceFormat.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::SMESHsurfaceFormat<Face>::SMESHsurfaceFormat()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::SMESHsurfaceFormat<Face>::write
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    writeHeader(os, surf.points(), faceLst.size());

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const Face& f = faceLst[faceIndex++];

            os << f.size();
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
    const UnsortedMeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();

    writeHeader(os, surf.points(), faceLst.size());

    labelList faceMap;
    List<surfGroup> patchLst = surf.sortedRegions(faceMap);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const Face& f = faceLst[faceMap[faceIndex++]];

            os << f.size();
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
