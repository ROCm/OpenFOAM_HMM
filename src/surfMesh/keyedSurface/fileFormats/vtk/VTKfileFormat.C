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

#include "VTKfileFormat.H"
#include "clock.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

addNamedToMemberFunctionSelectionTable
(
    keyedSurface,
    VTKfileFormat,
    write,
    fileExtension,
    vtk
);

addNamedToMemberFunctionSelectionTable
(
    meshedSurface,
    VTKfileFormat,
    write,
    fileExtension,
    vtk
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::VTKfileFormat::writeHead
(
    Ostream& os,
    const pointField& pointLst,
    const List<face>& faceLst
)
{
    // Write header
    os  << "# vtk DataFile Version 2.0" << nl
        << "surface written " << clock::dateTime().c_str() << nl
        << "ASCII" << nl
        << nl
        << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << pointLst.size() << " float" << nl;
    forAll(pointLst, ptI)
    {
        os  << pointLst[ptI].x() << ' '
            << pointLst[ptI].y() << ' '
            << pointLst[ptI].z() << nl;
    }

    label nNodes = 0;
    forAll(faceLst, faceI)
    {
        nNodes += faceLst[faceI].size();
    }

    os  << nl
        << "POLYGONS " << faceLst.size() << ' '
        << faceLst.size() + nNodes << nl;
}


void Foam::fileFormats::VTKfileFormat::writeTail
(
    Ostream& os,
    const List<surfacePatch>& patchLst
)
{
    label nFaces = 0;
    forAll(patchLst, patchI)
    {
        nFaces += patchLst[patchI].size();
    }

    // Print region numbers
    os  << nl
        << "CELL_DATA " << nFaces << nl
        << "FIELD attributes 1" << nl
        << "region 1 " << nFaces << " float" << nl;


    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            if (patchFaceI)
            {
                if ((patchFaceI % 20) == 0)
                {
                    os  << nl;
                }
                else
                {
                    os  << ' ';
                }
            }
            os  << patchI + 1;
        }
        os  << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::VTKfileFormat::VTKfileFormat()
:
    Foam::keyedSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fileFormats::VTKfileFormat::write
(
    Ostream& os,
    const keyedSurface& surf
)
{
    const List<face>& faceLst = surf.faces();

    writeHead(os, surf.points(), faceLst);

    labelList faceMap;
    List<surfacePatch> patchLst = surf.sortedRegions(faceMap);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const face& f = faceLst[faceMap[faceIndex++]];

            os  << f.size();
            forAll(f, fp)
            {
                os  << ' ' << f[fp];
            }
            os  << ' ' << nl;
        }
    }

    // Print region numbers
    writeTail(os, patchLst);
}


void Foam::fileFormats::VTKfileFormat::write
(
    Ostream& os,
    const meshedSurface& surf
)
{
    const List<face>& faceLst = surf.faces();
    const List<surfacePatch>& patchLst = surf.patches();

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
                os  << ' ' << f[fp];
            }
            os  << ' ' << nl;
        }
    }

    // Print region numbers
    writeTail(os, patchLst);
}

// ************************************************************************* //
