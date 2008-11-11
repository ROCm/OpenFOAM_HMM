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

#include "OFFfileFormat.H"
#include "triFace.H"
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

addNamedToRunTimeSelectionTable
(
    keyedSurface,
    OFFfileFormat,
    fileExtension,
    off
);

addNamedToMemberFunctionSelectionTable
(
    keyedSurface,
    OFFfileFormat,
    write,
    fileExtension,
    off
);

addNamedToMemberFunctionSelectionTable
(
    meshedSurface,
    OFFfileFormat,
    write,
    fileExtension,
    off
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::OFFfileFormat::writeHead
(
    Ostream& os,
    const pointField& pointLst,
    const List<FaceType>& faceLst,
    const List<surfGroup>& patchLst
)
{
    // Write header
    os  << "OFF" << endl
        << "# Geomview OFF file written " << clock::dateTime().c_str() << nl
        << nl
        << "# points : " << pointLst.size() << nl
        << "# faces  : " << faceLst.size() << nl
        << "# patches: " << patchLst.size() << nl;

    // Print patch names as comment
    forAll(patchLst, patchI)
    {
        os  << "#   " << patchI << "  " << patchLst[patchI].name()
            << "  (nFaces: " << patchLst[patchI].size() << ")" << nl;
    }

    os  << nl
        << "# nPoints  nFaces  nEdges" << nl
        << pointLst.size() << ' ' << faceLst.size() << ' ' << 0 << nl;

    os  << nl
        << "# <points count=\"" << pointLst.size() << "\">" << endl;

    // Write vertex coords
    forAll(pointLst, ptI)
    {
        os  << pointLst[ptI].x() << ' '
            << pointLst[ptI].y() << ' '
            << pointLst[ptI].z() << " #" << ptI << endl;
    }

    os  << "# </points>" << nl
        << nl
        << "# <faces count=\"" << faceLst.size() << "\">" << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::OFFfileFormat::OFFfileFormat()
:
    keyedSurface()
{}


Foam::fileFormats::OFFfileFormat::OFFfileFormat
(
    const fileName& fName,
    const bool triangulate
)
:
    keyedSurface()
{
    read(fName,triangulate);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::OFFfileFormat::read
(
    const fileName& fName,
    const bool triangulate
)
{
    clear();
    IFstream is(fName);

    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::OFFfileFormat::read(const fileName&)"
        )
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    // Read header
    string hdr = getLineNoComment(is);
    if (hdr != "OFF")
    {
        FatalErrorIn
        (
            "fileFormats::OFFfileFormat::read(const fileName&)"
        )
            << "OFF file " << fName << " does not start with 'OFF'"
            << exit(FatalError);
    }


    // get dimensions
    label nPoints, nEdges, nElems;

    string line = getLineNoComment(is);
    {
        IStringStream lineStream(line);
        lineStream >> nPoints >> nElems >> nEdges;
    }

    // Read points
    pointField pointLst(nPoints);
    forAll(pointLst, pointI)
    {
        scalar x, y, z;
        line = getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream >> x >> y >> z;
        }
        pointLst[pointI] = point(x, y, z);
    }

    // Read faces - ignore optional region information
    // use a DynamicList for possible on-the-fly triangulation
    DynamicList<FaceType>  faceLst(nElems);

    forAll(faceLst, faceI)
    {
        line = getLineNoComment(is);
        {
            IStringStream lineStream(line);

            label nVerts;
            lineStream >> nVerts;

            face f(nVerts);

            forAll(f, fp)
            {
                lineStream >> f[fp];
            }

            if (triangulate && f.size() > 3)
            {
                triFace fTri;

                // simple face triangulation about f[0].
                // cannot use face::triangulation since points are incomplete
                fTri[0] = f[0];
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    label fp2 = (fp1 + 1) % f.size();

                    fTri[1] = f[fp1];
                    fTri[2] = f[fp2];

                    faceLst.append(fTri);
                }
            }
            else
            {
                faceLst.append(f);
            }
        }
    }

    // transfer to normal lists
    points().transfer(pointLst);
    faces().transfer(faceLst);

    // no region information
    regions().setSize(size());
    regions() = 0;

    setPatches(0);
    return true;
}


void Foam::fileFormats::OFFfileFormat::write
(
    Ostream& os,
    const keyedSurface& surf
)
{
    const List<FaceType>& faceLst = surf.faces();

    labelList faceMap;
    List<surfGroup> patchLst = surf.sortedRegions(faceMap);

    writeHead(os, surf.points(), faceLst, patchLst);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        os << "# <patch name=\"" << patchLst[patchI].name() << "\">" << endl;

        forAll(patchLst[patchI], patchFaceI)
        {
            const face& f = faceLst[faceMap[faceIndex++]];

            os  << f.size();
            forAll(f, fp)
            {
                os << ' ' << f[fp];
            }

            // add optional region information
            os << ' ' << patchI << endl;
        }
        os << "# </patch>" << endl;
    }
    os  << "# </faces>" << endl;
}


void Foam::fileFormats::OFFfileFormat::write
(
    Ostream& os,
    const meshedSurface& surf
)
{
    const List<FaceType>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    writeHead(os, surf.points(), faceLst, patchLst);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        os << "# <patch name=\"" << patchLst[patchI].name() << "\">" << endl;

        forAll(patchLst[patchI], patchFaceI)
        {
            const face& f = faceLst[faceIndex++];

            os  << f.size();
            forAll(f, fp)
            {
                os << ' ' << f[fp];
            }

            // add optional region information
            os << ' ' << patchI << endl;
        }
        os << "# </patch>" << endl;
    }
    os  << "# </faces>" << endl;
}

// ************************************************************************* //
