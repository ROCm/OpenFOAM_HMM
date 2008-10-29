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
    IFstream is(fName);

    if (!is.good())
    {
        FatalErrorIn("fileFormats::OFFfileFormat(const fileName&)")
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    // Read header
    string hdr = getLineNoComment(is);
    if (hdr != "OFF")
    {
        FatalErrorIn("fileFormats::OFFfileFormat(const fileName&)")
            << "OFF file " << fName
            << " does not start with 'OFF'"
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
    DynamicList<keyedFace> faceLst(nElems);

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
                face fTri(3);

                // simple face triangulation about f[0].
                // cannot use face::triangulation since points are incomplete
                fTri[0] = f[0];
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    label fp2 = (fp1 + 1) % f.size();

                    fTri[1] = f[fp1];
                    fTri[2] = f[fp2];

                    faceLst.append(keyedFace(fTri, 0));
                }
            }
            else
            {
                faceLst.append(keyedFace(f, 0));
            }
        }
    }

    // no region information
    points().transfer(pointLst);
    faces().transfer(faceLst);
    setPatches(0);
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fileFormats::OFFfileFormat::write
(
    Ostream& os,
    const keyedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<keyedFace>& faceLst = surf.faces();

    labelList faceMap;
    List<surfacePatch> patchLst = surf.sortedRegions(faceMap);

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
    const pointField& pointLst = surf.points();
    const List<face>& faceLst = surf.faces();
    const List<surfacePatch>& patchLst = surf.patches();

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


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// ************************************************************************* //
