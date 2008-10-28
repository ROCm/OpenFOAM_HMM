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

#include "TRIfileFormat.H"
#include "clock.H"
#include "IFstream.H"
#include "IOmanip.H"
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
    TRIfileFormat,
    fileExtension,
    tri
);

addNamedToMemberFunctionSelectionTable
(
    keyedSurface,
    TRIfileFormat,
    write,
    fileExtension,
    tri
);

addNamedToMemberFunctionSelectionTable
(
    meshedSurface,
    TRIfileFormat,
    write,
    fileExtension,
    tri
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::fileFormats::TRIfileFormat::writeShell
(
    Ostream& os,
    const pointField& pointLst,
    const face& f,
    const label patchI
)
{
    // simple triangulation about f[0].
    // better triangulation should have been done before
    const point& p0 = pointLst[f[0]];
    for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
    {
        label fp2 = (fp1 + 1) % f.size();

        const point& p1 = pointLst[f[fp1]];
        const point& p2 = pointLst[f[fp2]];

        os  << p0.x() << ' ' << p0.y() << ' ' << p0.z() << ' '
            << p1.x() << ' ' << p1.y() << ' ' << p1.z() << ' '
            << p2.x() << ' ' << p2.y() << ' ' << p2.z() << ' '
            // region as colour
            << "0x" << hex << patchI << dec << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::TRIfileFormat::TRIfileFormat()
:
    Foam::keyedSurface()
{}


Foam::fileFormats::TRIfileFormat::TRIfileFormat
(
    const fileName& fName,
    const bool triangulate
)
:
    Foam::keyedSurface()
{
    IFstream is(fName);

    if (!is.good())
    {
        FatalErrorIn("fileFormats::TRIfileFormat(const fileName&)")
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    // uses similar structure as STL, just some points
    DynamicList<point> pointLst;
    DynamicList<label> regionLst;
    HashTable<label>   groupToPatch;

    // leave faces that didn't have a group in 0
    label nUngrouped = 0;
    label groupID = 0;
    label maxGroupID = -1;

    while (is.good())
    {
        string line = getLineNoComment(is);

        // handle continuations ?
        //          if (line[line.size()-1] == '\\')
        //          {
        //              line.substr(0, line.size()-1);
        //              line += getLineNoComment(is);
        //          }

        IStringStream lineStream(line);

        point p
        (
            readScalar(lineStream),
            readScalar(lineStream),
            readScalar(lineStream)
        );

        if (!lineStream) break;

        pointLst.append(p);
        pointLst.append
        (
            point
            (
                readScalar(lineStream),
                readScalar(lineStream),
                readScalar(lineStream)
            )
        );
        pointLst.append
        (
            point
            (
                readScalar(lineStream),
                readScalar(lineStream),
                readScalar(lineStream)
            )
        );

        // Region/colour in .tri file starts with 0x. Skip.
        // ie, instead of having 0xFF, skip 0 and leave xFF to
        // get read as a word and name it "patchFF"

        char zero;
        lineStream >> zero;

        word rawName(lineStream);
        word groupName("patch" + rawName(1, rawName.size()-1));

        HashTable<label>::const_iterator findGroup =
            groupToPatch.find(groupName);

        if (findGroup != groupToPatch.end())
        {
            groupID = findGroup();
        }
        else
        {
            // special treatment if any initial faces were not in a group
            if (maxGroupID == -1 && regionLst.size())
            {
                groupToPatch.insert("patch0", 0);
                nUngrouped = regionLst.size();
                maxGroupID = 0;
            }
            groupID = ++maxGroupID;
            groupToPatch.insert(groupName, groupID);
        }

        regionLst.append(groupID);
    }

    // transfer to normal list
    points().transfer(pointLst);
    regions().transfer(regionLst);

    // make our triangles directly
    List<face>& faceLst = faces();
    faceLst.setSize(regions().size());

    label ptI = 0;
    forAll(faceLst, faceI)
    {
        face& fTri = faceLst[faceI];
        fTri.setSize(3);

        fTri[0] = ptI++;
        fTri[1] = ptI++;
        fTri[2] = ptI++;
    }

    setPatches(groupToPatch);
    stitchFaces(SMALL);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fileFormats::TRIfileFormat::write
(
    Ostream& os,
    const keyedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<face>& faceLst  = surf.faces();

    labelList faceMap;
    List<surfacePatch> patchLst = surf.sortedRegions(faceMap);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const face& f = faceLst[faceMap[faceIndex++]];
            writeShell(os, pointLst, f, patchI);
        }
    }
}


void Foam::fileFormats::TRIfileFormat::write
(
    Ostream& os,
    const meshedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<face>& faceLst  = surf.faces();
    const List<surfacePatch>& patchLst = surf.patches();

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const face& f = faceLst[faceIndex++];
            writeShell(os, pointLst, f, patchI);
        }
    }
}

// ************************************************************************* //
