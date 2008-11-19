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

#include "TRIsurfaceFormat.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::TRIsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const pointField& pointLst,
    const Face& f,
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

template<class Face>
Foam::fileFormats::TRIsurfaceFormat<Face>::TRIsurfaceFormat
(
    const fileName& fName
)
{
    read(fName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::TRIsurfaceFormat<Face>::read
(
    const fileName& fName
)
{
    this->clear();

    IFstream is(fName);
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::TRIsurfaceFormat::read(const fileName&)"
        )
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
        string line = this->getLineNoComment(is);

        // handle continuations ?
        //          if (line[line.size()-1] == '\\')
        //          {
        //              line.substr(0, line.size()-1);
        //              line += this->getLineNoComment(is);
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

    // make our triangles directly
    List<Face>& faceLst = this->storedFaces();
    faceLst.setSize(regionLst.size());

    // transfer to normal list
    this->storedPoints().transfer(pointLst);
    this->storedRegions().transfer(regionLst);

    label ptI = 0;
    forAll(faceLst, faceI)
    {
        const label startPt = 3 * faceI;
        faceLst[faceI] = triFace(startPt, startPt+1, startPt+2);
    }

    this->setPatches(groupToPatch);
    this->stitchFaces(SMALL);
    return true;
}


template<class Face>
void Foam::fileFormats::TRIsurfaceFormat<Face>::write
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>& faceLst  = surf.faces();

    bool doSort = false;
    // a single region needs no sorting
    if (surf.patches().size() == 1)
    {
        doSort = false;
    }

    if (doSort)
    {
        labelList faceMap;
        List<surfGroup> patchLst = surf.sortedRegions(faceMap);

        label faceIndex = 0;
        forAll(patchLst, patchI)
        {
            forAll(patchLst[patchI], patchFaceI)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];
                writeShell(os, pointLst, f, patchI);
            }
        }
    }
    else
    {
        const List<label>& regionLst  = surf.regions();

        forAll(faceLst, faceI)
        {
            writeShell(os, pointLst, faceLst[faceI], regionLst[faceI]);
        }
    }
}


template<class Face>
void Foam::fileFormats::TRIsurfaceFormat<Face>::write
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>& faceLst  = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const Face& f = faceLst[faceIndex++];
            writeShell(os, pointLst, f, patchI);
        }
    }
}

// ************************************************************************* //
