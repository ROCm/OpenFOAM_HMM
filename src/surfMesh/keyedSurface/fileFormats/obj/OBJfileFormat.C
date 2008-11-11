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

#include "OBJfileFormat.H"
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
    OBJfileFormat,
    fileExtension,
    obj
);

addNamedToMemberFunctionSelectionTable
(
    keyedSurface,
    OBJfileFormat,
    write,
    fileExtension,
    obj
);

addNamedToMemberFunctionSelectionTable
(
    meshedSurface,
    OBJfileFormat,
    write,
    fileExtension,
    obj
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::OBJfileFormat::writeHead
(
    Ostream& os,
    const pointField& pointLst,
    const List<FaceType>& faceLst,
    const List<surfGroup>& patchLst
)
{
    os  << "# Wavefront OBJ file written " << clock::dateTime().c_str() << nl
        << "o " << os.name().lessExt().name() << nl
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
        << "# <points count=\"" << pointLst.size() << "\">" << endl;

    // Write vertex coords
    forAll(pointLst, ptI)
    {
        os  << "v " << pointLst[ptI].x()
            << ' '  << pointLst[ptI].y()
            << ' '  << pointLst[ptI].z() << nl;
    }

    os  << "# </points>" << nl
        << nl
        << "# <faces count=\"" << faceLst.size() << "\">" << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::OBJfileFormat::OBJfileFormat()
:
    Foam::keyedSurface()
{}


Foam::fileFormats::OBJfileFormat::OBJfileFormat
(
    const fileName& fName,
    const bool triangulate
)
:
    Foam::keyedSurface()
{
    read(fName,triangulate);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::OBJfileFormat::read
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
            "fileFormats::OBJfileFormat::read(const fileName&)"
        )
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    DynamicList<point>    pointLst;
    DynamicList<FaceType> faceLst;
    DynamicList<label>    regionLst;
    HashTable<label>      groupToPatch;

    // leave faces that didn't have a group in 0
    label groupID = 0;
    label maxGroupID = -1;

    while (is.good())
    {
        string line = getLineNoComment(is);

        // handle continuations
        if (line[line.size()-1] == '\\')
        {
            line.substr(0, line.size()-1);
            line += getLineNoComment(is);
        }

        // Read first word
        IStringStream lineStream(line);
        word cmd;
        lineStream >> cmd;

        if (cmd == "v")
        {
            scalar x, y, z;
            lineStream >> x >> y >> z;
            pointLst.append(point(x, y, z));
        }
        else if (cmd == "g")
        {
            word groupName;
            lineStream >> groupName;

            HashTable<label>::const_iterator findGroup =
                groupToPatch.find(groupName);

            if (findGroup != groupToPatch.end())
            {
                groupID = findGroup();
            }
            else
            {
                // special treatment if any initial faces were not in a group
                if (maxGroupID == -1 && faceLst.size())
                {
                    groupToPatch.insert("patch0", 0);
                    maxGroupID = 0;
                }
                groupID = ++maxGroupID;
                groupToPatch.insert(groupName, groupID);
            }
        }
        else if (cmd == "f")
        {
            DynamicList<label> dynVertices;

            // Assume 'f' is followed by space.
            string::size_type endNum = 1;

            while (true)
            {
                string::size_type startNum =
                    line.find_first_not_of(' ', endNum);

                if (startNum == string::size_type(string::npos))
                {
                    break;
                }

                endNum = line.find(' ', startNum);

                string vertexSpec;
                if (endNum != string::size_type(string::npos))
                {
                    vertexSpec = line.substr(startNum, endNum-startNum);
                }
                else
                {
                    vertexSpec = line.substr(startNum, line.size() - startNum);
                }

                string::size_type slashPos = vertexSpec.find('/');

                label vertI = 0;
                if (slashPos != string::size_type(string::npos))
                {
                    IStringStream intStream(vertexSpec.substr(0, slashPos));

                    intStream >> vertI;
                }
                else
                {
                    IStringStream intStream(vertexSpec);

                    intStream >> vertI;
                }
                dynVertices.append(vertI - 1);
            }
            dynVertices.shrink();

            FaceType f(SubList<label>(dynVertices, dynVertices.size()));

            if (triangulate && f.size() > 3)
            {
                triFace fTri;

                // simple face triangulation about f[0].
                // Cannot use face::triangulation since points are incomplete
                fTri[0] = f[0];
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    label fp2 = (fp1 + 1) % f.size();

                    fTri[1] = f[fp1];
                    fTri[2] = f[fp2];

                    faceLst.append(fTri);
                    regionLst.append(groupID);
                }
            }
            else
            {
                faceLst.append(f);
                regionLst.append(groupID);
            }
        }
    }

    // transfer to normal lists
    points().transfer(pointLst);
    faces().transfer(faceLst);
    regions().transfer(regionLst);

    setPatches(groupToPatch);

    return true;
}


void Foam::fileFormats::OBJfileFormat::write
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
        // Print all faces belonging to this region
        const surfGroup& patch = patchLst[patchI];

        os << "g " << patch.name() << endl;

        forAll(patch, patchFaceI)
        {
            const face& f = faceLst[faceMap[faceIndex++]];

            os  << 'f';
            forAll(f, fp)
            {
                os << ' ' << f[fp] + 1;
            }
            os << endl;
        }
    }

    os  << "# </faces>" << endl;
}


void Foam::fileFormats::OBJfileFormat::write
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
        const surfGroup& patch = patchLst[patchI];

        os << "g " << patch.name() << endl;

        forAll(patch, patchFaceI)
        {
            const face& f = faceLst[faceIndex++];

            os  << 'f';
            forAll(f, fp)
            {
                os << ' ' << f[fp] + 1;
            }
            os << endl;
        }
    }
    os  << "# </faces>" << endl;
}

// ************************************************************************* //
