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

#include "OBJsurfaceFormat.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::OBJsurfaceFormat<Face>::OBJsurfaceFormat
(
    const fileName& fName
)
:
    ParentType()
{
    read(fName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::OBJsurfaceFormat<Face>::read
(
    const fileName& fName
)
{
    ParentType::clear();
    const bool mustTriangulate = ParentType::isTri();

    IFstream is(fName);
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::OBJsurfaceFormat::read(const fileName&)"
        )
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    DynamicList<point>    pointLst;
    DynamicList<Face>     faceLst;
    DynamicList<label>    regionLst;
    HashTable<label>      groupToPatch;

    // leave faces that didn't have a group in 0
    label groupID = 0;
    label maxGroupID = -1;

    while (is.good())
    {
        string line = ParentType::getLineNoComment(is);

        // handle continuations
        if (line[line.size()-1] == '\\')
        {
            line.substr(0, line.size()-1);
            line += ParentType::getLineNoComment(is);
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

            UList<label>& f = static_cast<UList<label>&>(dynVertices);

            if (mustTriangulate && f.size() > 3)
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
                faceLst.append(Face(f));
                regionLst.append(groupID);
            }
        }
    }

    // transfer to normal lists
    ParentType::storedPoints().transfer(pointLst);
    ParentType::storedFaces().transfer(faceLst);
    ParentType::storedRegions().transfer(regionLst);

    ParentType::setPatches(groupToPatch);
    return true;
}


template<class Face>
void Foam::fileFormats::OBJsurfaceFormat<Face>::write
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();

    labelList faceMap;
    List<surfGroup> patchLst = surf.sortedRegions(faceMap);

    writeHeader(os, surf.points(), faceLst.size(), patchLst);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        // Print all faces belonging to this region
        const surfGroup& patch = patchLst[patchI];

        os << "g " << patch.name() << endl;

        forAll(patch, patchFaceI)
        {
            const Face& f = faceLst[faceMap[faceIndex++]];

            os << 'f';
            forAll(f, fp)
            {
                os << ' ' << f[fp] + 1;
            }
            os << endl;
        }
    }

    os << "# </faces>" << endl;
}


template<class Face>
void Foam::fileFormats::OBJsurfaceFormat<Face>::write
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    writeHeader(os, surf.points(), faceLst.size(), patchLst);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        const surfGroup& patch = patchLst[patchI];

        os << "g " << patch.name() << endl;

        forAll(patch, patchFaceI)
        {
            const Face& f = faceLst[faceIndex++];

            os << 'f';
            forAll(f, fp)
            {
                os << ' ' << f[fp] + 1;
            }
            os << endl;
        }
    }
    os << "# </faces>" << endl;
}

// ************************************************************************* //
