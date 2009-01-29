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

#include "OBJsurfaceFormat.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::OBJsurfaceFormat<Face>::OBJsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::OBJsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    const bool mustTriangulate = this->isTri();
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::OBJsurfaceFormat::read(const fileName&)"
        )
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // assume that the groups are not intermixed
    bool sorted = true;

    DynamicList<point> dynPoints;
    DynamicList<Face>  dynFaces;
    DynamicList<label> dynRegions;
    DynamicList<word>  dynNames;
    DynamicList<label> dynSizes;
    HashTable<label>   lookup;

    // place faces without a group in region0
    label regionI = 0;
    lookup.insert("region0", regionI);
    dynNames.append("region0");
    dynSizes.append(0);

    while (is.good())
    {
        string line = this->getLineNoComment(is);

        // handle continuations
        if (line[line.size()-1] == '\\')
        {
            line.substr(0, line.size()-1);
            line += this->getLineNoComment(is);
        }

        // Read first word
        IStringStream lineStream(line);
        word cmd;
        lineStream >> cmd;

        if (cmd == "v")
        {
            scalar x, y, z;
            lineStream >> x >> y >> z;
            dynPoints.append(point(x, y, z));
        }
        else if (cmd == "g")
        {
            word name;
            lineStream >> name;

            HashTable<label>::const_iterator fnd = lookup.find(name);
            if (fnd != lookup.end())
            {
                if (regionI != fnd())
                {
                    // group appeared out of order
                    sorted = false;
                }
                regionI = fnd();
            }
            else
            {
                regionI = dynSizes.size();
                lookup.insert(name, regionI);
                dynNames.append(name);
                dynSizes.append(0);
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

                if (startNum == string::npos)
                {
                    break;
                }

                endNum = line.find(' ', startNum);

                string vertexSpec;
                if (endNum != string::npos)
                {
                    vertexSpec = line.substr(startNum, endNum-startNum);
                }
                else
                {
                    vertexSpec = line.substr(startNum, line.size() - startNum);
                }

                string::size_type slashPos = vertexSpec.find('/');

                label vertI = 0;
                if (slashPos != string::npos)
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
                // simple face triangulation about f[0]
                // points may be incomplete
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    label fp2 = (fp1 + 1) % f.size();

                    dynFaces.append(triFace(f[0], f[fp1], f[fp2]));
                    dynRegions.append(regionI);
                    dynSizes[regionI]++;
                }
            }
            else
            {
                dynFaces.append(Face(f));
                dynRegions.append(regionI);
                dynSizes[regionI]++;
            }
        }
    }


    // transfer to normal lists
    this->storedPoints().transfer(dynPoints);

    sortFacesAndStore(dynFaces.xfer(), dynRegions.xfer(), sorted);

    // add regions, culling empty ones
    this->addRegions(dynSizes, dynNames, true);
    return true;
}


template<class Face>
void Foam::fileFormats::OBJsurfaceFormat<Face>::write
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();
    const List<surfRegion>& regionLst = surf.regions();

    writeHeader(os, surf.points(), faceLst.size(), regionLst);

    label faceIndex = 0;
    forAll(regionLst, regionI)
    {
        const surfRegion& reg = regionLst[regionI];

        os << "g " << reg.name() << endl;

        forAll(reg, localFaceI)
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


template<class Face>
void Foam::fileFormats::OBJsurfaceFormat<Face>::write
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();

    labelList faceMap;
    List<surfRegion> regionLst = surf.sortedRegions(faceMap);

    writeHeader(os, surf.points(), faceLst.size(), regionLst);

    label faceIndex = 0;
    forAll(regionLst, regionI)
    {
        // Print all faces belonging to this region
        const surfRegion& reg = regionLst[regionI];

        os << "g " << reg.name() << endl;

        forAll(reg, localFaceI)
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

// ************************************************************************* //
