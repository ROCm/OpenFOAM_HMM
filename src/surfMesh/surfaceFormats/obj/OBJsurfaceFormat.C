/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "OBJsurfaceFormat.H"
#include "clock.H"
#include "Fstream.H"
#include "StringStream.H"
#include "ListOps.H"
#include "faceTraits.H"

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
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // assume that the groups are not intermixed
    bool sorted = true;

    DynamicList<point> dynPoints;
    DynamicList<Face>  dynFaces;
    DynamicList<label> dynZones;
    DynamicList<word>  dynNames;
    DynamicList<label> dynSizes;
    HashTable<label>   lookup;

    // place faces without a group in zone0
    label zoneI = 0;
    lookup.insert("zone0", zoneI);
    dynNames.append("zone0");
    dynSizes.append(0);

    while (is.good())
    {
        string line = this->getLineNoComment(is);

        // Handle continuations
        if (line.removeEnd("\\"))
        {
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
                if (zoneI != fnd())
                {
                    // group appeared out of order
                    sorted = false;
                }
                zoneI = fnd();
            }
            else
            {
                zoneI = dynSizes.size();
                lookup.insert(name, zoneI);
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
                    vertexSpec = line.substr(startNum);
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

            labelUList& f = static_cast<labelUList&>(dynVertices);

            if (faceTraits<Face>::isTri() && f.size() > 3)
            {
                // simple face triangulation about f[0]
                // points may be incomplete
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    const label fp2 = f.fcIndex(fp1);

                    dynFaces.append(Face{f[0], f[fp1], f[fp2]});
                    dynZones.append(zoneI);
                    dynSizes[zoneI]++;
                }
            }
            else
            {
                dynFaces.append(Face(f));
                dynZones.append(zoneI);
                dynSizes[zoneI]++;
            }
        }
    }


    // transfer to normal lists
    this->storedPoints().transfer(dynPoints);

    this->sortFacesAndStore(dynFaces.xfer(), dynZones.xfer(), sorted);
    this->addZones(dynSizes, dynNames, true); // add zones, cull empty ones
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::OBJsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    const dictionary& options
)
{
    const UList<point>& pointLst = surf.points();
    const UList<Face>&  faceLst  = surf.surfFaces();
    const UList<label>& faceMap  = surf.faceMap();

    // for no zones, suppress the group name
    const UList<surfZone>& zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst, "")
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    os  << "# Wavefront OBJ file written " << clock::dateTime().c_str() << nl
        << "o " << os.name().nameLessExt() << nl
        << nl
        << "# points : " << pointLst.size() << nl
        << "# faces  : " << faceLst.size() << nl
        << "# zones  : " << zones.size() << nl;

    // Print zone names as comment
    forAll(zones, zonei)
    {
        os  << "#   " << zonei << "  " << zones[zonei].name()
            << "  (nFaces: " << zones[zonei].size() << ")" << nl;
    }

    os  << nl
        << "# <points count=\"" << pointLst.size() << "\">" << nl;

    // Write vertex coords
    for (const point& pt : pointLst)
    {
        os  << "v " << pt.x() << ' '  << pt.y() << ' '  << pt.z() << nl;
    }

    os  << "# </points>" << nl
        << nl
        << "# <faces count=\"" << faceLst.size() << "\">" << nl;


    label faceIndex = 0;
    for (const surfZone& zone : zones)
    {
        if (zone.name().size())
        {
            os << "g " << zone.name() << nl;
        }

        const label nLocalFaces = zone.size();

        if (useFaceMap)
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                os << 'f';
                forAll(f, fp)
                {
                    os << ' ' << f[fp] + 1;
                }
                os << nl;
            }
        }
        else
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceIndex++];

                os << 'f';
                forAll(f, fp)
                {
                    os << ' ' << f[fp] + 1;
                }
                os << nl;
            }
        }
    }
    os << "# </faces>" << nl;
}


// ************************************************************************* //
