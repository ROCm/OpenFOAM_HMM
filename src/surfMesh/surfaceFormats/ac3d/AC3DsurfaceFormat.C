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

#include "AC3DsurfaceFormat.H"
#include "StringStream.H"
#include "PrimitivePatch.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::AC3DsurfaceFormat<Face>::AC3DsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::AC3DsurfaceFormat<Face>::read
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

    string line, cmd, args;

    is.getLine(line);

    // Verify version
    {
        const string version = line.substr(4);

        if (version != "b")
        {
            WarningInFunction
                << "When reading AC3D file " << filename
                << " read header " << line << " with version "
                << version << endl
                << "Only tested reading with version 'b'."
                << " This might give problems" << endl;
        }
    }


    if (!cueTo(is, "OBJECT", args) || args != "world")
    {
        FatalErrorInFunction
            << "Cannot find \"OBJECT world\" in file " << filename
            << exit(FatalError);
    }

    // Number of kids is the number of zones
    args = cueToOrDie(is, "kids");
    const label nZones = parse<int>(args);

    // Start of vertices for object/zones
    label vertexOffset = 0;

    DynamicList<point> dynPoints;
    DynamicList<Face>  dynFaces;
    List<word>         names(nZones);
    List<label>        sizes(nZones, 0);

    for (label zoneI = 0; zoneI < nZones; ++zoneI)
    {
        names[zoneI] = word("zone") + Foam::name(zoneI);

        args = cueToOrDie(is, "OBJECT", "while reading " + names[zoneI]);

        // number of vertices for this zone
        label  nZonePoints = 0;
        vector location(Zero);
        // tensor rotation(I);

        // Read all info for current zone
        while (is.good())
        {
            // Read line and get first word. If end of file break since
            // zone should always end with 'kids' command ?not sure.
            if (!readCmd(is, cmd, args))
            {
                FatalErrorInFunction
                    << "Did not read up to \"kids 0\" while reading zone "
                    << zoneI << " from file " << filename
                    << exit(FatalError);
            }

            if (cmd == "name")
            {
                // name %s
                const string str = parse<string>(args);
                names[zoneI] = word::validate(str);
            }
            else if (cmd == "rot")
            {
                // rot  %f %f %f  %f %f %f  %f %f %f

                // IStringStream lineStream(args);
                //
                // lineStream
                //     >> rotation.xx() >> rotation.xy() >> rotation.xz()
                //     >> rotation.yx() >> rotation.yy() >> rotation.yz()
                //     >> rotation.zx() >> rotation.zy() >> rotation.zz();

                WarningInFunction
                    << "rot (rotation tensor) command not implemented"
                    << "Line:" << cmd << ' ' << args << endl
                    << "while reading zone " << zoneI << endl;
            }
            else if (cmd == "loc")
            {
                // loc  %f %f %f
                IStringStream lineStream(args);

                lineStream
                    >> location.x()
                    >> location.y()
                    >> location.z();
            }
            else if (cmd == "numvert")
            {
                // numvert  %d
                nZonePoints = parse<int>(args);

                for (label vertI = 0; vertI < nZonePoints; ++vertI)
                {
                    is.getLine(line);
                    IStringStream lineStream(line);

                    point pt;
                    lineStream
                        >> pt.x() >> pt.y() >> pt.z();

                    // Offset with current translation vector
                    dynPoints.append(location + pt);
                }
            }
            else if (cmd == "numsurf")
            {
                const label nFaces = parse<int>(args);

                for (label facei = 0; facei < nFaces; ++facei)
                {
                    const string errorMsg =
                        string(" while reading face ")
                            + Foam::name(facei) + " on zone "
                            + Foam::name(zoneI)
                            + " from file " + filename;

                    cueToOrDie(is, "SURF", errorMsg);
                    cueToOrDie(is, "mat", errorMsg);
                    args = cueToOrDie(is, "refs", errorMsg);

                    const label nVert = parse<int>(args);

                    List<label> verts(nVert);
                    forAll(verts, vertI)
                    {
                        is.getLine(line);
                        verts[vertI] = vertexOffset + parse<int>(line);
                    }

                    const labelUList& f = static_cast<const labelUList&>(verts);

                    if (faceTraits<Face>::isTri() && f.size() > 3)
                    {
                        // simple face triangulation about f[0]
                        // points may be incomplete
                        for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
                        {
                            label fp2 = f.fcIndex(fp1);

                            dynFaces.append(Face{f[0], f[fp1], f[fp2]});
                            sizes[zoneI]++;
                        }
                    }
                    else
                    {
                        dynFaces.append(Face(f));
                        sizes[zoneI]++;
                    }
                }

                // Done the current zone.
                // Increment the offset vertices are stored at
                vertexOffset += nZonePoints;
            }
            else if (cmd == "kids")
            {
                // 'kids' denotes the end of the current zone.
                const label nKids = parse<int>(args);

                if (nKids != 0)
                {
                    FatalErrorInFunction
                        << "Can only read objects without kids."
                        << " Encountered " << nKids << " kids when"
                        << " reading zone " << zoneI
                        << exit(FatalError);
                }

                // Done reading current zone
                break;
            }
        }
    }

    // transfer to normal lists
    this->storedPoints().transfer(dynPoints);
    this->storedFaces().transfer(dynFaces);

    // add zones, culling empty ones
    this->addZones(sizes, names, true);
    this->addZonesToFaces(); // for labelledTri
    this->stitchFaces(SMALL);

    return true;
}


namespace Foam
{
// file-scope writing of a patch of faces
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
static void writeZone
(
    Ostream& os,
    const PrimitivePatch<Face, FaceList, PointField, PointType>& patch,
    const word& name,
    const label zoneI
)
{
    // An isolated surface region (patch).
    os  << "OBJECT poly" << nl
        << "name \"" << name << "\"" << nl;

    os << "numvert " << patch.nPoints() << nl;

    for (const point& pt : patch.localPoints())
    {
        os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }

    os << "numsurf " << patch.size() << nl;

    for (const Face& f : patch.localFaces())
    {
        os  << "SURF 0x20" << nl          // polygon
            << "mat " << zoneI << nl
            << "refs " << f.size() << nl;

        for (const label verti : f)
        {
            os << verti << " 0 0" << nl;
        }
    }

    os << "kids 0" << endl;
}
}

template<class Face>
void Foam::fileFormats::AC3DsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    const dictionary&
)
{
    const pointField& pointLst = surf.points();
    const UList<Face>& faceLst = surf.surfFaces();

    const UList<surfZone>& zones =
    (
        surf.surfZones().size()
      ? surf.surfZones()
      : surfaceFormatsCore::oneZone(faceLst)
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    writeHeader(os, zones);

    if (zones.size() == 1)
    {
        PrimitivePatch<Face, UList, const pointField&> patch
        (
            faceLst, pointLst
        );

        writeZone(os, patch, zones[0].name(), 0);
        return;
    }

    label zoneIndex = 0;
    for (const surfZone& zone : zones)
    {
        if (useFaceMap)
        {
            SubList<label> zoneMap(surf.faceMap(), zone.size(), zone.start());
            PrimitivePatch<Face, UIndirectList, const pointField&> patch
            (
                UIndirectList<Face>(faceLst, zoneMap),
                pointLst
            );

            writeZone(os, patch, zone.name(), zoneIndex);
        }
        else
        {
            PrimitivePatch<Face, UList, const pointField&> patch
            (
                SubList<Face>(faceLst, zone.size(), zone.start()),
                pointLst
            );

            writeZone(os, patch, zone.name(), zoneIndex);
        }

        ++zoneIndex;
    }
}


template<class Face>
void Foam::fileFormats::AC3DsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf,
    const dictionary&
)
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    labelList faceMap;
    List<surfZone> zoneLst = surf.sortedZones(faceMap);

    if (zoneLst.size() <= 1)
    {
        const List<surfZone>& zones =
        (
            zoneLst.size()
          ? zoneLst
          : surfaceFormatsCore::oneZone(surf.surfFaces())
        );

        writeHeader(os, zones);
        writeZone(os, surf, zones[0].name(), 0);
        return;
    }

    writeHeader(os, zoneLst);

    label zoneIndex = 0;
    for (const surfZone& zone : zoneLst)
    {
        SubList<label> zoneMap(faceMap, zone.size(), zone.start());
        PrimitivePatch<Face, UIndirectList, const pointField&> patch
        (
            UIndirectList<Face>(surf.surfFaces(), zoneMap),
            surf.points()
        );

        writeZone(os, patch, zone.name(), zoneIndex);

        ++zoneIndex;
    }
}


// ************************************************************************* //
