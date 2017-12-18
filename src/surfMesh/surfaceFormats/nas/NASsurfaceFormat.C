/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "NASsurfaceFormat.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline Foam::label Foam::fileFormats::NASsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const Face& f,
    const label groupId,
    label elementId
)
{
    const label n = f.size();

    if (n == 3)
    {
        os  << "CTRIA3" << ','
            << ++elementId << ','
            << (groupId + 1) << ','
            << (f[0] + 1) << ','
            << (f[1] + 1) << ','
            << (f[2] + 1) << nl;
    }
    else if (n == 4)
    {
        os  << "CTRIA3" << ','
            << ++elementId << ','
            << (groupId + 1) << ','
            << (f[0] + 1) << ','
            << (f[1] + 1) << ','
            << (f[2] + 1) << ','
            << (f[3] + 1) << nl;
    }
    else
    {
        // simple triangulation about f[0].
        // better triangulation should have been done before
        for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
        {
            const label fp2 = f.fcIndex(fp1);

            os  << "CTRIA3" << ','
                << ++elementId << ','
                << (groupId + 1) << ','
                << (f[0] + 1) << ','
                << (f[fp1] + 1) << ','
                << (f[fp2] + 1) << nl;
        }
    }

    return elementId;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::NASsurfaceFormat<Face>::NASsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::NASsurfaceFormat<Face>::read
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

    // Nastran index of points
    DynamicList<label>  pointId;
    DynamicList<point>  dynPoints;
    DynamicList<Face>   dynFaces;
    DynamicList<label>  dynZones;
    DynamicList<label>  dynSizes;
    Map<label>          zoneLookup;

    // assume the types are not intermixed
    // leave faces that didn't have a group in 0
    bool sorted = true;
    label zoneI = 0;

    // Name for face group
    Map<word> nameLookup;

    // Ansa tags. Denoted by $ANSA_NAME.
    // These will appear just before the first use of a type.
    // We read them and store the PSHELL types which are used to name
    // the zones.
    label ansaId = -1;
    word  ansaType, ansaName;

    // A single warning per unrecognized command
    HashSet<word> unhandledCmd;

    while (is.good())
    {
        string::size_type linei = 0;  // parsing position within current line
        string line;
        is.getLine(line);

        // ANSA extension
        // line 1: $ANSA_NAME;<int>;<word>;
        // line 2: $partName
        if (line.startsWith("$ANSA_NAME"))
        {
            const auto sem0 = line.find(';', 0);
            const auto sem1 = line.find(';', sem0+1);
            const auto sem2 = line.find(';', sem1+1);

            if
            (
                sem0 != std::string::npos
             && sem1 != std::string::npos
             && sem2 != std::string::npos
            )
            {
                ansaId = readLabel(line.substr(sem0+1, sem1-sem0-1));
                ansaType = line.substr(sem1+1, sem2-sem1-1);

                string rawName;
                is.getLine(rawName);
                rawName.removeEnd("\r");  // Possible CR-NL
                ansaName = word::validate(rawName.substr(1));

                // Info<< "ANSA tag for NastranID:" << ansaId
                //     << " of type " << ansaType
                //     << " name " << ansaName << endl;
            }
        }


        // HYPERMESH extension
        // $HMNAME COMP                   1"partName"
        if (line.startsWith("$HMNAME COMP") && line.find('"') != string::npos)
        {
            label groupId = readLabel(line.substr(16, 16));

            // word::validate automatically removes quotes too
            const word groupName = word::validate(line.substr(32));

            nameLookup.insert(groupId, groupName);
            // Info<< "group " << groupId << " => " << groupName << endl;
        }

        if (line.empty() || line[0] == '$')
        {
            continue; // Skip empty or comment
        }

        // Check if character 72 is continuation
        if (line.size() > 72 && line[72] == '+')
        {
            line.resize(72);

            while (true)
            {
                string buf;
                is.getLine(buf);

                if (buf.size() > 72 && buf[72] == '+')
                {
                    line += buf.substr(8, 64);
                }
                else
                {
                    line += buf.substr(8);
                    break;
                }
            }
        }

        // First word (column 0-8)
        const word cmd(word::validate(nextNasField(line, linei, 8)));

        if (cmd == "CTRIA3")
        {
            (void) nextNasField(line, linei, 8); // 8-16
            label groupId = readLabel(nextNasField(line, linei, 8)); // 16-24
            label a = readLabel(nextNasField(line, linei, 8)); // 24-32
            label b = readLabel(nextNasField(line, linei, 8)); // 32-40
            label c = readLabel(nextNasField(line, linei, 8)); // 40-48

            // Convert groupId into zoneId
            const auto fnd = zoneLookup.cfind(groupId);
            if (fnd.found())
            {
                if (zoneI != fnd())
                {
                    // pshell types are intermixed
                    sorted = false;
                }
                zoneI = fnd();
            }
            else
            {
                zoneI = dynSizes.size();
                zoneLookup.insert(groupId, zoneI);
                dynSizes.append(0);
                // Info<< "zone" << zoneI << " => group " << groupId <<endl;
            }

            dynFaces.append(Face{a, b, c});
            dynZones.append(zoneI);
            dynSizes[zoneI]++;
        }
        else if (cmd == "CQUAD4")
        {
            (void) nextNasField(line, linei, 8); // 8-16
            label groupId = readLabel(nextNasField(line, linei, 8)); // 16-24
            label a = readLabel(nextNasField(line, linei, 8)); // 24-32
            label b = readLabel(nextNasField(line, linei, 8)); // 32-40
            label c = readLabel(nextNasField(line, linei, 8)); // 40-48
            label d = readLabel(nextNasField(line, linei, 8)); // 48-56

            // Convert groupID into zoneId
            const auto fnd = zoneLookup.cfind(groupId);
            if (fnd.found())
            {
                if (zoneI != fnd())
                {
                    // pshell types are intermixed
                    sorted = false;
                }
                zoneI = fnd();
            }
            else
            {
                zoneI = dynSizes.size();
                zoneLookup.insert(groupId, zoneI);
                dynSizes.append(0);
                // Info<< "zone" << zoneI << " => group " << groupId <<endl;
            }

            if (faceTraits<Face>::isTri())
            {
                dynFaces.append(Face{a, b, c});
                dynFaces.append(Face{c, d, a});
                dynZones.append(zoneI);
                dynZones.append(zoneI);
                dynSizes[zoneI] += 2;
            }
            else
            {
                dynFaces.append(Face{a,b,c,d});
                dynZones.append(zoneI);
                dynSizes[zoneI]++;
            }
        }
        else if (cmd == "GRID")
        {
            label index = readLabel(nextNasField(line, linei, 8)); // 8-16
            (void) nextNasField(line, linei, 8); // 16-24
            scalar x = readNasScalar(nextNasField(line, linei, 8)); // 24-32
            scalar y = readNasScalar(nextNasField(line, linei, 8)); // 32-40
            scalar z = readNasScalar(nextNasField(line, linei, 8)); // 40-48

            pointId.append(index);
            dynPoints.append(point(x, y, z));
        }
        else if (cmd == "GRID*")
        {
            // Long format is on two lines with '*' continuation symbol
            // on start of second line.
            // Typical line (spaces compacted)
            // GRID*      126   0 -5.55999875E+02 -5.68730474E+02
            // *         2.14897901E+02

            label index = readLabel(nextNasField(line, linei, 16)); // 8-24
            (void) nextNasField(line, linei, 16); // 24-40
            scalar x = readNasScalar(nextNasField(line, linei, 16)); // 40-56
            scalar y = readNasScalar(nextNasField(line, linei, 16)); // 56-72

            linei = 0; // restart at index 0
            is.getLine(line);
            if (line[0] != '*')
            {
                FatalErrorInFunction
                    << "Expected continuation symbol '*' when reading GRID*"
                    << " (double precision coordinate) format" << nl
                    << "Read:" << line << nl
                    << "File:" << is.name() << " line:" << is.lineNumber()
                    << exit(FatalError);
            }
            (void) nextNasField(line, linei, 8); // 0-8
            scalar z = readNasScalar(nextNasField(line, linei, 16)); // 8-16

            pointId.append(index);
            dynPoints.append(point(x, y, z));
        }
        else if (cmd == "PSHELL")
        {
            // Read shell type since group gives patchnames (ANSA extension)
            label groupId = readLabel(nextNasField(line, linei, 8)); // 8-16
            if (groupId == ansaId && ansaType == "PSHELL")
            {
                const word groupName = word::validate(ansaName);
                nameLookup.insert(groupId, groupName);
                // Info<< "group " << groupId << " => " << groupName << endl;
            }
        }
        else if (unhandledCmd.insert(cmd))
        {
            Info<< "Unhandled Nastran command " << line << nl
                << "File:" << is.name() << " line:" << is.lineNumber()
                << endl;
        }
    }

    //    Info<< "Read faces:" << dynFaces.size()
    //        << " points:" << dynPoints.size()
    //        << endl;

    // transfer to normal lists
    this->storedPoints().transfer(dynPoints);

    pointId.shrink();
    dynFaces.shrink();

    // Build inverse mapping (NASTRAN pointId -> index)
    Map<label> mapPointId(2*pointId.size());
    forAll(pointId, i)
    {
        mapPointId.insert(pointId[i], i);
    }

    // Relabel faces
    // ~~~~~~~~~~~~~
    for (Face& f : dynFaces)
    {
        forAll(f, fp)
        {
            f[fp] = mapPointId[f[fp]];
        }
    }
    pointId.clearStorage();
    mapPointId.clear();


    // create default zone names, or from ANSA/Hypermesh information
    List<word> names(dynSizes.size());
    forAllConstIters(zoneLookup, iter)
    {
        const label groupId = iter.key();
        const label zoneId  = iter.object();

        const auto fnd = nameLookup.cfind(groupId);
        if (fnd.found())
        {
            names[zoneId] = fnd();
        }
        else
        {
            names[zoneId] = word("zone") + ::Foam::name(zoneId);
        }
    }

    this->sortFacesAndStore(dynFaces.xfer(), dynZones.xfer(), sorted);
    this->addZones(dynSizes, names, true); // add zones, cull empty ones
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::NASsurfaceFormat<Face>::write
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

    // For simplicity, use fieldFormat::FREE throughout
    fileFormats::NASCore::setPrecision(os, fieldFormat::FREE);

    os  << "CEND" << nl
        << "TITLE = " << os.name().nameLessExt() << nl;

    // Print zone names as comment
    forAll(zones, zonei)
    {
        // HYPERMESH extension
        os  << "$HMNAME COMP" << setw(20) << (zonei+1)
            << '"' << zones[zonei].name() << '"' << nl;
    }

    // Write vertex coords with 1-based point Id
    os  << "$ GRID POINTS" << nl
        << "BEGIN BULK" << nl;

    label pointId = 0;
    for (const point& pt : pointLst)
    {
        os  << "GRID" << ','
            << ++pointId << ','
            << 0 << ','  // global coordinate system
            << pt.x() << ',' << pt.y() << ',' << pt.z() << nl;
    }

    os << "$ ELEMENTS" << nl;

    label faceIndex = 0;
    label zoneIndex = 0;
    label elementId = 0;
    for (const surfZone& zone : zones)
    {
        const label nLocalFaces = zone.size();

        if (useFaceMap)
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];
                elementId = writeShell(os, f, zoneIndex, elementId);
            }
        }
        else
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceIndex++];
                elementId = writeShell(os, f, zoneIndex, elementId);
            }
        }

        ++zoneIndex;
    }

    os << "ENDDATA" << nl;
}


// ************************************************************************* //
