/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "STARCDsurfaceFormat.H"
#include "ListOps.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::STARCDsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const Face& f,
    const label cellId,
    const label cellTableId
)
{
    os  << cellId                    // includes 1 offset
        << ' ' << starcdShell        // 3(shell) shape
        << ' ' << f.size()
        << ' ' << cellTableId
        << ' ' << starcdShellType;   // 4(shell)

    // Primitives have <= 8 vertices, but prevent overrun anyhow
    // indent following lines for ease of reading
    label count = 0;
    for (const label pointi : f)
    {
        if ((count % 8) == 0)
        {
            os  << nl << "  " << cellId;
        }
        os  << ' ' << pointi + 1;
        ++count;
    }
    os  << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::STARCDsurfaceFormat<Face>::STARCDsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STARCDsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    this->clear();

    fileName baseName = filename.lessExt();

    // Read cellTable names (if possible)
    Map<word> cellTableLookup = readInpCellTable
    (
        IFstream(starFileName(baseName, STARCDCore::INP_FILE))()
    );


    // STARCD index of points
    List<label> pointId;

    // read points from .vrt file
    readPoints
    (
        IFstream(starFileName(baseName, STARCDCore::VRT_FILE))(),
        this->storedPoints(),
        pointId
    );

    // Build inverse mapping (STARCD pointId -> index)
    Map<label> mapPointId(2*pointId.size());
    forAll(pointId, i)
    {
        mapPointId.insert(pointId[i], i);
    }
    pointId.clear();


    // Read .cel file
    // ~~~~~~~~~~~~~~
    IFstream is(starFileName(baseName, STARCDCore::CEL_FILE));
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    readHeader(is, STARCDCore::HEADER_CEL);

    DynamicList<Face>  dynFaces;
    DynamicList<label> dynZones;
    DynamicList<word>  dynNames;
    DynamicList<label> dynSizes;
    Map<label> lookup;

    // assume the cellTableIds are not intermixed
    bool sorted = true;
    label zoneId = 0;

    label ignoredLabel, shapeId, nLabels, cellTableId, typeId;
    DynamicList<label> vertexLabels(64);

    token tok;

    while (is.read(tok).good() && tok.isLabel())
    {
        // const label starCellId = tok.labelToken();
        is  >> shapeId
            >> nLabels
            >> cellTableId
            >> typeId;

        vertexLabels.clear();
        vertexLabels.reserve(nLabels);

        // Read indices - max 8 per line
        for (label i = 0; i < nLabels; ++i)
        {
            label vrtId;
            if ((i % 8) == 0)
            {
                is >> ignoredLabel; // Skip cellId for continuation lines
            }
            is >> vrtId;

            // Convert original vertex id to point label
            vertexLabels.append(mapPointId[vrtId]);
        }

        if (typeId == starcdShellType)
        {
            // Convert cellTableId to zoneId
            const auto iterGroup = lookup.cfind(cellTableId);
            if (iterGroup.found())
            {
                if (zoneId != *iterGroup)
                {
                    // cellTableIds are intermixed
                    sorted = false;
                }
                zoneId = *iterGroup;
            }
            else
            {
                zoneId = dynSizes.size();
                lookup.insert(cellTableId, zoneId);

                const auto iterTableName = cellTableLookup.cfind(cellTableId);

                if (iterTableName.found())
                {
                    dynNames.append(*iterTableName);
                }
                else
                {
                    dynNames.append("cellTable_" + ::Foam::name(cellTableId));
                }

                dynSizes.append(0);
            }

            SubList<label> vertices(vertexLabels, vertexLabels.size());
            if (faceTraits<Face>::isTri() && nLabels > 3)
            {
                // face needs triangulation
                face f(vertices);

                faceList trias(f.nTriangles());
                label nTri = 0;
                f.triangles(this->points(), nTri, trias);

                for (const face& tri : trias)
                {
                    // A triangular 'face', convert to 'triFace' etc
                    dynFaces.append(Face(tri));
                    dynZones.append(zoneId);
                    dynSizes[zoneId]++;
                }
            }
            else if (nLabels >= 3)
            {
                dynFaces.append(Face(vertices));
                dynZones.append(zoneId);
                dynSizes[zoneId]++;
            }
        }
    }
    mapPointId.clear();

    this->sortFacesAndStore(dynFaces, dynZones, sorted);

    // Add zones (retaining empty ones)
    this->addZones(dynSizes, dynNames);
    this->addZonesToFaces(); // for labelledTri

    return true;
}


template<class Face>
void Foam::fileFormats::STARCDsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    const dictionary&
)
{
    const UList<point>& pointLst = surf.points();
    const UList<Face>&  faceLst  = surf.surfFaces();
    const UList<label>& faceMap  = surf.faceMap();

    const surfZoneList zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst)
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    fileName baseName = filename.lessExt();

    writePoints
    (
        OFstream(starFileName(baseName, STARCDCore::VRT_FILE))(),
        pointLst
    );
    OFstream os(starFileName(baseName, STARCDCore::CEL_FILE));
    writeHeader(os, STARCDCore::HEADER_CEL);

    label faceIndex = 0;
    forAll(zones, zonei)
    {
        const surfZone& zone = zones[zonei];
        const label nLocalFaces = zone.size();

        if (useFaceMap)
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];
                writeShell(os, f, faceIndex, zonei + 1);
            }
        }
        else
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceIndex++];
                writeShell(os, f, faceIndex, zonei + 1);
            }
        }
    }

    // Write simple .inp file
    writeCase
    (
        OFstream(starFileName(baseName, STARCDCore::INP_FILE))(),
        pointLst,
        faceLst.size(),
        zones
    );
}


// ************************************************************************* //
