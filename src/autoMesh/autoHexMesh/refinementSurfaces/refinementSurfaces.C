/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "refinementSurfaces.H"
#include "orientedSurface.H"
#include "Time.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileNameList Foam::refinementSurfaces::extractFileNames
(
    const PtrList<dictionary>& surfaceDicts
)
{
    fileNameList surfaceFileNames(surfaceDicts.size());

    forAll(surfaceDicts, i)
    {
        surfaceFileNames[i] = fileName(surfaceDicts[i].lookup("file"));
    }
    return surfaceFileNames;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::refinementSurfaces::refinementSurfaces
(
    const IOobject& io,
    const PtrList<dictionary>& surfaceDicts
)
:
    triSurfaceMeshes(io, extractFileNames(surfaceDicts)),
    names_(surfaceDicts.size()),
    closed_(surfaceDicts.size(), true),
    faceZoneNames_(surfaceDicts.size()),
    cellZoneNames_(surfaceDicts.size()),
    zoneInside_(surfaceDicts.size()),
    regionOffset_(surfaceDicts.size())
{
    labelList globalMinLevel(surfaceDicts.size(), 0);
    labelList globalMaxLevel(surfaceDicts.size(), 0);
    List<HashTable<label> > regionMinLevel(surfaceDicts.size());
    List<HashTable<label> > regionMaxLevel(surfaceDicts.size());

    labelList globalSurfLayers(surfaceDicts.size());
    List<HashTable<label> > regionSurfLayers(surfaceDicts.size());

    wordList globalPatchType(surfaceDicts.size());
    List<HashTable<word> > regionPatchType(surfaceDicts.size());
    List<HashTable<word> > regionPatchName(surfaceDicts.size());

    forAll(surfaceDicts, surfI)
    {
        const dictionary& dict = surfaceDicts[surfI];

        names_[surfI] = word(dict.lookup("name"));

        // Global refinement level
        globalMinLevel[surfI] = readLabel(dict.lookup("minRefinementLevel"));
        globalMaxLevel[surfI] = readLabel(dict.lookup("maxRefinementLevel"));

        // Global number of layers
        globalSurfLayers[surfI] = readLabel(dict.lookup("surfaceLayers"));

        // Global zone names per surface
        if (dict.found("faceZone"))
        {
            dict.lookup("faceZone") >> faceZoneNames_[surfI];
            dict.lookup("cellZone") >> cellZoneNames_[surfI];
            dict.lookup("zoneInside") >> zoneInside_[surfI];
        }

        // Global patch name per surface
        if (dict.found("patchType"))
        {
            dict.lookup("patchType") >> globalPatchType[surfI];
        }


        if (dict.found("regions"))
        {
            PtrList<dictionary> regionDicts(dict.lookup("regions"));

            forAll(regionDicts, dictI)
            {
                const dictionary& regionDict = regionDicts[dictI];

                const word regionName(regionDict.lookup("name"));
                label min = readLabel(regionDict.lookup("minRefinementLevel"));
                label max = readLabel(regionDict.lookup("maxRefinementLevel"));

                regionMinLevel[surfI].insert(regionName, min);
                regionMaxLevel[surfI].insert(regionName, max);

                label nLayers = readLabel(regionDict.lookup("surfaceLayers"));
                regionSurfLayers[surfI].insert(regionName, nLayers);

                if (regionDict.found("patchType"))
                {
                    regionPatchType[surfI].insert
                    (
                        regionName,
                        regionDict.lookup("patchType")
                    );
                    regionPatchName[surfI].insert
                    (
                        regionName,
                        regionDict.lookup("patchName")
                    );
                }
            }
        }
    }


    // Check for duplicate surface or region names
    {
        HashTable<label> surfaceNames(names_.size());

        forAll(names_, surfI)
        {
            if (!surfaceNames.insert(names_[surfI], surfI))
            {
                FatalErrorIn
                (
                    "refinementSurfaces::refinementSurfaces"
                    "(const IOobject&, const PtrList<dictionary>&)"
                )   << "Duplicate surface name " << names_[surfI] << endl
                    << "Previous occurrence of name at surface "
                    << surfaceNames[names_[surfI]]
                    << exit(FatalError);
            }

            // Check for duplicate region names
            const geometricSurfacePatchList& patches =
                operator[](surfI).patches();

            HashTable<label> regionNames(patches.size());
            forAll(patches, i)
            {
                if (!regionNames.insert(patches[i].name(), i))
                {
                    FatalErrorIn
                    (
                        "refinementSurfaces::refinementSurfaces"
                        "(const IOobject&, const PtrList<dictionary>&)"
                    )   << "Duplicate region name " << patches[i].name()
                        << " on surface " << names_[surfI] << endl
                        << "Previous occurrence of region at index "
                        << regionNames[patches[i].name()]
                        << exit(FatalError);
                }
            }
        }
    }


    // Calculate closedness
    forAll(closed_, surfI)
    {
        const triSurface& s = operator[](surfI);

        const labelListList& edgeFaces = s.edgeFaces();

        forAll(edgeFaces, edgeI)
        {
            if (edgeFaces[edgeI].size() != 2)
            {
                closed_[surfI] = false;
                break;
            }
        }
    }


    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaceDicts, surfI)
    {
        regionOffset_[surfI] = nRegions;

        nRegions += operator[](surfI).patches().size();
    }

    // Rework surface specific information into information per global region
    minLevel_.setSize(nRegions);
    minLevel_ = 0;
    maxLevel_.setSize(nRegions);
    maxLevel_ = 0;
    numLayers_.setSize(nRegions);
    numLayers_ = 0;
    patchName_.setSize(nRegions);
    patchType_.setSize(nRegions);

    forAll(surfaceDicts, surfI)
    {
        const geometricSurfacePatchList& regions = operator[](surfI).patches();

        // Initialise to global (i.e. per surface)
        forAll(regions, i)
        {
            minLevel_[regionOffset_[surfI] + i] = globalMinLevel[surfI];
            maxLevel_[regionOffset_[surfI] + i] = globalMaxLevel[surfI];
            numLayers_[regionOffset_[surfI] + i] = globalSurfLayers[surfI];
            patchType_[regionOffset_[surfI] + i] = globalPatchType[surfI];
        }

        // Get the region names
        wordList regionNames(regions.size());
        forAll(regions, regionI)
        {
            regionNames[regionI] = regions[regionI].name();
        }

        // Overwrite with region specific information
        forAllConstIter(HashTable<label>, regionMinLevel[surfI], iter)
        {
            // Find the local region number.
            label regionI = findIndex(regionNames, iter.key());

            if (regionI == -1)
            {
                FatalErrorIn
                (
                    "refinementSurfaces::refinementSurfaces"
                    "(const IOobject&, const PtrList<dictionary>&)"
                )   << "Cannot find region " << iter.key()
                    << " in surface " << names_[surfI]
                    << " which has regions " << regionNames
                    << abort(FatalError);
            }

            label globalRegionI = regionOffset_[surfI] + regionI;

            minLevel_[globalRegionI] = iter();
            maxLevel_[globalRegionI] = regionMaxLevel[surfI][iter.key()];
            numLayers_[globalRegionI] = regionSurfLayers[surfI][iter.key()];

            // Check validity
            if
            (
                minLevel_[globalRegionI] < 0
             || maxLevel_[globalRegionI] < minLevel_[globalRegionI]
             || numLayers_[globalRegionI] < 0
            )
            {
                FatalErrorIn
                (
                    "refinementSurfaces::refinementSurfaces"
                    "(const IOobject&, const PtrList<dictionary>&)"
                )   << "Illegal level or layer specification for surface "
                    << names_[surfI]
                    << " : minLevel:" << minLevel_[globalRegionI]
                    << " maxLevel:" << maxLevel_[globalRegionI]
                    << " numLayers:" << numLayers_[globalRegionI]
                    << exit(FatalError);
            }
        }

        // Optional patch names and patch types
        forAllConstIter(HashTable<word>, regionPatchName[surfI], iter)
        {
            label regionI = findIndex(regionNames, iter.key());
            label globalRegionI = regionOffset_[surfI] + regionI;

            patchName_[globalRegionI] = iter();
            patchType_[globalRegionI] = regionPatchType[surfI][iter.key()];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Get indices of named surfaces (surfaces with cellZoneNam)
Foam::labelList Foam::refinementSurfaces::getNamedSurfaces() const
{
   labelList namedSurfaces(cellZoneNames_.size());

    label namedI = 0;
    forAll(cellZoneNames_, surfI)
    {
        if (cellZoneNames_[surfI].size() > 0)
        {
            namedSurfaces[namedI++] = surfI;
        }
    }
    namedSurfaces.setSize(namedI);

    return namedSurfaces;
}


bool Foam::refinementSurfaces::isSurfaceClosed(const triSurface& s)
{
    if (s.nEdges() == s.nInternalEdges())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Orient surface (if they're closed) before any searching is done.
void Foam::refinementSurfaces::orientSurface
(
    const point& keepPoint,
    triSurfaceMesh& s
)
{
    // Flip surface so keepPoint is outside.
    bool anyFlipped = orientedSurface::orient(s, keepPoint, true);

    if (anyFlipped)
    {
        // orientedSurface will have done a clearOut of the surface.
        // we could do a clearout of the triSurfaceMeshes::trees()
        // but these aren't affected by orientation (except for cached
        // sideness which should not be set at this point. !!Should
        // check!)

        Info<< "orientSurfaces : Flipped orientation of surface "
            << s.searchableSurface::name()
            << " so point " << keepPoint << " is outside." << endl;
    }
}


// Count number of triangles per surface region
Foam::labelList Foam::refinementSurfaces::countRegions(const triSurface& s)
{
    const geometricSurfacePatchList& regions = s.patches();

    labelList nTris(regions.size(), 0);

    forAll(s, triI)
    {
        nTris[s[triI].region()]++;
    }
    return nTris;
}


// Find intersections of edge. Return -1 or first surface with higher minLevel
// number.
Foam::label Foam::refinementSurfaces::findHigherIntersection
(
    const point& start,
    const point& end,
    const label currentLevel,   // current cell refinement level

    pointIndexHit& hit
) const
{
    forAll(*this, surfI)
    {
        hit = operator[](surfI).findLineAny(start, end);

        if (hit.hit())
        {
            if (currentLevel == -1)
            {
                // Return any.
                return surfI;
            }
            else
            {
                //const triSurface& s = trees()[surfI].shapes().surface();
                //label region = globalRegion(surfI, s[hit.index()].region());
                //if (minLevel_[region] > currentLevel)
                if (minLevelFields_[surfI][hit.index()] > currentLevel)
                {
                    return surfI;
                }
            }
        }
    }

    hit.setMiss();

    return -1;
}


// Find intersection with max of edge. Return -1 or the surface
// with the highest maxLevel above currentLevel
Foam::label Foam::refinementSurfaces::findHighestIntersection
(
    const point& start,
    const point& end,
    const label currentLevel,   // current cell refinement level

    pointIndexHit& maxHit
) const
{
    // surface with highest maxlevel
    label maxSurface = -1;
    // maxLevel of maxSurface
    label maxLevel = currentLevel;

    forAll(*this, surfI)
    {
        pointIndexHit hit = operator[](surfI).findLineAny(start, end);

        if (hit.hit())
        {
            const triSurface& s = operator[](surfI);

            label region = globalRegion(surfI, s[hit.index()].region());

            if (maxLevel_[region] > maxLevel)
            {
                maxSurface = surfI;
                maxLevel = maxLevel_[region];
                maxHit = hit;
            }
        }
    }

    if (maxSurface == -1)
    {
        // maxLevel unchanged. No interesting surface hit.
        maxHit.setMiss();
    }

    return maxSurface;
}


Foam::label Foam::refinementSurfaces::insideZone
(
    const labelList& surfaces,
    const point& pt
) const
{
    forAll(surfaces, i)
    {
        label surfI = surfaces[i];

        if (closed_[surfI])
        {
            const triSurfaceMesh& s = operator[](surfI);

            triSurfaceMesh::volumeType t = s.getVolumeType(pt);

            if
            (
                (t == triSurfaceMesh::INSIDE && zoneInside_[surfI])
             || (t == triSurfaceMesh::OUTSIDE && !zoneInside_[surfI])
            )
            {
                return surfI;
            }
        }
    }
    return -1;
}


Foam::label Foam::refinementSurfaces::markInsidePoints
(
    const pointField& points,
    PackedList<1>& isInside
) const
{
    isInside.setSize(points.size());
    isInside = 0u;

    label nPointInside = 0;

    forAll(points, pointI)
    {
        forAll(*this, surfI)
        {
            if (closed()[surfI])
            {
                const triSurfaceMesh& s = operator[](surfI);

                triSurfaceMesh::volumeType t = s.getVolumeType(points[pointI]);

                if (t == triSurfaceMesh::INSIDE)
                {
                    isInside.set(pointI, 1u);
                    nPointInside++;
                    break;
                }
            }
        }
    }

    return nPointInside;
}


void Foam::refinementSurfaces::setMinLevelFields
(
    const PtrList<searchableSurface>& shells,
    const labelList& shellLevels,
    const boolList& shellRefineInside
)
{
    typedef Foam::DimensionedField<label, triSurfaceGeoMesh>
        triSurfaceLabelField;

    minLevelFields_.setSize(size());

    forAll(*this, surfI)
    {
        const triSurfaceMesh& surfMesh = operator[](surfI);

        minLevelFields_.set
        (
            surfI,
            new triSurfaceLabelField
            (
                IOobject
                (
                    surfMesh.IOobject::name(),
                    surfMesh.time().constant(), // directory
                    "triSurface",               // instance
                    surfMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                surfMesh,
                dimless
            )
        );
    }

    // Initialise fields to region wise minLevel
    forAll(*this, surfI)
    {
        const triSurface& s = operator[](surfI);

        triSurfaceLabelField& minLevelField = minLevelFields_[surfI];

        forAll(s, i)
        {
            minLevelField[i] = minLevel_[globalRegion(surfI, s[i].region())];
        }
    }

    // Adapt for triangles inside shells
    forAll(*this, surfI)
    {
        const triSurface& s = operator[](surfI);
        triSurfaceLabelField& minLevelField = minLevelFields_[surfI];

        forAll(s, i)
        {
            point fc = s[i].centre(s.points());

            forAll(shells, shellI)
            {
                // Which side of shell is to be refined
                searchableSurface::volumeType refSide =
                (
                    shellRefineInside[shellI]
                  ? searchableSurface::INSIDE
                  : searchableSurface::OUTSIDE
                );

                // Find whether point is inside or outside shell
                searchableSurface::volumeType t =
                    shells[shellI].getVolumeType(fc);

                if (t == refSide)
                {
                    minLevelField[i] = max
                    (
                        minLevelField[i],
                        shellLevels[shellI]
                    );
                }
            }
        }
    }
}


// ************************************************************************* //
