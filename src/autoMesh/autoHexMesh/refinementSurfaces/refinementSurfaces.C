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

Class
    refinementSurfaces

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

    forAll(surfaceDicts, i)
    {
        const dictionary& dict = surfaceDicts[i];

        names_[i] = word(dict.lookup("name"));

        globalMinLevel[i] = readLabel(dict.lookup("minRefinementLevel"));
        globalMaxLevel[i] = readLabel(dict.lookup("maxRefinementLevel"));

        if
        (
            globalMinLevel[i] < 0
         || globalMaxLevel[i] < globalMinLevel[i]
        )
        {
            FatalErrorIn
            (
                "refinementSurfaces::refinementSurfaces"
                "(const IOobject&, const PtrList<dictionary>&)"
            )   << "Illegal level specification for surface " << names_[i]
                << " : minLevel:" << globalMinLevel[i]
                << " maxLevel:" << globalMaxLevel[i]
                << exit(FatalError);
        }

        if (dict.found("faceZone"))
        {
            dict.lookup("faceZone") >> faceZoneNames_[i];
            dict.lookup("cellZone") >> cellZoneNames_[i];
            dict.lookup("zoneInside") >> zoneInside_[i];
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

                regionMinLevel[i].insert(regionName, min);
                regionMaxLevel[i].insert(regionName, max);
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

    // From global region number to refinement level
    minLevel_.setSize(nRegions);
    minLevel_ = 0;
    maxLevel_.setSize(nRegions);
    maxLevel_ = 0;

    forAll(surfaceDicts, surfI)
    {
        const geometricSurfacePatchList& regions = operator[](surfI).patches();

        forAll(regions, i)
        {
            minLevel_[regionOffset_[surfI] + i] = globalMinLevel[surfI];
            maxLevel_[regionOffset_[surfI] + i] = globalMaxLevel[surfI];
        }

        forAllConstIter(HashTable<label>, regionMinLevel[surfI], iter)
        {
            // Find the patch
            forAll(regions, regionI)
            {
                if (regions[regionI].name() == iter.key())
                {
                    label globalRegionI = regionOffset_[surfI] + regionI;

                    minLevel_[globalRegionI] = iter();
                    maxLevel_[globalRegionI] =
                        regionMaxLevel[surfI][iter.key()];
                    break;
                }
            }
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
