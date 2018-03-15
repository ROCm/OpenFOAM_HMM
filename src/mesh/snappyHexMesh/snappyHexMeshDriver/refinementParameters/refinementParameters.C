/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "refinementParameters.H"
#include "unitConversion.H"
#include "polyMesh.H"
#include "globalIndex.H"
#include "Tuple2.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementParameters::refinementParameters(const dictionary& dict)
:
    maxGlobalCells_(readLabel(dict.lookup("maxGlobalCells"))),
    maxLocalCells_(readLabel(dict.lookup("maxLocalCells"))),
    minRefineCells_(readLabel(dict.lookup("minRefinementCells"))),
    planarAngle_
    (
        dict.lookupOrDefault
        (
            "planarAngle",
            readScalar(dict.lookup("resolveFeatureAngle"))
        )
    ),
    nBufferLayers_(readLabel(dict.lookup("nCellsBetweenLevels"))),
    locationsOutsideMesh_
    (
        dict.lookupOrDefault
        (
            "locationsOutsideMesh",
            pointField(0)
        )
    ),
    faceZoneControls_(dict.subOrEmptyDict("faceZoneControls")),
    allowFreeStandingZoneFaces_(dict.lookup("allowFreeStandingZoneFaces")),
    useTopologicalSnapDetection_
    (
        dict.lookupOrDefault<bool>("useTopologicalSnapDetection", true)
    ),
    maxLoadUnbalance_(dict.lookupOrDefault<scalar>("maxLoadUnbalance", 0)),
    handleSnapProblems_
    (
        dict.lookupOrDefault<Switch>("handleSnapProblems", true)
    ),
    interfaceRefine_
    (
        dict.lookupOrDefault<Switch>("interfaceRefine", true)
    ),
    nErodeCellZone_(dict.lookupOrDefault<label>("nCellZoneErodeIter", 0))
{
    point locationInMesh;
    if (dict.readIfPresent("locationInMesh", locationInMesh))
    {
        locationsInMesh_.append(locationInMesh);
        zonesInMesh_.append("none");    // special name for no cellZone

        if (dict.found("locationsInMesh"))
        {
            FatalIOErrorInFunction(dict)
                << "Cannot both specify 'locationInMesh' and 'locationsInMesh'"
                << exit(FatalIOError);
        }
    }

    List<Tuple2<point, word>> pointsToZone;
    if (dict.readIfPresent("locationsInMesh", pointsToZone))
    {
        label nZones = locationsInMesh_.size();
        locationsInMesh_.setSize(nZones+pointsToZone.size());
        zonesInMesh_.setSize(locationsInMesh_.size());

        forAll(pointsToZone, i)
        {
            locationsInMesh_[nZones] = pointsToZone[i].first();
            zonesInMesh_[nZones] = pointsToZone[i].second();
            if (zonesInMesh_[nZones] == word::null)
            {
                zonesInMesh_[nZones] = "none";
            }
            nZones++;
        }
    }


    scalar featAngle(readScalar(dict.lookup("resolveFeatureAngle")));

    if (featAngle < 0 || featAngle > 180)
    {
        curvature_ = -GREAT;
    }
    else
    {
        curvature_ = Foam::cos(degToRad(featAngle));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary Foam::refinementParameters::getZoneInfo
(
    const word& fzName,
    surfaceZonesInfo::faceZoneType& faceType
) const
{
    dictionary patchInfo;
    patchInfo.add("type", wallPolyPatch::typeName);
    faceType = surfaceZonesInfo::INTERNAL;

    if (faceZoneControls_.found(fzName))
    {
        const dictionary& fzDict = faceZoneControls_.subDict(fzName);

        if (fzDict.found("patchInfo"))
        {
            patchInfo = fzDict.subDict("patchInfo");
        }

        word faceTypeName;
        if (fzDict.readIfPresent("faceType", faceTypeName))
        {
            faceType = surfaceZonesInfo::faceZoneTypeNames[faceTypeName];
        }
    }
    return patchInfo;
}


Foam::labelList Foam::refinementParameters::addCellZonesToMesh
(
    polyMesh& mesh
) const
{
    labelList zoneIDs(zonesInMesh_.size(), -1);
    forAll(zonesInMesh_, i)
    {
        if (zonesInMesh_[i] != word::null && zonesInMesh_[i] != "none")
        {
            zoneIDs[i] = surfaceZonesInfo::addCellZone
            (
                zonesInMesh_[i],    // name
                labelList(0),       // addressing
                mesh
            );
        }
    }
    return zoneIDs;
}


Foam::labelList Foam::refinementParameters::findCells
(
    const bool checkInsideMesh,
    const polyMesh& mesh,
    const pointField& locations
)
{
    // Force calculation of tet-diag decomposition (for use in findCell)
    (void)mesh.tetBasePtIs();

    // Global calculation engine
    globalIndex globalCells(mesh.nCells());

    // Cell label per point
    labelList cellLabels(locations.size());

    forAll(locations, i)
    {
        const point& location = locations[i];

        label localCellI = mesh.findCell(location);

        label globalCellI = -1;

        if (localCellI != -1)
        {
            globalCellI = globalCells.toGlobal(localCellI);
        }

        reduce(globalCellI, maxOp<label>());

        if (checkInsideMesh && globalCellI == -1)
        {
            FatalErrorInFunction
                << "Point " << location
                << " is not inside the mesh or on a face or edge." << nl
                << "Bounding box of the mesh:" << mesh.bounds()
                << exit(FatalError);
        }


        label procI = globalCells.whichProcID(globalCellI);
        label procCellI = globalCells.toLocal(procI, globalCellI);

        Info<< "Found point " << location << " in cell " << procCellI
            << " on processor " << procI << endl;

        if (globalCells.isLocal(globalCellI))
        {
            cellLabels[i] = localCellI;
        }
        else
        {
            cellLabels[i] = -1;
        }
    }
    return cellLabels;
}


Foam::labelList Foam::refinementParameters::zonedLocations
(
    const wordList& zonesInMesh
)
{
    DynamicList<label> indices(zonesInMesh.size());

    forAll(zonesInMesh, i)
    {
        if
        (
            zonesInMesh[i] != word::null
         && zonesInMesh[i] != "none"
        )
        {
            indices.append(i);
        }
    }
    return indices;
}


Foam::labelList Foam::refinementParameters::unzonedLocations
(
    const wordList& zonesInMesh
)
{
    DynamicList<label> indices(0);

    forAll(zonesInMesh, i)
    {
        if
        (
            zonesInMesh[i] == word::null
         || zonesInMesh[i] == "none"
        )
        {
            indices.append(i);
        }
    }
    return indices;
}


// ************************************************************************* //
