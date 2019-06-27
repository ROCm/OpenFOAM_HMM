/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "meshRefinement.H"
#include "fvMesh.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "removeCells.H"
#include "unitConversion.H"
#include "bitSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::label Foam::meshRefinement::markFakeGapRefinement
//(
//    const scalar planarCos,
//
//    const label nAllowRefine,
//    const labelList& neiLevel,
//    const pointField& neiCc,
//
//    labelList& refineCell,
//    label& nRefine
//) const
//{
//    label oldNRefine = nRefine;
//
//
//    // Collect candidate faces (i.e. intersecting any surface and
//    // owner/neighbour not yet refined.
//    const labelList testFaces(getRefineCandidateFaces(refineCell));
//
//    // Collect segments
//    pointField start(testFaces.size());
//    pointField end(testFaces.size());
//    labelList minLevel(testFaces.size());
//
//    calcCellCellRays
//    (
//        neiCc,
//        neiLevel,
//        testFaces,
//        start,
//        end,
//        minLevel
//    );
//
//
//    // Re-use the gap shooting methods. This needs:
//    //  - shell gapLevel : faked. Set to 0,labelMax
//    //  - surface gapLevel : faked by overwriting
//
//
//    List<FixedList<label, 3>>& surfGapLevel = const_cast
//    <
//        List<FixedList<label, 3>>&
//    >(surfaces_.extendedGapLevel());
//
//    List<volumeType>& surfGapMode = const_cast
//    <
//        List<volumeType>&
//    >(surfaces_.extendedGapMode());
//
//    const List<FixedList<label, 3>> surfOldLevel(surfGapLevel);
//    const List<volumeType> surfOldMode(surfGapMode);
//
//    // Set the extended gap levels
//    forAll(surfaces_.gapLevel(), regioni)
//    {
//        surfGapLevel[regioni] = FixedList<label, 3>
//        ({
//            3,
//           -1,
//            surfaces_.gapLevel()[regioni]+1
//        });
//    }
//    surfGapMode = volumeType::MIXED;
//
//Pout<< "gapLevel was:" << surfOldLevel << endl;
//Pout<< "gapLevel now:" << surfGapLevel << endl;
//Pout<< "gapMode was:" << surfOldMode << endl;
//Pout<< "gapMode now:" << surfGapMode << endl;
//Pout<< "nRefine was:" << oldNRefine << endl;
//
//
//
//    List<List<FixedList<label, 3>>>& shellGapLevel = const_cast
//    <
//        List<List<FixedList<label, 3>>>&
//    >(shells_.extendedGapLevel());
//
//    List<List<volumeType>>& shellGapMode = const_cast
//    <
//        List<List<volumeType>>&
//    >(shells_.extendedGapMode());
//
//    const List<List<FixedList<label, 3>>> shellOldLevel(shellGapLevel);
//    const List<List<volumeType>> shellOldMode(shellGapMode);
//
//    // Set the extended gap levels
//    forAll(shellGapLevel, shelli)
//    {
//        shellGapLevel[shelli] =  FixedList<label, 3>({3, -1, labelMax});
//        shellGapMode[shelli] = volumeType::MIXED;
//    }
//Pout<< "shellLevel was:" << shellOldLevel << endl;
//Pout<< "shellLevel now:" << shellGapLevel << endl;
//
//    const label nAdditionalRefined = markSurfaceGapRefinement
//    (
//        planarCos,
//
//        nAllowRefine,
//        neiLevel,
//        neiCc,
//
//        refineCell,
//        nRefine
//    );
//
//Pout<< "nRefine now:" << nRefine << endl;
//
//    // Restore
//    surfGapLevel = surfOldLevel;
//    surfGapMode = surfOldMode;
//    shellGapLevel = shellOldLevel;
//    shellGapMode = shellOldMode;
//
//    return nAdditionalRefined;
//}


void Foam::meshRefinement::markOutsideFaces
(
    const labelList& cellLevel,
    const labelList& neiLevel,
    const labelList& refineCell,
    bitSet& isOutsideFace
) const
{
    // Get faces:
    // - on outside of cell set
    // - inbetween same cell level (i.e. quads)

    isOutsideFace.setSize(mesh_.nFaces());
    isOutsideFace = Zero;

    forAll(mesh_.faceNeighbour(), facei)
    {
        label own = mesh_.faceOwner()[facei];
        label nei = mesh_.faceNeighbour()[facei];
        if
        (
            (cellLevel[own] == cellLevel[nei])
         && (
                (refineCell[own] != -1)
             != (refineCell[nei] != -1)
            )
        )
        {
            isOutsideFace.set(facei);
        }
    }
    {

        const label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();

        labelList neiRefineCell(nBnd);
        syncTools::swapBoundaryCellList(mesh_, refineCell, neiRefineCell);
        for (label bFacei = 0; bFacei < nBnd; ++bFacei)
        {
            label facei = mesh_.nInternalFaces()+bFacei;
            label own = mesh_.faceOwner()[facei];

            if
            (
                (cellLevel[own] == neiLevel[bFacei])
             && (
                    (refineCell[own] != -1)
                 != (neiRefineCell[bFacei] != -1)
                )
            )
            {
                isOutsideFace.set(facei);
            }
        }
    }
}


Foam::label Foam::meshRefinement::countFaceDirs
(
    const bitSet& isOutsideFace,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];
    const vectorField& faceAreas = mesh_.faceAreas();

    Vector<bool> haveDirs(vector::uniform(false));

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];

        if (isOutsideFace[facei])
        {
            const vector& n = faceAreas[facei];
            scalar magSqrN = magSqr(n);

            if (magSqrN > ROOTVSMALL)
            {
                for
                (
                    direction dir = 0;
                    dir < pTraits<vector>::nComponents;
                    dir++
                )
                {
                    if (Foam::sqr(n[dir]) > 0.99*magSqrN)
                    {
                        haveDirs[dir] = true;
                        break;
                    }
                }
            }
        }
    }

    label nDirs = 0;
    forAll(haveDirs, dir)
    {
        if (haveDirs[dir])
        {
            nDirs++;
        }
    }
    return nDirs;
}


void Foam::meshRefinement::growSet
(
    const labelList& neiLevel,
    const bitSet& isOutsideFace,
    labelList& refineCell,
    label& nRefine
) const
{
    // Get cells with three or more outside faces
    const cellList& cells = mesh_.cells();
    forAll(cells, celli)
    {
        if (refineCell[celli] == -1)
        {
            if (countFaceDirs(isOutsideFace, celli) == 3)
            {
                // Mark cell with any value
                refineCell[celli] = 0;
                nRefine++;
            }
        }
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::removeGapCells
(
    const scalar planarAngle,
    const labelList& minSurfaceLevel,
    const labelList& globalToMasterPatch,
    const label growIter
)
{
    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList refineCell(mesh_.nCells(), -1);
    label nRefine = 0;
    markProximityRefinement
    (
        Foam::cos(degToRad(planarAngle)),

        minSurfaceLevel,                                // surface min level
        labelList(minSurfaceLevel.size(), labelMax),    // surfaceGapLevel

        labelMax/Pstream::nProcs(), //nAllowRefine,
        neiLevel,
        neiCc,

        refineCell,
        nRefine
    );


    //// Mark big-gap refinement
    //markFakeGapRefinement
    //(
    //    Foam::cos(degToRad(planarAngle)),
    //
    //    labelMax/Pstream::nProcs(), //nAllowRefine,
    //    neiLevel,
    //    neiCc,
    //
    //    refineCell,
    //    nRefine
    //);


    Info<< "Marked for blocking due to close opposite surfaces         : "
        << returnReduce(nRefine, sumOp<label>()) << " cells." << endl;

    // Remove outliers, i.e. cells with all points exposed
    if (growIter)
    {
        labelList oldRefineCell(refineCell);

        // Pass1: extend the set to fill in gaps
        bitSet isOutsideFace;
        for (label iter = 0; iter < growIter; iter++)
        {
            // Get outside faces
            markOutsideFaces
            (
                meshCutter_.cellLevel(),
                neiLevel,
                refineCell,
                isOutsideFace
            );
            // Extend with cells with three outside faces
            growSet(neiLevel, isOutsideFace, refineCell, nRefine);
        }


        // Pass2: erode back to original set if pass1 didn't help
        for (label iter = 0; iter < growIter; iter++)
        {
            // Get outside faces. Ignore cell level.
            markOutsideFaces
            (
                labelList(mesh_.nCells(), 0),
                labelList(neiLevel.size(), 0),
                refineCell,
                isOutsideFace
            );

            // Unmark cells with three or more outside faces
            for (label celli = 0; celli < mesh_.nCells(); celli++)
            {
                if (refineCell[celli] != -1 && oldRefineCell[celli] == -1)
                {
                    if (countFaceDirs(isOutsideFace, celli) >= 3)
                    {
                        refineCell[celli] = -1;
                        --nRefine;
                    }
                }
            }
        }

        Info<< "Marked for blocking after filtering                        : "
            << returnReduce(nRefine, sumOp<label>()) << " cells." << endl;
    }


    // Determine patch for every mesh face
    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    labelList unnamedSurfaces(surfaceZonesInfo::getUnnamedSurfaces(surfZones));
    const label defaultRegion(surfaces_.globalRegion(unnamedSurfaces[0], 0));

    const labelList nearestRegion
    (
        nearestIntersection
        (
            unnamedSurfaces,
            defaultRegion
        )
    );

    // Pack
    labelList cellsToRemove(nRefine);
    nRefine = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell[cellI] != -1)
        {
            cellsToRemove[nRefine++] = cellI;
        }
    }

    // Remove cells
    removeCells cellRemover(mesh_);
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));

    labelList exposedPatches(exposedFaces.size());
    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        exposedPatches[i] = globalToMasterPatch[nearestRegion[facei]];
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


// ************************************************************************* //
