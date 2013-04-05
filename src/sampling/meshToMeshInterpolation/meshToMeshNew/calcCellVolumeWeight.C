/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "meshToMeshNew.H"
#include "tetOverlapVolume.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshToMeshNew::calcCellVolumeWeight
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcSeedI,
    const label tgtSeedI,
    const labelList& srcCellIDs,
    boolList& mapFlag,
    label& startSeedI
)
{
    label srcCellI = srcSeedI;
    label tgtCellI = tgtSeedI;

    List<DynamicList<label> > srcToTgtAddr(src.nCells());
    List<DynamicList<scalar> > srcToTgtWght(src.nCells());

    List<DynamicList<label> > tgtToSrcAddr(tgt.nCells());
    List<DynamicList<scalar> > tgtToSrcWght(tgt.nCells());

    // list of tgt cell neighbour cells
    DynamicList<label> nbrTgtCells(10);

    // list of tgt cells currently visited for srcCellI to avoid multiple hits
    DynamicList<label> visitedTgtCells(10);

    // list to keep track of tgt cells used to seed src cells
    labelList seedCells(src.nCells(), -1);
    seedCells[srcCellI] = tgtCellI;

    const scalarField& srcVol = src.cellVolumes();

    do
    {
        nbrTgtCells.clear();
        visitedTgtCells.clear();

        // append initial target cell and neighbours
        nbrTgtCells.append(tgtCellI);
        appendNbrCells(tgtCellI, tgt, visitedTgtCells, nbrTgtCells);

        do
        {
            tgtCellI = nbrTgtCells.remove();
            visitedTgtCells.append(tgtCellI);

            scalar vol = interVol(src, tgt, srcCellI, tgtCellI);

            // accumulate addressing and weights for valid intersection
            if (vol/srcVol[srcCellI] > tolerance_)
            {
                // store src/tgt cell pair
                srcToTgtAddr[srcCellI].append(tgtCellI);
                srcToTgtWght[srcCellI].append(vol);

                tgtToSrcAddr[tgtCellI].append(srcCellI);
                tgtToSrcWght[tgtCellI].append(vol);

                appendNbrCells(tgtCellI, tgt, visitedTgtCells, nbrTgtCells);

                // accumulate intersection volume
                V_ += vol;
            }
        }
        while (!nbrTgtCells.empty());

        mapFlag[srcCellI] = false;

        // find new source seed cell
        setNextCells
        (
            startSeedI,
            srcCellI,
            tgtCellI,
            src,
            tgt,
            srcCellIDs,
            mapFlag,
            visitedTgtCells,
            seedCells
        );
    }
    while (srcCellI != -1);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr_, i)
    {
        srcToTgtCellAddr_[i].transfer(srcToTgtAddr[i]);
        srcToTgtCellWght_[i].transfer(srcToTgtWght[i]);
    }

    forAll(tgtToSrcCellAddr_, i)
    {
        tgtToSrcCellAddr_[i].transfer(tgtToSrcAddr[i]);
        tgtToSrcCellWght_[i].transfer(tgtToSrcWght[i]);
    }
}


void Foam::meshToMeshNew::setNextCells
(
    label& startSeedI,
    label& srcCellI,
    label& tgtCellI,
    const polyMesh& src,
    const polyMesh& tgt,
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const DynamicList<label>& visitedCells,
    labelList& seedCells
) const
{
    const labelList& srcNbrCells = src.cellCells()[srcCellI];

    // set possible seeds for later use by querying all src cell neighbours
    // with all visited target cells
    bool valuesSet = false;
    forAll(srcNbrCells, i)
    {
        label cellS = srcNbrCells[i];

        if (mapFlag[cellS] && seedCells[cellS] == -1)
        {
            forAll(visitedCells, j)
            {
                label cellT = visitedCells[j];

                if (intersect(src, tgt, cellS, cellT))
                {
                    seedCells[cellS] = cellT;

                    if (!valuesSet)
                    {
                        srcCellI = cellS;
                        tgtCellI = cellT;
                        valuesSet = true;
                    }
                }
            }
        }
    }

    // set next src and tgt cells if not set above
    if (valuesSet)
    {
        return;
    }
    else
    {
        // try to use existing seed
        bool foundNextSeed = false;
        for (label i = startSeedI; i < srcCellIDs.size(); i++)
        {
            label cellS = srcCellIDs[i];

            if (mapFlag[cellS])
            {
                if (!foundNextSeed)
                {
                    startSeedI = i;
                    foundNextSeed = true;
                }

                if (seedCells[cellS] != -1)
                {
                    srcCellI = cellS;
                    tgtCellI = seedCells[cellS];

                    return;
                }
            }
        }

        // perform new search to find match
        if (debug)
        {
            Pout<< "Advancing front stalled: searching for new "
                << "target cell" << endl;
        }

        bool restart =
            findInitialSeeds
            (
                src,
                tgt,
                srcCellIDs,
                mapFlag,
                startSeedI,
                srcCellI,
                tgtCellI
            );

        if (restart)
        {
            // successfully found new starting seed-pair
            return;
        }
    }

    // if we have got to here, there are no more src/tgt cell intersections
    srcCellI = -1;
    tgtCellI = -1;
}


bool Foam::meshToMeshNew::intersect
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcCellI,
    const label tgtCellI
) const
{
    scalar threshold = tolerance_*src.cellVolumes()[srcCellI];

    tetOverlapVolume overlapEngine;

    treeBoundBox bbTgtCell
    (
        pointField
        (
            tgt.points(),
            tgt.cellPoints()[tgtCellI]
        )
    );

    return overlapEngine.cellCellOverlapMinDecomp
    (
        src,
        srcCellI,
        tgt,
        tgtCellI,
        bbTgtCell,
        threshold
    );
}


Foam::scalar Foam::meshToMeshNew::interVol
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcCellI,
    const label tgtCellI
) const
{
    tetOverlapVolume overlapEngine;

    treeBoundBox bbTgtCell
    (
        pointField
        (
            tgt.points(),
            tgt.cellPoints()[tgtCellI]
        )
    );

    scalar vol = overlapEngine.cellCellOverlapVolumeMinDecomp
    (
        src,
        srcCellI,
        tgt,
        tgtCellI,
        bbTgtCell
    );

    return vol;
}


// ************************************************************************* //
