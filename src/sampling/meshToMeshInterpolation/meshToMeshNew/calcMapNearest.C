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
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshToMeshNew::calcMapNearest
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
    List<DynamicList<label> > srcToTgt(src.nCells());
    List<DynamicList<label> > tgtToSrc(tgt.nCells());

    const scalarField& srcVc = src.cellVolumes();
    const scalarField& tgtVc = tgt.cellVolumes();

    label srcCellI = srcSeedI;
    label tgtCellI = tgtSeedI;

    boolList tgtMapFlag(tgt.nCells(), true);

    do
    {
        // find nearest tgt cell
        findNearestCell(src, tgt, srcCellI, tgtCellI);

        // store src/tgt cell pair
        srcToTgt[srcCellI].append(tgtCellI);
        tgtToSrc[tgtCellI].append(srcCellI);

        // mark source cell srcCellI and tgtCellI as matched
        mapFlag[srcCellI] = false;
        tgtMapFlag[tgtCellI] = false;

        // accumulate intersection volume
        V_ += srcVc[srcCellI];

        // find new source cell
        setNextNearestCells
        (
            startSeedI,
            srcCellI,
            tgtCellI,
            mapFlag,
            src,
            tgt,
            srcCellIDs
        );
    }
    while (srcCellI >= 0);

    // If there are more target cells than source cells, some target cells
    // will not yet be mapped
    forAll(tgtMapFlag, tgtCellI)
    {
        if (tgtMapFlag[tgtCellI])
        {
            label srcCellI = findMappedSrcCell(tgt, tgtCellI, tgtToSrc);

            findNearestCell(tgt, src, tgtCellI, srcCellI);

            tgtToSrc[tgtCellI].append(srcCellI);
        }
    }

    // transfer addressing into persistent storage
    // note: always 1 target cell per source cell (srcToTgt)
    //       can be multiple source cells per target cell (tgtToSrc)
    forAll(srcToTgtCellAddr_, i)
    {
        scalar v = srcVc[i];
        srcToTgtCellWght_[i] = scalarList(srcToTgt[i].size(), v);
        srcToTgtCellAddr_[i].transfer(srcToTgt[i]);
    }

    forAll(tgtToSrcCellAddr_, i)
    {
        scalar v = tgtVc[i];
        tgtToSrcCellWght_[i] = scalarList(tgtToSrc[i].size(), v);
        tgtToSrcCellAddr_[i].transfer(tgtToSrc[i]);
    }
}


void Foam::meshToMeshNew::findNearestCell
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcCellI,
    label& tgtCellI
)
{
    const vectorField& srcC = src.cellCentres();
    const vectorField& tgtC = tgt.cellCentres();

    const vector& srcP = srcC[srcCellI];

    DynamicList<label> tgtCells(10);
    tgtCells.append(tgtCellI);

    DynamicList<label> visitedCells(10);

    scalar d = GREAT;

    do
    {
        label tgtI = tgtCells.remove();
        visitedCells.append(tgtI);

        scalar dTest = magSqr(tgtC[tgtI] - srcP);
        if (dTest < d)
        {
            tgtCellI = tgtI;
            d = dTest;
            appendNbrCells(tgtCellI, tgt, visitedCells, tgtCells);
        }

    } while (tgtCells.size() > 0);
}


void Foam::meshToMeshNew::setNextNearestCells
(
    label& startSeedI,
    label& srcCellI,
    label& tgtCellI,
    boolList& mapFlag,
    const polyMesh& src,
    const polyMesh& tgt,
    const labelList& srcCellIDs
)
{
    const labelList& srcNbr = src.cellCells()[srcCellI];

    srcCellI = -1;
    forAll(srcNbr, i)
    {
        label cellI = srcNbr[i];
        if (mapFlag[cellI])
        {
            srcCellI = cellI;
            startSeedI = cellI + 1;

            return;
        }
    }

    (void)findInitialSeeds
    (
        src,
        tgt,
        srcCellIDs,
        mapFlag,
        startSeedI,
        srcCellI,
        tgtCellI
    );
}


Foam::label Foam::meshToMeshNew::findMappedSrcCell
(
    const polyMesh& tgt,
    const label tgtCellI,
    const List<DynamicList<label> >& tgtToSrc
) const
{
    DynamicList<label> testCells(10);
    DynamicList<label> visitedCells(10);

    testCells.append(tgtCellI);

    do
    {
        // search target tgtCellI neighbours for match with source cell
        label tgtI = testCells.remove();

        if (findIndex(visitedCells, tgtI) == -1)
        {
            visitedCells.append(tgtI);

            if (tgtToSrc[tgtI].size())
            {
                return tgtToSrc[tgtI][0];
            }
            else
            {
                const labelList& nbrCells = tgt.cellCells()[tgtI];

                forAll(nbrCells, i)
                {
                    if (findIndex(visitedCells, nbrCells[i]) == -1)
                    {
                        testCells.append(nbrCells[i]);
                    }
                }
            }
        }
    } while (testCells.size());

    // did not find any match - should not be possible to get here!
    return -1;
}


// ************************************************************************* //
