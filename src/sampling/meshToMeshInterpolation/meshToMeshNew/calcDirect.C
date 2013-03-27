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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshToMeshNew::calcDirect
(
    const polyMesh& src,
    const polyMesh& tgt,
    const label srcSeedI,
    const label tgtSeedI
)
{
    // store a list of src cells already mapped
    boolList srcSeedFlag(src.nCells(), true);
    labelList srcTgtSeed(src.nCells(), -1);

    List<DynamicList<label> > srcToTgt(src.nCells());
    List<DynamicList<label> > tgtToSrc(tgt.nCells());

    DynamicList<label> srcSeeds;

    const scalarField& srcVc = src.cellVolumes();
    const scalarField& tgtVc = tgt.cellVolumes();

    label srcCellI = srcSeedI;
    label tgtCellI = tgtSeedI;

    do
    {
        // store src/tgt cell pair
        srcToTgt[srcCellI].append(tgtCellI);
        tgtToSrc[tgtCellI].append(srcCellI);

        // mark source cell srcSeedI as matched
        srcSeedFlag[srcCellI] = false;

        // accumulate intersection volume
        V_ += srcVc[srcCellI];

        // find new source seed cell
        appendToDirectSeeds
        (
            src,
            tgt,
            srcSeedFlag,
            srcTgtSeed,
            srcSeeds,
            srcCellI,
            tgtCellI
        );
    }
    while (srcCellI >= 0);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr_, i)
    {
        scalar v = srcVc[i];
        srcToTgtCellAddr_[i].transfer(srcToTgt[i]);
        srcToTgtCellWght_[i] = scalarList(srcToTgtCellAddr_[i].size(), v);
    }

    forAll(tgtToSrcCellAddr_, i)
    {
        scalar v = tgtVc[i];
        tgtToSrcCellAddr_[i].transfer(tgtToSrc[i]);
        tgtToSrcCellWght_[i] = scalarList(tgtToSrcCellAddr_[i].size(), v);
    }
}


void Foam::meshToMeshNew::appendToDirectSeeds
(
    const polyMesh& src,
    const polyMesh& tgt,
    boolList& mapFlag,
    labelList& srcTgtSeed,
    DynamicList<label>& srcSeeds,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const labelList& srcNbr = src.cellCells()[srcSeedI];
    const labelList& tgtNbr = tgt.cellCells()[tgtSeedI];

    const vectorField& srcCentre = src.cellCentres();

    forAll(srcNbr, i)
    {
        label srcI = srcNbr[i];

        if (mapFlag[srcI] && (srcTgtSeed[srcI] == -1))
        {
            // source cell srcI not yet mapped

            // identfy if target cell exists for source cell srcI
            bool found = false;
            forAll(tgtNbr, j)
            {
                label tgtI = tgtNbr[j];

                if (tgt.pointInCell(srcCentre[srcI], tgtI))
                {
                    // new match - append to lists
                    found = true;

                    srcTgtSeed[srcI] = tgtI;
                    srcSeeds.append(srcI);

                    break;
                }
            }

            if (!found)
            {
                // no match available for source cell srcI
                mapFlag[srcI] = false;
            }
        }
    }

    if (srcSeeds.size())
    {
        srcSeedI = srcSeeds.remove();
        tgtSeedI = srcTgtSeed[srcSeedI];
    }
    else
    {
        srcSeedI = -1;
        tgtSeedI = -1;
    }
}


// ************************************************************************* //
