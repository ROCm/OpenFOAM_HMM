/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "shortestPathSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "topoDistanceData.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shortestPathSet, 0);
    addToRunTimeSelectionTable(sampledSet, shortestPathSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::shortestPathSet::findMinFace
(
    const polyMesh& mesh,
    const label cellI,
    const List<topoDistanceData>& allFaceInfo,
    const point& origin
)
{
    const cell& cFaces2 = mesh.cells()[cellI];

    // 1. Get topologically nearest face

    label minDist = labelMax;
    label minFaceI = -1;
    forAll(cFaces2, i)
    {
        label faceI = cFaces2[i];
        const topoDistanceData& info = allFaceInfo[faceI];
        if (info.distance() < minDist)
        {
            minDist = info.distance();
            minFaceI = faceI;
        }
    }

    // 2. Check all faces with minDist for minimum distance to origin
    scalar minDist2 = ROOTVGREAT;
    forAll(cFaces2, i)
    {
        label faceI = cFaces2[i];
        if (allFaceInfo[faceI].distance() == minDist)
        {
            scalar d2 = magSqr(mesh.faceCentres()[faceI]-origin);
            if (d2 < minDist2)
            {
                minDist2 = d2;
                minFaceI  = faceI;
            }
        }
    }

    return minFaceI;
}


void Foam::shortestPathSet::genSamples(const polyMesh& mesh)
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    forAll(insidePoints_, pointI)
    {
        label cell1I = mesh.findCell(insidePoints_[pointI]);

        //
        // Pass1: Set initial changed faces from cell1 (seed)
        //

        List<topoDistanceData> faceDist;
        labelList cFaces1;

        if (cell1I != -1)
        {
            cFaces1 = mesh.cells()[cell1I];
            faceDist.setSize
            (
                cFaces1.size(),
                topoDistanceData(123, 0)
            );
        }

        List<topoDistanceData> allFaceInfo(mesh.nFaces());
        List<topoDistanceData> allCellInfo(mesh.nCells());

        // Walk through face-cell wave till all cells are reached
        FaceCellWave
        <
            topoDistanceData
        > wallDistCalc
        (
            mesh,
            cFaces1,
            faceDist,
            allFaceInfo,
            allCellInfo,
            mesh.globalData().nTotalCells()+1   // max iterations
        );

        // Pass2: walk from outside points backwards. Note: could be done using
        //        FaceCellWave as well but is overly complex since does
        //        not allow logic comparing all faces of a cell.

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        // Get the target point
        label cell2I = mesh.findCell(outsidePoints_[pointI]);

        // The number of cells between cell1 and cell2 is the max number of
        // iterations to search backward
        label nPathPoints = 0;

        if (cell2I != -1)
        {
            if (!allCellInfo[cell2I].valid(wallDistCalc.data()))
            {
                WarningInFunction
                    << "Point " << outsidePoints_[pointI]
                    << " not reachable by walk. Probably mesh has "
                    << " island/regions. Skipped route detection." << endl;
                return;
            }

            nPathPoints = allCellInfo[cell2I].distance();
        }
        reduce(nPathPoints, maxOp<label>());

        // Start with given target cell and walk back
        label frontCellI = cell2I;

        while (nPathPoints--)
        {
            label frontFaceI = -1;

            // Work within same processor
            if (frontCellI != -1)
            {
                // Find face with lowest distance from seed
                //   x  |  x  2  1  2  2  |  x  |  x
                //  --- + --- + -1- + -2- + --- + ---
                //   x  |  1  1  0  1  1  |  x  |  x
                //  --- + --- + -1- + -2- + --- + ---
                //   x  |  x  2  1  2  2  3  3  4  4
                //  --- + --- + --- + -3- + -4- + -5-
                //   x  |  x  3  2  3  3  4  4  5  5
                // e.g. if we start from cell with value = 4, we have neighbour
                // faces 4, 4, 5, 5. Choose 4 (least distance to seed)
                // and continue...

                frontFaceI = findMinFace
                (
                    mesh,
                    frontCellI,
                    allFaceInfo,
                    outsidePoints_[pointI]
                );

                // Loop until we hit a boundary face
                while (mesh.isInternalFace(frontFaceI))
                {
                    // Step to neighbouring cell
                    label nbrCellI = mesh.faceOwner()[frontFaceI];
                    if (nbrCellI == frontCellI)
                    {
                        nbrCellI = mesh.faceNeighbour()[frontFaceI];
                    }

                    if (nbrCellI == cell1I)
                    {
                        // Pout<< " Found connection seed cell!" << endl;
                        frontCellI = -1;
                        break;
                    }

                    frontCellI = nbrCellI;

                    // Pick best face on cell
                    frontFaceI = findMinFace
                    (
                        mesh,
                        frontCellI,
                        allFaceInfo,
                        outsidePoints_[pointI]
                    );

                    // Set the sampling point
                    samplingPts.append(mesh.cellCentres()[frontCellI]);
                    samplingCells.append(frontCellI);
                    samplingFaces.append(-1);
                    samplingSegments.append(pointI);
                    //Check if mag of distance is useful
                    samplingCurveDist.append(nPathPoints);
                }
            }

            // Situation 1: we found the destination cell (do nothing)
            if (!returnReduce(frontCellI != -1, orOp<bool>()))
            {
                break;
            }

            // Situation 2: we're on a coupled patch and might need to
            //              switch processor/cell
            boolList isFront(mesh.nFaces()-mesh.nInternalFaces(), false);

            if (frontFaceI != -1)
            {
                isFront[frontFaceI-mesh.nInternalFaces()] = true;
            }
            syncTools::swapBoundaryFaceList(mesh, isFront);

            frontCellI = -1;
            forAll(pbm, patchI)
            {
                const polyPatch& pp = pbm[patchI];
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;
                    if (isFront[faceI-mesh.nInternalFaces()])
                    {
                        frontCellI = pp.faceCells()[i];
                        break;
                    }
                }

                if (frontCellI != -1)
                {
                    samplingPts.append(mesh.cellCentres()[frontCellI]);
                    samplingCells.append(frontCellI);
                    samplingFaces.append(-1);
                    samplingSegments.append(pointI);
                    samplingCurveDist.append(nPathPoints);
                    break;
                }
            }
        }
    }
    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shortestPathSet::shortestPathSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const pointField& insidePoints,
    const pointField& outsidePoints
)
:
    sampledSet(name, mesh, searchEngine, axis),
    insidePoints_(insidePoints),
    outsidePoints_(outsidePoints)
{
    genSamples(mesh);

    if (debug)
    {
        write(Info);
    }
}


Foam::shortestPathSet::shortestPathSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    insidePoints_(dict.lookup("insidePoints")),
    outsidePoints_(dict.lookup("outsidePoints"))
{
    genSamples(mesh);

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shortestPathSet::~shortestPathSet()
{}


// ************************************************************************* //
