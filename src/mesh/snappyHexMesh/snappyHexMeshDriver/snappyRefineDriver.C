/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "snappyRefineDriver.H"
#include "meshRefinement.H"
#include "fvMesh.H"
#include "Time.H"
#include "cellSet.H"
#include "syncTools.H"
#include "refinementParameters.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "mapDistributePolyMesh.H"
#include "unitConversion.H"
#include "snapParameters.H"
#include "localPointRegion.H"
#include "IOmanip.H"
#include "labelVector.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(snappyRefineDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::snappyRefineDriver::snappyRefineDriver
(
    meshRefinement& meshRefiner,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
)
:
    meshRefiner_(meshRefiner),
    decomposer_(decomposer),
    distributor_(distributor),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::snappyRefineDriver::featureEdgeRefine
(
    const refinementParameters& refineParams,
    const label maxIter,
    const label minRefine
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(edge, "snappyHexMesh::refine::edge");
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    if (meshRefiner_.features().size() && maxIter > 0)
    {
        for (; iter < maxIter; iter++)
        {
            Info<< nl
                << "Feature refinement iteration " << iter << nl
                << "------------------------------" << nl
                << endl;

            labelList candidateCells
            (
                meshRefiner_.refineCandidates
                (
                    refineParams.locationsInMesh(),
                    refineParams.curvature(),
                    refineParams.planarAngle(),

                    true,               // featureRefinement
                    false,              // featureDistanceRefinement
                    false,              // internalRefinement
                    false,              // surfaceRefinement
                    false,              // curvatureRefinement
                    false,              // smallFeatureRefinement
                    false,              // gapRefinement
                    false,              // bigGapRefinement
                    false,              // spreadGapSize
                    refineParams.maxGlobalCells(),
                    refineParams.maxLocalCells()
                )
            );
            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );
            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;



            label nCellsToRefine = cellsToRefine.size();
            reduce(nCellsToRefine, sumOp<label>());

            Info<< "Selected for feature refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            if (nCellsToRefine <= minRefine)
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }


            if (debug > 0)
            {
                const_cast<Time&>(mesh.time())++;
            }


            if
            (
                returnReduce
                (
                    (mesh.nCells() >= refineParams.maxLocalCells()),
                    orOp<bool>()
                )
            )
            {
                meshRefiner_.balanceAndRefine
                (
                    "feature refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance()
                );
            }
            else
            {
                meshRefiner_.refineAndBalance
                (
                    "feature refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance()
                );
            }
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::smallFeatureRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(feature, "snappyHexMesh::refine::smallFeature");
    const fvMesh& mesh = meshRefiner_.mesh();


    label iter = 0;

    // See if any surface has an extendedGapLevel
    labelList surfaceMaxLevel(meshRefiner_.surfaces().maxGapLevel());
    labelList shellMaxLevel(meshRefiner_.shells().maxGapLevel());

    if (max(surfaceMaxLevel) == 0 && max(shellMaxLevel) == 0)
    {
        return iter;
    }

    for (; iter < maxIter; iter++)
    {
        Info<< nl
            << "Small surface feature refinement iteration " << iter << nl
            << "--------------------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                true,               // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false,              // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if (nCellsToRefine == 0)
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "small feature refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "small feature refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::surfaceOnlyRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(surface, "snappyHexMesh::refine::surface");
    const fvMesh& mesh = meshRefiner_.mesh();

    // Determine the maximum refinement level over all surfaces. This
    // determines the minumum number of surface refinement iterations.
    label overallMaxLevel = max(meshRefiner_.surfaces().maxLevel());

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Surface refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                true,               // surfaceRefinement
                true,               // curvatureRefinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false,              // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );
        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::gapOnlyRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    const fvMesh& mesh = meshRefiner_.mesh();

    // Determine the maximum refinement level over all surfaces. This
    // determines the minumum number of surface refinement iterations.

    label maxIncrement = 0;
    const labelList& maxLevel = meshRefiner_.surfaces().maxLevel();
    const labelList& gapLevel = meshRefiner_.surfaces().gapLevel();

    forAll(maxLevel, i)
    {
        maxIncrement = max(maxIncrement, gapLevel[i]-maxLevel[i]);
    }

    label iter = 0;

    if (maxIncrement == 0)
    {
        return iter;
    }

    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Gap refinement iteration " << iter << nl
            << "--------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // smallFeatureRefinement
                true,               // gapRefinement
                false,              // bigGapRefinement
                false,              // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromGap." << endl;
            cellSet c(mesh, "candidateCellsFromGap", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        // Grow by one layer to make sure we're covering the gap
        {
            boolList isCandidateCell(mesh.nCells(), false);
            forAll(candidateCells, i)
            {
                isCandidateCell[candidateCells[i]] = true;
            }

            for (label i=0; i<1; i++)
            {
                boolList newIsCandidateCell(isCandidateCell);

                // Internal faces
                for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
                {
                    label own = mesh.faceOwner()[facei];
                    label nei = mesh.faceNeighbour()[facei];

                    if (isCandidateCell[own] != isCandidateCell[nei])
                    {
                        newIsCandidateCell[own] = true;
                        newIsCandidateCell[nei] = true;
                    }
                }

                // Get coupled boundary condition values
                boolList neiIsCandidateCell;
                syncTools::swapBoundaryCellList
                (
                    mesh,
                    isCandidateCell,
                    neiIsCandidateCell
                );

                // Boundary faces
                for
                (
                    label facei = mesh.nInternalFaces();
                    facei < mesh.nFaces();
                    facei++
                )
                {
                    label own = mesh.faceOwner()[facei];
                    label bFacei = facei-mesh.nInternalFaces();

                    if (isCandidateCell[own] != neiIsCandidateCell[bFacei])
                    {
                        newIsCandidateCell[own] = true;
                    }
                }

                isCandidateCell.transfer(newIsCandidateCell);
            }

            label n = 0;
            forAll(isCandidateCell, celli)
            {
                if (isCandidateCell[celli])
                {
                    n++;
                }
            }
            candidateCells.setSize(n);
            n = 0;
            forAll(isCandidateCell, celli)
            {
                if (isCandidateCell[celli])
                {
                    candidateCells[n++] = celli;
                }
            }
        }


        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromGapPlusBuffer." << endl;
            cellSet c(mesh, "candidateCellsFromGapPlusBuffer", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }


        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= maxIncrement
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::bigGapOnlyRefine
(
    const refinementParameters& refineParams,
    const bool spreadGapSize,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    // See if any surface has an extendedGapLevel
    labelList surfaceMaxLevel(meshRefiner_.surfaces().maxGapLevel());
    labelList shellMaxLevel(meshRefiner_.shells().maxGapLevel());

    label overallMaxLevel(max(max(surfaceMaxLevel), max(shellMaxLevel)));

    if (overallMaxLevel == 0)
    {
        return iter;
    }


    for (; iter < maxIter; iter++)
    {
        Info<< nl
            << "Big gap refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                true,               // bigGapRefinement
                spreadGapSize,      // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );


        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromBigGap." << endl;
            cellSet c(mesh, "candidateCellsFromBigGap", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "big gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "big gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::danglingCellRefine
(
    const refinementParameters& refineParams,
    const label nFaces,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(dangling, "snappyHexMesh::refine::danglingCell");
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Dangling coarse cells refinement iteration " << iter << nl
            << "--------------------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        const cellList& cells = mesh.cells();
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        labelList candidateCells;
        {
            cellSet candidateCellSet(mesh, "candidateCells", cells.size()/1000);

            forAll(cells, celli)
            {
                const cell& cFaces = cells[celli];

                label nIntFaces = 0;
                forAll(cFaces, i)
                {
                    label bFacei = cFaces[i]-mesh.nInternalFaces();
                    if (bFacei < 0)
                    {
                        nIntFaces++;
                    }
                    else
                    {
                        label patchi = pbm.patchID()[bFacei];
                        if (pbm[patchi].coupled())
                        {
                            nIntFaces++;
                        }
                    }
                }

                if (nIntFaces == nFaces)
                {
                    candidateCellSet.insert(celli);
                }
            }

            if (debug&meshRefinement::MESH)
            {
                Pout<< "Dumping " << candidateCellSet.size()
                    << " cells to cellSet candidateCellSet." << endl;
                candidateCellSet.instance() = meshRefiner_.timeName();
                candidateCellSet.write();
            }
            candidateCells = candidateCellSet.toc();
        }



        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine. After a few iterations check if too
        // few cells
        if
        (
            nCellsToRefine == 0
         || (
                iter >= 1
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "coarse cell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "coarse cell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
    }
    return iter;
}


// Detect cells with opposing intersected faces of differing refinement
// level and refine them.
Foam::label Foam::snappyRefineDriver::refinementInterfaceRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(interface, "snappyHexMesh::refine::transition");
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    if (refineParams.interfaceRefine())
    {
        for (;iter < maxIter; iter++)
        {
            Info<< nl
                << "Refinement transition refinement iteration " << iter << nl
                << "--------------------------------------------" << nl
                << endl;

            const labelList& surfaceIndex = meshRefiner_.surfaceIndex();
            const hexRef8& cutter = meshRefiner_.meshCutter();
            const vectorField& fA = mesh.faceAreas();
            const labelList& faceOwner = mesh.faceOwner();


            // Determine cells to refine
            // ~~~~~~~~~~~~~~~~~~~~~~~~~

            const cellList& cells = mesh.cells();

            labelList candidateCells;
            {
                // Pass1: pick up cells with differing face level

                cellSet transitionCells
                (
                    mesh,
                    "transitionCells",
                    cells.size()/100
                );

                forAll(cells, celli)
                {
                    const cell& cFaces = cells[celli];
                    label cLevel = cutter.cellLevel()[celli];

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        if (surfaceIndex[facei] != -1)
                        {
                            label fLevel = cutter.faceLevel(facei);
                            if (fLevel != cLevel)
                            {
                                transitionCells.insert(celli);
                            }
                        }
                    }
                }


                cellSet candidateCellSet
                (
                    mesh,
                    "candidateCells",
                    cells.size()/1000
                );

                // Pass2: check for oppositeness

                //forAllConstIter(cellSet, transitionCells, iter)
                //{
                //    label celli = iter.key();
                //    const cell& cFaces = cells[celli];
                //    const point& cc = cellCentres[celli];
                //    const scalar rCVol = pow(cellVolumes[celli], -5.0/3.0);
                //
                //    // Determine principal axes of cell
                //    symmTensor R(Zero);
                //
                //    forAll(cFaces, i)
                //    {
                //        label facei = cFaces[i];
                //
                //        const point& fc = faceCentres[facei];
                //
                //        // Calculate face-pyramid volume
                //        scalar pyrVol = 1.0/3.0 * fA[facei] & (fc-cc);
                //
                //        if (faceOwner[facei] != celli)
                //        {
                //            pyrVol = -pyrVol;
                //        }
                //
                //        // Calculate face-pyramid centre
                //        vector pc = (3.0/4.0)*fc + (1.0/4.0)*cc;
                //
                //        R += pyrVol*sqr(pc-cc)*rCVol;
                //    }
                //
                //    //- MEJ: Problem: truncation errors cause complex evs
                //    vector lambdas(eigenValues(R));
                //    const tensor axes(eigenVectors(R, lambdas));
                //
                //
                //    // Check if this cell has
                //    // - opposing sides intersected
                //    // - which are of different refinement level
                //    // - plus the inbetween face
                //
                //    labelVector plusFaceLevel(labelVector(-1, -1, -1));
                //    labelVector minFaceLevel(labelVector(-1, -1, -1));
                //
                //    forAll(cFaces, cFacei)
                //    {
                //        label facei = cFaces[cFacei];
                //
                //        if (surfaceIndex[facei] != -1)
                //        {
                //            label fLevel = cutter.faceLevel(facei);
                //
                //            // Get outwards pointing normal
                //            vector n = fA[facei]/mag(fA[facei]);
                //            if (faceOwner[facei] != celli)
                //            {
                //                n = -n;
                //            }
                //
                //            // What is major direction and sign
                //            direction cmpt = vector::X;
                //            scalar maxComp = (n&axes.x());
                //
                //            scalar yComp = (n&axes.y());
                //            scalar zComp = (n&axes.z());
                //
                //            if (mag(yComp) > mag(maxComp))
                //            {
                //                maxComp = yComp;
                //                cmpt = vector::Y;
                //            }
                //
                //            if (mag(zComp) > mag(maxComp))
                //            {
                //                maxComp = zComp;
                //                cmpt = vector::Z;
                //            }
                //
                //            if (maxComp > 0)
                //            {
                //                plusFaceLevel[cmpt] = max
                //                (
                //                    plusFaceLevel[cmpt],
                //                    fLevel
                //                );
                //            }
                //            else
                //            {
                //                minFaceLevel[cmpt] = max
                //                (
                //                    minFaceLevel[cmpt],
                //                    fLevel
                //                );
                //            }
                //        }
                //    }
                //
                //    // Check if we picked up any opposite differing level
                //    for (direction dir = 0; dir < vector::nComponents; dir++)
                //    {
                //        if
                //        (
                //            plusFaceLevel[dir] != -1
                //         && minFaceLevel[dir] != -1
                //         && plusFaceLevel[dir] != minFaceLevel[dir]
                //        )
                //        {
                //            candidateCellSet.insert(celli);
                //        }
                //    }
                //}

                const scalar oppositeCos = Foam::cos(degToRad(135.0));

                forAllConstIter(cellSet, transitionCells, iter)
                {
                    label celli = iter.key();
                    const cell& cFaces = cells[celli];
                    label cLevel = cutter.cellLevel()[celli];

                    // Detect opposite intersection
                    bool foundOpposite = false;

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        if
                        (
                            surfaceIndex[facei] != -1
                         && cutter.faceLevel(facei) > cLevel
                        )
                        {
                            // Get outwards pointing normal
                            vector n = fA[facei]/mag(fA[facei]);
                            if (faceOwner[facei] != celli)
                            {
                                n = -n;
                            }

                            // Check for any opposite intersection
                            forAll(cFaces, cFaceI2)
                            {
                                label face2i = cFaces[cFaceI2];

                                if
                                (
                                    face2i != facei
                                 && surfaceIndex[face2i] != -1
                                )
                                {
                                    // Get outwards pointing normal
                                    vector n2 = fA[face2i]/mag(fA[face2i]);
                                    if (faceOwner[face2i] != celli)
                                    {
                                        n2 = -n2;
                                    }


                                    if ((n&n2) < oppositeCos)
                                    {
                                        foundOpposite = true;
                                        break;
                                    }
                                }
                            }

                            if (foundOpposite)
                            {
                                break;
                            }
                        }
                    }


                    if (foundOpposite)
                    {
                        candidateCellSet.insert(celli);
                    }
                }

                if (debug&meshRefinement::MESH)
                {
                    Pout<< "Dumping " << candidateCellSet.size()
                        << " cells to cellSet candidateCellSet." << endl;
                    candidateCellSet.instance() = meshRefiner_.timeName();
                    candidateCellSet.write();
                }
                candidateCells = candidateCellSet.toc();
            }



            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );
            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;


            label nCellsToRefine = cellsToRefine.size();
            reduce(nCellsToRefine, sumOp<label>());

            Info<< "Selected for refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            // Stop when no cells to refine. After a few iterations check if too
            // few cells
            if
            (
                nCellsToRefine == 0
             || (
                    iter >= 1
                 && nCellsToRefine <= refineParams.minRefineCells()
                )
            )
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }


            if (debug)
            {
                const_cast<Time&>(mesh.time())++;
            }


            if
            (
                returnReduce
                (
                    (mesh.nCells() >= refineParams.maxLocalCells()),
                    orOp<bool>()
                )
            )
            {
                meshRefiner_.balanceAndRefine
                (
                    "interface cell refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance()
                );
            }
            else
            {
                meshRefiner_.refineAndBalance
                (
                    "interface cell refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance()
                );
            }
        }
    }
    return iter;
}


void Foam::snappyRefineDriver::removeInsideCells
(
    const refinementParameters& refineParams,
    const label nBufferLayers
)
{
    Info<< nl
        << "Removing mesh beyond surface intersections" << nl
        << "------------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (debug)
    {
       const_cast<Time&>(mesh.time())++;
    }

    meshRefiner_.splitMesh
    (
        nBufferLayers,                  // nBufferLayers
        refineParams.nErodeCellZone(),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.locationsInMesh(),
        refineParams.zonesInMesh(),
        refineParams.locationsOutsideMesh()
    );

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;

        Pout<< "Writing subsetted mesh to time "
            << meshRefiner_.timeName() << endl;
        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner_.timeName()
        );
        Pout<< "Dumped mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


Foam::label Foam::snappyRefineDriver::shellRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(shell, "snappyHexMesh::refine::shell");
    const fvMesh& mesh = meshRefiner_.mesh();

    // Mark current boundary faces with 0. Have meshRefiner maintain them.
    meshRefiner_.userFaceData().setSize(1);

    // mark list to remove any refined faces
    meshRefiner_.userFaceData()[0].first() = meshRefinement::REMOVE;
    meshRefiner_.userFaceData()[0].second() = createWithValues<labelList>
    (
        mesh.nFaces(),
        -1,
        meshRefiner_.intersectedFaces(),
        0
    );

    // Determine the maximum refinement level over all volume refinement
    // regions. This determines the minumum number of shell refinement
    // iterations.
    label overallMaxShellLevel = meshRefiner_.shells().maxLevel();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Shell refinement iteration " << iter << nl
            << "----------------------------" << nl
            << endl;

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                true,               // featureDistanceRefinement
                true,               // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false,              // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromShells." << endl;

            cellSet c(mesh, "candidateCellsFromShells", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        // Problem choosing starting faces for bufferlayers (bFaces)
        //  - we can't use the current intersected boundary faces
        //    (intersectedFaces) since this grows indefinitely
        //  - if we use 0 faces we don't satisfy bufferLayers from the
        //    surface.
        //  - possibly we want to have bFaces only the initial set of faces
        //    and maintain the list while doing the refinement.
        labelList bFaces
        (
            findIndices(meshRefiner_.userFaceData()[0].second(), 0)
        );

        //Info<< "Collected boundary faces : "
        //    << returnReduce(bFaces.size(), sumOp<label>()) << endl;

        labelList cellsToRefine;

        if (refineParams.nBufferLayers() <= 2)
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentSlowRefinement
            (
                refineParams.nBufferLayers(),
                candidateCells,                     // cells to refine
                bFaces,                             // faces for nBufferLayers
                1,                                  // point difference
                meshRefiner_.intersectedPoints()    // points to check
            );
        }
        else
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentSlowRefinement2
            (
                refineParams.nBufferLayers(),
                candidateCells,                 // cells to refine
                bFaces                          // faces for nBufferLayers
            );
        }

        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for internal refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxShellLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "shell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "shell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance()
            );
        }
    }
    meshRefiner_.userFaceData().clear();

    return iter;
}


void Foam::snappyRefineDriver::baffleAndSplitMesh
(
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool handleSnapProblems,
    const dictionary& motionDict
)
{
    addProfiling(split, "snappyHexMesh::refine::splitting");
    Info<< nl
        << "Splitting mesh at surface intersections" << nl
        << "---------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (debug)
    {
       const_cast<Time&>(mesh.time())++;
    }

    // Introduce baffles at surface intersections. Note:
    // meshRefiment::surfaceIndex() will
    // be like boundary face from now on so not coupled anymore.
    meshRefiner_.baffleAndSplitMesh
    (
        handleSnapProblems,             // detect&remove potential snap problem

        // Snap problem cell detection
        snapParams,
        refineParams.useTopologicalSnapDetection(),
        false,                          // perpendicular edge connected cells
        scalarField(0),                 // per region perpendicular angle
        refineParams.nErodeCellZone(),

        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.locationsInMesh(),
        refineParams.zonesInMesh(),
        refineParams.locationsOutsideMesh()
    );


    if (!handleSnapProblems) // merge free standing baffles?
    {
        meshRefiner_.mergeFreeStandingBaffles
        (
            snapParams,
            refineParams.useTopologicalSnapDetection(),
            false,                  // perpendicular edge connected cells
            scalarField(0),         // per region perpendicular angle
            refineParams.planarAngle(),
            motionDict,
            const_cast<Time&>(mesh.time()),
            globalToMasterPatch_,
            globalToSlavePatch_,
            refineParams.locationsInMesh(),
            refineParams.locationsOutsideMesh()
        );
    }
}


void Foam::snappyRefineDriver::zonify
(
    const refinementParameters& refineParams,
    wordPairHashTable& zonesToFaceZone
)
{
    // Mesh is at its finest. Do zoning
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This puts all faces with intersection across a zoneable surface
    // into that surface's faceZone. All cells inside faceZone get given the
    // same cellZone.

    const labelList namedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(meshRefiner_.surfaces().surfZones());

    if
    (
        namedSurfaces.size()
     || refineParams.zonesInMesh().size()
    )
    {
        Info<< nl
            << "Introducing zones for interfaces" << nl
            << "--------------------------------" << nl
            << endl;

        const fvMesh& mesh = meshRefiner_.mesh();

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        meshRefiner_.zonify
        (
            refineParams.allowFreeStandingZoneFaces(),
            refineParams.nErodeCellZone(),
            refineParams.locationsInMesh(),
            refineParams.zonesInMesh(),
            zonesToFaceZone
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing zoned mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner_.timeName()
            );
        }

        // Check that all faces are synced
        meshRefinement::checkCoupledFaceZones(mesh);
    }
}


void Foam::snappyRefineDriver::splitAndMergeBaffles
(
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool handleSnapProblems,
    const dictionary& motionDict
)
{
    Info<< nl
        << "Handling cells with snap problems" << nl
        << "---------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    // Introduce baffles and split mesh
    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    const scalarField& perpAngle = meshRefiner_.surfaces().perpendicularAngle();

    meshRefiner_.baffleAndSplitMesh
    (
        handleSnapProblems,

        // Snap problem cell detection
        snapParams,
        refineParams.useTopologicalSnapDetection(),
        handleSnapProblems,                 // remove perp edge connected cells
        perpAngle,                          // perp angle
        refineParams.nErodeCellZone(),

        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.locationsInMesh(),
        refineParams.zonesInMesh(),
        refineParams.locationsOutsideMesh()
    );

    // Merge free-standing baffles always
    meshRefiner_.mergeFreeStandingBaffles
    (
        snapParams,
        refineParams.useTopologicalSnapDetection(),
        handleSnapProblems,
        perpAngle,
        refineParams.planarAngle(),
        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.locationsInMesh(),
        refineParams.locationsOutsideMesh()
    );

    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    // Duplicate points on baffles that are on more than one cell
    // region. This will help snapping pull them to separate surfaces.
    meshRefiner_.dupNonManifoldPoints();


    // Merge all baffles that are still remaining after duplicating points.
    List<labelPair> couples(localPointRegion::findDuplicateFacePairs(mesh));

    label nCouples = returnReduce(couples.size(), sumOp<label>());

    Info<< "Detected unsplittable baffles : " << nCouples << endl;

    if (nCouples > 0)
    {
        // Actually merge baffles. Note: not exactly parallellized. Should
        // convert baffle faces into processor faces if they resulted
        // from them.
        meshRefiner_.mergeBaffles(couples, Map<label>(0));

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            meshRefiner_.checkData();
        }

        // Remove any now dangling parts
        meshRefiner_.splitMeshRegions
        (
            globalToMasterPatch_,
            globalToSlavePatch_,
            refineParams.locationsInMesh(),
            refineParams.locationsOutsideMesh()
        );

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            meshRefiner_.checkData();
        }

        Info<< "Merged free-standing baffles in = "
            << mesh.time().cpuTimeIncrement() << " s." << endl;
    }

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing handleProblemCells mesh to time "
            << meshRefiner_.timeName() << endl;
        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner_.timeName()
        );
    }
}


void Foam::snappyRefineDriver::addFaceZones
(
    meshRefinement& meshRefiner,
    const refinementParameters& refineParams,
    const HashTable<Pair<word>>& faceZoneToPatches
)
{
    if (faceZoneToPatches.size())
    {
        Info<< nl
            << "Adding patches for face zones" << nl
            << "-----------------------------" << nl
            << endl;

        Info<< setf(ios_base::left)
            << setw(6) << "Patch"
            << setw(20) << "Type"
            << setw(30) << "Name"
            << setw(30) << "FaceZone"
            << setw(10) << "FaceType"
            << nl
            << setw(6) << "-----"
            << setw(20) << "----"
            << setw(30) << "----"
            << setw(30) << "--------"
            << setw(10) << "--------"
            << endl;

        const polyMesh& mesh = meshRefiner.mesh();

        // Add patches for added inter-region faceZones
        forAllConstIter(HashTable<Pair<word>>, faceZoneToPatches, iter)
        {
            const word& fzName = iter.key();
            const Pair<word>& patchNames = iter();

            // Get any user-defined faceZone data
            surfaceZonesInfo::faceZoneType fzType;
            dictionary patchInfo = refineParams.getZoneInfo(fzName, fzType);

            const word& masterName = fzName;
            //const word slaveName = fzName + "_slave";
            //const word slaveName = czNames.second()+"_to_"+czNames.first();
            const word& slaveName = patchNames.second();

            label mpi = meshRefiner.addMeshedPatch(masterName, patchInfo);

            Info<< setf(ios_base::left)
                << setw(6) << mpi
                << setw(20) << mesh.boundaryMesh()[mpi].type()
                << setw(30) << masterName
                << setw(30) << fzName
                << setw(10) << surfaceZonesInfo::faceZoneTypeNames[fzType]
                << nl;


            label sli = meshRefiner.addMeshedPatch(slaveName, patchInfo);

            Info<< setf(ios_base::left)
                << setw(6) << sli
                << setw(20) << mesh.boundaryMesh()[sli].type()
                << setw(30) << slaveName
                << setw(30) << fzName
                << setw(10) << surfaceZonesInfo::faceZoneTypeNames[fzType]
                << nl;

            meshRefiner.addFaceZone(fzName, masterName, slaveName, fzType);
        }

        Info<< endl;
    }
}


void Foam::snappyRefineDriver::mergePatchFaces
(
    const bool geometricMerge,
    const refinementParameters& refineParams,
    const dictionary& motionDict
)
{
    addProfiling(merge, "snappyHexMesh::refine::merge");
    Info<< nl
        << "Merge refined boundary faces" << nl
        << "----------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (geometricMerge)
    {
        meshRefiner_.mergePatchFacesUndo
        (
            Foam::cos(degToRad(45.0)),
            Foam::cos(degToRad(45.0)),
            meshRefiner_.meshedPatches(),
            motionDict,
            labelList(mesh.nFaces(), -1)
        );
    }
    else
    {
        // Still merge refined boundary faces if all four are on same patch
        meshRefiner_.mergePatchFaces
        (
            Foam::cos(degToRad(45.0)),
            Foam::cos(degToRad(45.0)),
            4,          // only merge faces split into 4
            meshRefiner_.meshedPatches()
        );
    }

    if (debug)
    {
        meshRefiner_.checkData();
    }

    meshRefiner_.mergeEdgesUndo(Foam::cos(degToRad(45.0)), motionDict);

    if (debug)
    {
        meshRefiner_.checkData();
    }
}


void Foam::snappyRefineDriver::doRefine
(
    const dictionary& refineDict,
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool prepareForSnapping,
    const bool doMergePatchFaces,
    const dictionary& motionDict
)
{
    addProfiling(refine, "snappyHexMesh::refine");
    Info<< nl
        << "Refinement phase" << nl
        << "----------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    // Check that all the keep points are inside the mesh.
    refineParams.findCells(true, mesh, refineParams.locationsInMesh());

    // Refine around feature edges
    featureEdgeRefine
    (
        refineParams,
        100,    // maxIter
        0       // min cells to refine
    );


    if
    (
        max(meshRefiner_.surfaces().maxGapLevel()) > 0
     || max(meshRefiner_.shells().maxGapLevel()) > 0
    )
    {
        // In case we use automatic gap level refinement do some pre-refinement
        // (fast) since it is so slow.

        // Refine based on surface
        surfaceOnlyRefine
        (
            refineParams,
            20     // maxIter
        );

        // Refine cells that contain a gap
        smallFeatureRefine
        (
            refineParams,
            100     // maxIter
        );
    }


    // Refine based on surface
    surfaceOnlyRefine
    (
        refineParams,
        100     // maxIter
    );

    gapOnlyRefine
    (
        refineParams,
        100     // maxIter
    );

    // Remove cells (a certain distance) beyond surface intersections
    removeInsideCells
    (
        refineParams,
        1       // nBufferLayers
    );

    // Refine consistently across narrow gaps (a form of shell refinement)
    bigGapOnlyRefine
    (
        refineParams,
        true,   // spreadGapSize
        100     // maxIter
    );

    // Internal mesh refinement
    shellRefine
    (
        refineParams,
        100    // maxIter
    );

    // Refine any hexes with 5 or 6 faces refined to make smooth edges
    danglingCellRefine
    (
        refineParams,
        21,     // 1 coarse face + 5 refined faces
        100     // maxIter
    );
    danglingCellRefine
    (
        refineParams,
        24,     // 0 coarse faces + 6 refined faces
        100     // maxIter
    );

    // Refine any cells with differing refinement level on either side
    refinementInterfaceRefine
    (
        refineParams,
        10      // maxIter
    );

    // Introduce baffles at surface intersections. Remove sections unreachable
    // from keepPoint.
    baffleAndSplitMesh
    (
        refineParams,
        snapParams,
        prepareForSnapping,
        motionDict
    );

    // Mesh is at its finest. Do optional zoning (cellZones and faceZones)
    wordPairHashTable zonesToFaceZone;
    zonify(refineParams, zonesToFaceZone);

    // Create pairs of patches for faceZones
    {
        HashTable<Pair<word>> faceZoneToPatches(zonesToFaceZone.size());

        //    Note: zonesToFaceZone contains the same data on different
        //          processors but in different order. We could sort the
        //          contents but instead just loop in sortedToc order.
        List<Pair<word>> czs(zonesToFaceZone.sortedToc());

        forAll(czs, i)
        {
            const Pair<word>& czNames = czs[i];
            const word& fzName = zonesToFaceZone[czNames];

            const word& masterName = fzName;
            const word slaveName = czNames.second() + "_to_" + czNames.first();
            Pair<word> patches(masterName, slaveName);
            faceZoneToPatches.insert(fzName, patches);
        }
        addFaceZones(meshRefiner_, refineParams, faceZoneToPatches);
    }

    // Pull baffles apart
    splitAndMergeBaffles
    (
        refineParams,
        snapParams,
        prepareForSnapping,
        motionDict
    );

    // Do something about cells with refined faces on the boundary
    if (prepareForSnapping)
    {
        mergePatchFaces(doMergePatchFaces, refineParams, motionDict);
    }


    if (Pstream::parRun())
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;

        // Do final balancing. Keep zoned faces on one processor since the
        // snap phase will convert them to baffles and this only works for
        // internal faces.
        meshRefiner_.balance
        (
            true,                           // keepZoneFaces
            false,                          // keepBaffles
            scalarField(mesh.nCells(), 1),  // cellWeights
            decomposer_,
            distributor_
        );
    }
}


// ************************************************************************* //
