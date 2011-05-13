/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "parallelHierarchicalDensityWeightedStochastic.H"
#include "addToRunTimeSelectionTable.H"

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(parallelHierarchicalDensityWeightedStochastic, 0);
addToRunTimeSelectionTable
(
    initialPointsMethod,
    parallelHierarchicalDensityWeightedStochastic,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void parallelHierarchicalDensityWeightedStochastic::printMeshData
(
    const polyMesh& mesh
) const
{
    // Collect all data on master

    globalIndex globalCells(mesh.nCells());
    labelListList patchNeiProcNo(Pstream::nProcs());
    labelListList patchSize(Pstream::nProcs());
    const labelList& pPatches = mesh.globalData().processorPatches();
    patchNeiProcNo[Pstream::myProcNo()].setSize(pPatches.size());
    patchSize[Pstream::myProcNo()].setSize(pPatches.size());
    forAll(pPatches, i)
    {
        const processorPolyPatch& ppp = refCast<const processorPolyPatch>
        (
            mesh.boundaryMesh()[pPatches[i]]
        );
        patchNeiProcNo[Pstream::myProcNo()][i] = ppp.neighbProcNo();
        patchSize[Pstream::myProcNo()][i] = ppp.size();
    }
    Pstream::gatherList(patchNeiProcNo);
    Pstream::gatherList(patchSize);


    // Print stats

    globalIndex globalBoundaryFaces(mesh.nFaces()-mesh.nInternalFaces());

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        Info<< endl
            << "Processor " << procI << nl
            << "    Number of cells = " << globalCells.localSize(procI)
            << endl;

        label nProcFaces = 0;

        const labelList& nei = patchNeiProcNo[procI];

        forAll(patchNeiProcNo[procI], i)
        {
            Info<< "    Number of faces shared with processor "
                << patchNeiProcNo[procI][i] << " = " << patchSize[procI][i]
                << endl;

            nProcFaces += patchSize[procI][i];
        }

        Info<< "    Number of processor patches = " << nei.size() << nl
            << "    Number of processor faces = " << nProcFaces << nl
            << "    Number of boundary faces = "
            << globalBoundaryFaces.localSize(procI) << endl;
    }
}


bool parallelHierarchicalDensityWeightedStochastic::estimateCellWeight
(
    const polyMesh& mesh,
    label cellI,
    label volType,
    scalar& weightEstimate
) const
{
    weightEstimate = 10*scalar(Pstream::myProcNo() + 1.0);

    return false;
}


labelList parallelHierarchicalDensityWeightedStochastic::selectRefinementCells
(
    const hexRef8& meshCutter,
    labelList& volumeStatus,
    volScalarField& cellWeights
) const
{
    labelHashSet cellsToRefine;

    const polyMesh& mesh = meshCutter.mesh();

    // Determine/update the status of each cell
    forAll(volumeStatus, cellI)
    {
        if (volumeStatus[cellI] == searchableSurface::MIXED)
        {
            if (meshCutter.cellLevel()[cellI] <= minLevels_)
            {
                cellsToRefine.insert(cellI);
            }
        }

        if
        (
            estimateCellWeight
            (
                mesh,
                cellI,
                volumeStatus[cellI],
                cellWeights.internalField()[cellI]
            )
        )
        {
            cellsToRefine.insert(cellI);
        }
    }

    return cellsToRefine.toc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parallelHierarchicalDensityWeightedStochastic::
parallelHierarchicalDensityWeightedStochastic
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh
)
:
    initialPointsMethod(typeName, initialPointsDict, cvMesh),
    globalTrialPoints_(0),
    minCellSizeLimit_
    (
        detailsDict().lookupOrDefault<scalar>("minCellSizeLimit", 0.0)
    ),
    minLevels_(readLabel(detailsDict().lookup("minLevels"))),
    maxSizeRatio_(readScalar(detailsDict().lookup("maxSizeRatio"))),
    volRes_(readLabel(detailsDict().lookup("sampleResolution"))),
    surfRes_
    (
        detailsDict().lookupOrDefault<label>("surfaceSampleResolution", volRes_)
    )
{
    if (maxSizeRatio_ <= 1.0)
    {
        maxSizeRatio_ = 2.0;

        WarningIn
        (
            "parallelHierarchicalDensityWeightedStochastic::"
            "parallelHierarchicalDensityWeightedStochastic"
            "("
                "const dictionary& initialPointsDict,"
                "const conformalVoronoiMesh& cvMesh"
            ")"
        )
            << "The maxSizeRatio must be greater than one to be sensible, "
            << "setting to " << maxSizeRatio_
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::vector<Vb::Point>
parallelHierarchicalDensityWeightedStochastic::initialPoints() const
{
    std::vector<Vb::Point> initialPoints;

    if (Pstream::parRun())
    {
        fvMesh mesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                cvMesh_.time().timeName(),
                cvMesh_.time(),
                IOobject::MUST_READ
            )
        );

        //Read decomposePar dictionary
        IOdictionary decomposeDict
        (
            IOobject
            (
                "decomposeParDict",
                cvMesh_.time().system(),
                cvMesh_.time(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        scalar mergeDist = 1e-6*mesh.bounds().mag();

        volScalarField cellWeights
        (
            IOobject
            (
                "cellWeights",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("one", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        );

        const conformationSurfaces& geometry = cvMesh_.geometryToConformTo();

        // Decomposition
        autoPtr<decompositionMethod> decomposerPtr
        (
            decompositionMethod::New(decomposeDict)
        );

        decompositionMethod& decomposer = decomposerPtr();

        if (!decomposer.parallelAware())
        {
            FatalErrorIn
            (
                "parallelHierarchicalDensityWeightedStochastic"
            )
                << "You have selected decomposition method "
                << decomposer.typeName
                << " which is not parallel aware." << endl
                << exit(FatalError);
        }

        hexRef8 meshCutter
        (
            mesh,
            labelList(mesh.nCells(), 0),
            labelList(mesh.nPoints(), 0)
        );

        // For each cell in the mesh has it been determined if it is fully
        // inside, outside, or overlaps the surface
        labelList volumeStatus
        (
            mesh.nCells(),
            searchableSurface::UNKNOWN
        );

        // Surface refinement
        {
            while (true)
            {
                // Determine/update the status of each cell
                forAll(volumeStatus, cellI)
                {
                    if (volumeStatus[cellI] == searchableSurface::UNKNOWN)
                    {
                        treeBoundBox cellBb
                        (
                            mesh.cells()[cellI].points
                            (
                                mesh.faces(),
                                mesh.points()
                            )
                        );

                        if (geometry.overlaps(cellBb))
                        {
                            volumeStatus[cellI] = searchableSurface::MIXED;
                        }
                        else if (geometry.inside(cellBb.midpoint()))
                        {
                            volumeStatus[cellI] = searchableSurface::INSIDE;
                        }
                        else
                        {
                            volumeStatus[cellI] =
                            searchableSurface::OUTSIDE;
                        }
                    }
                }

                labelList refCells = selectRefinementCells
                (
                    meshCutter,
                    volumeStatus,
                    cellWeights
                );

                // Maintain 2:1 ratio
                labelList newCellsToRefine
                (
                    meshCutter.consistentRefinement
                    (
                        refCells,
                        true                  // extend set
                    )
                );

                forAll(newCellsToRefine, nCTRI)
                {
                    label cellI = newCellsToRefine[nCTRI];

                    if (volumeStatus[cellI] == searchableSurface::MIXED)
                    {
                        volumeStatus[cellI] = searchableSurface::UNKNOWN;
                    }

                    cellWeights.internalField()[cellI] = max
                    (
                        1.0,
                        cellWeights.internalField()[cellI]/8.0
                    );
                }

                if (newCellsToRefine.empty())
                {
                    break;
                }

                // Mesh changing engine.
                polyTopoChange meshMod(mesh);

                // Play refinement commands into mesh changer.
                meshCutter.setRefinement(newCellsToRefine, meshMod);

                // Create mesh, return map from old to new mesh.
                autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

                // Update fields
                mesh.updateMesh(map);

                // Update numbering of cells/vertices.
                meshCutter.updateMesh(map);

                {
                    // Map volumeStatus

                    const labelList& cellMap = map().cellMap();

                    labelList newVolumeStatus(cellMap.size());

                    forAll(cellMap, newCellI)
                    {
                        label oldCellI = cellMap[newCellI];

                        if (oldCellI == -1)
                        {
                            newVolumeStatus[newCellI] =
                                searchableSurface::UNKNOWN;;
                        }
                        else
                        {
                            newVolumeStatus[newCellI] = volumeStatus[oldCellI];
                        }
                    }

                    volumeStatus.transfer(newVolumeStatus);
                }

                Info<< "Refined from "
                    << returnReduce(map().nOldCells(), sumOp<label>())
                    << " to " << mesh.globalData().nTotalCells() << " cells."
                    << nl << endl;

                const_cast<Time&>(mesh.time())++;
                meshCutter.write();
                mesh.write();
                cellWeights.write();

                labelList newDecomp = decomposer.decompose
                (
                    mesh,
                    mesh.cellCentres(),
                    cellWeights
                );

                fvMeshDistribute distributor(mesh, mergeDist);

                autoPtr<mapDistributePolyMesh> mapDist =
                    distributor.distribute(newDecomp);

                meshCutter.distribute(mapDist);

                mapDist().distributeCellData(volumeStatus);

                printMeshData(mesh);

                const_cast<Time&>(mesh.time())++;
                meshCutter.write();
                mesh.write();
                cellWeights.write();
            }
        }

        // const_cast<Time&>(mesh.time())++;
        // cellWeights.write();
        // mesh.write();

        // Pout<< "crash now" << endl;
        // Pout<< acos(2.0) << endl;
    }

    return initialPoints;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
