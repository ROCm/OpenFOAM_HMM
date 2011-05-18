/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "backgroundMeshDecomposition.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(backgroundMeshDecomposition, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::backgroundMeshDecomposition::initialRefinement()
{
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

    scalar mergeDist = 1e-6*mesh_.bounds().mag();

    volScalarField cellWeights
    (
        IOobject
        (
            "cellWeights",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
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
            "void Foam::backgroundMeshDecomposition::initialRefinement() const"
        )
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which is not parallel aware." << endl
            << exit(FatalError);
    }

    hexRef8 meshCutter
    (
        mesh_,
        labelList(mesh_.nCells(), 0),
        labelList(mesh_.nPoints(), 0)
    );

    // For each cell in the mesh has it been determined if it is fully
    // inside, outside, or overlaps the surface
    labelList volumeStatus
    (
        mesh_.nCells(),
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
                        mesh_.cells()[cellI].points
                        (
                            mesh_.faces(),
                            mesh_.points()
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

            {
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


                if
                (
                    returnReduce(newCellsToRefine.size(), sumOp<label>())
                 == 0
                )
                {
                    break;
                }

                // Mesh changing engine.
                polyTopoChange meshMod(mesh_);

                // Play refinement commands into mesh changer.
                meshCutter.setRefinement(newCellsToRefine, meshMod);

                // Create mesh, return map from old to new mesh.
                autoPtr<mapPolyMesh> map = meshMod.changeMesh
                (
                    mesh_,
                    false,  // inflate
                    true,   // syncParallel
                    true,   // orderCells (to reduce cell motion in scotch)
                    false   // orderPoints
                );

                // Update fields
                mesh_.updateMesh(map);

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
                            newVolumeStatus[newCellI] =
                            volumeStatus[oldCellI];
                        }
                    }

                    volumeStatus.transfer(newVolumeStatus);
                }

                if (debug)
                {
                    Info<< "Refined from "
                        << returnReduce(map().nOldCells(), sumOp<label>())
                        << " to " << mesh_.globalData().nTotalCells()
                        << " cells." << endl;
                }
            }

            // Determine/update the status of each cell
            forAll(volumeStatus, cellI)
            {
                if (volumeStatus[cellI] == searchableSurface::UNKNOWN)
                {
                    treeBoundBox cellBb
                    (
                        mesh_.cells()[cellI].points
                        (
                            mesh_.faces(),
                            mesh_.points()
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

            bool removeOutsideCells = false;

            if (removeOutsideCells)
            {
                DynamicList<label> cellsToRemove;

                forAll(volumeStatus, cellI)
                {
                    if (volumeStatus[cellI] == searchableSurface::OUTSIDE)
                    {
                        cellsToRemove.append(cellI);
                    }
                }

                removeCells cellRemover(mesh_);

                // Mesh changing engine.
                polyTopoChange meshMod(mesh_);

                labelList exposedFaces = cellRemover.getExposedFaces
                (
                    cellsToRemove
                );

                // Play refinement commands into mesh changer.
                cellRemover.setRefinement
                (
                    cellsToRemove,
                    exposedFaces,
                    labelList(exposedFaces.size(), 0),  // patchID dummy
                    meshMod
                );

                // Create mesh, return map from old to new mesh.
                autoPtr<mapPolyMesh> map = meshMod.changeMesh
                (
                    mesh_,
                    false,  // inflate
                    true,   // syncParallel
                    true,   // orderCells (to reduce cell motion in scotch)
                    false   // orderPoints
                );

                // Update fields
                mesh_.updateMesh(map);

                // Update numbering of cells/vertices.
                meshCutter.updateMesh(map);
                cellRemover.updateMesh(map);

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
                            newVolumeStatus[newCellI] =
                                volumeStatus[oldCellI];
                        }
                    }

                    volumeStatus.transfer(newVolumeStatus);
                }

                Info<< "Removed "
                    << returnReduce(map().nOldCells(), sumOp<label>())
                     - mesh_.globalData().nTotalCells()
                    << " cells." << endl;
            }

            if (debug)
            {
                const_cast<Time&>(mesh_.time())++;
                meshCutter.write();
                mesh_.write();
                cellWeights.write();
            }

            labelList newDecomp = decomposer.decompose
            (
                mesh_,
                mesh_.cellCentres(),
                cellWeights
            );

            fvMeshDistribute distributor(mesh_, mergeDist);

            autoPtr<mapDistributePolyMesh> mapDist =
                distributor.distribute(newDecomp);

            meshCutter.distribute(mapDist);

            mapDist().distributeCellData(volumeStatus);

            if (debug)
            {
                printMeshData(mesh_);

                const_cast<Time&>(mesh_.time())++;
                meshCutter.write();
                mesh_.write();
                cellWeights.write();
            }
        }
    }

    if (debug)
    {
        const_cast<Time&>(mesh_.time())++;
        cellWeights.write();
        mesh_.write();
    }
}


void Foam::backgroundMeshDecomposition::printMeshData
(
    const polyMesh& mesh
) const
{
    // Collect all data on master

    globalIndex globalCells(mesh.nCells());
    // labelListList patchNeiProcNo(Pstream::nProcs());
    // labelListList patchSize(Pstream::nProcs());
    // const labelList& pPatches = mesh.globalData().processorPatches();
    // patchNeiProcNo[Pstream::myProcNo()].setSize(pPatches.size());
    // patchSize[Pstream::myProcNo()].setSize(pPatches.size());
    // forAll(pPatches, i)
    // {
    //     const processorPolyPatch& ppp = refCast<const processorPolyPatch>
    //     (
    //         mesh.boundaryMesh()[pPatches[i]]
    //     );
    //     patchNeiProcNo[Pstream::myProcNo()][i] = ppp.neighbProcNo();
    //     patchSize[Pstream::myProcNo()][i] = ppp.size();
    // }
    // Pstream::gatherList(patchNeiProcNo);
    // Pstream::gatherList(patchSize);


    // // Print stats

    // globalIndex globalBoundaryFaces(mesh.nFaces()-mesh.nInternalFaces());

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        Info<< "Processor " << procI << " "
            << "Number of cells = " << globalCells.localSize(procI)
            << endl;

        // label nProcFaces = 0;

        // const labelList& nei = patchNeiProcNo[procI];

        // forAll(patchNeiProcNo[procI], i)
        // {
        //     Info<< "    Number of faces shared with processor "
        //         << patchNeiProcNo[procI][i] << " = " << patchSize[procI][i]
        //         << endl;

        //     nProcFaces += patchSize[procI][i];
        // }

        // Info<< "    Number of processor patches = " << nei.size() << nl
        //     << "    Number of processor faces = " << nProcFaces << nl
        //     << "    Number of boundary faces = "
        //     << globalBoundaryFaces.localSize(procI) << endl;
    }
}


bool Foam::backgroundMeshDecomposition::refineCell
(
    const polyMesh& mesh,
    label cellI,
    label volType,
    scalar& weightEstimate
) const
{
    // Sample the box to find an estimate of the min size, and a volume
    // estimate when overlapping == true.

    const conformationSurfaces& geometry = cvMesh_.geometryToConformTo();

    treeBoundBox cellBb
    (
        mesh.cells()[cellI].points
        (
            mesh.faces(),
            mesh.points()
        )
    );

    weightEstimate = 1.0;

    if (volType == searchableSurface::MIXED)
    {
        // Assess the cell size at the nearest point on the surface for the
        // MIXED cells, if the cell is large with respect to the cell size,
        // then refine it.

        pointField samplePoints
        (
            volRes_*volRes_*volRes_,
            vector::zero
        );

        // scalar sampleVol = cellBb.volume()/samplePoints.size();

        vector delta = cellBb.span()/volRes_;

        label pI = 0;

        for (label i = 0; i < volRes_; i++)
        {
            for (label j = 0; j < volRes_; j++)
            {
                for (label k = 0; k < volRes_; k++)
                {
                    samplePoints[pI++] =
                        cellBb.min()
                      + vector
                        (
                            delta.x()*(i + 0.5),
                            delta.y()*(j + 0.5),
                            delta.z()*(k + 0.5)
                        );
                }
            }
        }

        List<pointIndexHit> hitInfo;
        labelList hitSurfaces;

        geometry.findSurfaceNearest
        (
            samplePoints,
            scalarField(samplePoints.size(), sqr(GREAT)),
            hitInfo,
            hitSurfaces
        );

        // weightEstimate = 0.0;

        scalar minCellSize = GREAT;

        forAll(samplePoints, i)
        {
            scalar s = cvMesh_.cellSizeControl().cellSize
            (
                hitInfo[i].hitPoint()
            );

            // Info<< "cellBb.midpoint() " << cellBb.midpoint() << nl
            //     << samplePoints[i] << nl
            //     << hitInfo[i] << nl
            //     << "cellBb.span() " << cellBb.span() << nl
            //     << "cellBb.mag() " << cellBb.mag() << nl
            //     << s << endl;

            if (s < minCellSize)
            {
                minCellSize = max(s, minCellSizeLimit_);
            }

            // Estimate the number of points in the cell by the surface size,
            // this is likely to be too small, so reduce.
            // weightEstimate += sampleVol/pow3(s);
        }

        if (sqr(spanScale_)*sqr(minCellSize) < magSqr(cellBb.span()))
        {
            return true;
        }
    }
    else if (volType == searchableSurface::INSIDE)
    {
        // scalar s = cvMesh_.cellSizeControl().cellSize(cellBb.midpoint());

        // Estimate the number of points in the cell by the size at the cell
        // midpoint
        // weightEstimate = cellBb.volume()/pow3(s);

        return false;
    }
    // else
    // {
    //     weightEstimate = 1.0;

    //     return false;
    // }

    return false;
}


Foam::labelList Foam::backgroundMeshDecomposition::selectRefinementCells
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
            if (meshCutter.cellLevel()[cellI] < minLevels_)
            {
                cellsToRefine.insert(cellI);
            }
        }

        if (volumeStatus[cellI] != searchableSurface::OUTSIDE)
        {
            if
            (
                refineCell
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
    }

    return cellsToRefine.toc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backgroundMeshDecomposition::backgroundMeshDecomposition
(
    const dictionary& coeffsDict,
    const conformalVoronoiMesh& cvMesh
)
:
    coeffsDict_(coeffsDict),
    cvMesh_(cvMesh),
    mesh_
    (
        IOobject
        (
            fvMesh::defaultRegion,
            cvMesh_.time().timeName(),
            cvMesh_.time(),
            IOobject::MUST_READ
        )
    ),
    spanScale_(readScalar(coeffsDict_.lookup("spanScale"))),
    minCellSizeLimit_
    (
        coeffsDict_.lookupOrDefault<scalar>("minCellSizeLimit", 0.0)
    ),
    minLevels_(readLabel(coeffsDict_.lookup("minLevels"))),
    volRes_(readLabel(coeffsDict_.lookup("sampleResolution")))
{
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "Foam::backgroundMeshDecomposition::backgroundMeshDecomposition"
            "("
                "const dictionary& coeffsDict, "
                "const conformalVoronoiMesh& cvMesh"
            ")"
        )
            << "This cannot be used when not running in parallel."
            << exit(FatalError);
    }

    Info<< nl << "Building initial background mesh decomposition" << endl;

    initialRefinement();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::backgroundMeshDecomposition::~backgroundMeshDecomposition()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::conformationSurfaces::positionOnThisProc(const point& pt) const
{

}


Foam::label Foam::conformationSurfaces::positionProc(const point& pt) const
{


    return -1;
}


// ************************************************************************* //
