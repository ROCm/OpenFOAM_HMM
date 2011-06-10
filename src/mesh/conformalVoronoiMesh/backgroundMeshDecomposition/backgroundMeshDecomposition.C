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
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(backgroundMeshDecomposition, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::backgroundMeshDecomposition::initialRefinement()
{
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

    decompositionMethod& decomposer = decomposerPtr_();

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
                    volumeStatus,
                    cellWeights
                );

                // Maintain 2:1 ratio
                labelList newCellsToRefine
                (
                    meshCutter_.consistentRefinement
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

                if (returnReduce(newCellsToRefine.size(), sumOp<label>()) == 0)
                {
                    break;
                }

                // Mesh changing engine.
                polyTopoChange meshMod(mesh_);

                // Play refinement commands into mesh changer.
                meshCutter_.setRefinement(newCellsToRefine, meshMod);

                // Create mesh, return map from old to new mesh.
                autoPtr<mapPolyMesh> map = meshMod.changeMesh
                (
                    mesh_,
                    false,  // inflate
                    true,   // syncParallel
                    true,   // orderCells (to reduce cell transfers)
                    false   // orderPoints
                );

                // Update fields
                mesh_.updateMesh(map);

                // Update numbering of cells/vertices.
                meshCutter_.updateMesh(map);

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
                                searchableSurface::UNKNOWN;
                        }
                        else
                        {
                            newVolumeStatus[newCellI] = volumeStatus[oldCellI];
                        }
                    }

                    volumeStatus.transfer(newVolumeStatus);
                }

                if (debug)
                {
                    Info<< "    Background mesh refined from "
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

            // Hard code switch for this stage for testing
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
                    true,   // orderCells (to reduce cell transfers)
                    false   // orderPoints
                );

                // Update fields
                mesh_.updateMesh(map);

                // Update numbering of cells/vertices.
                meshCutter_.updateMesh(map);
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
                Info<< "Time " << mesh_.time().timeName() << endl;
                meshCutter_.write();
                mesh_.write();
                cellWeights.write();
            }

            labelList newDecomp = decomposer.decompose
            (
                mesh_,
                mesh_.cellCentres(),
                cellWeights
            );

            fvMeshDistribute distributor(mesh_, mergeDist_);

            autoPtr<mapDistributePolyMesh> mapDist =
                distributor.distribute(newDecomp);

            meshCutter_.distribute(mapDist);

            mapDist().distributeCellData(volumeStatus);

            if (debug)
            {
                printMeshData(mesh_);

                const_cast<Time&>(mesh_.time())++;
                Info<< "Time " << mesh_.time().timeName() << endl;
                meshCutter_.write();
                mesh_.write();
                cellWeights.write();
            }
        }
    }

    if (debug)
    {
        const_cast<Time&>(mesh_.time())++;
        Info<< "Time " << mesh_.time().timeName() << endl;
        cellWeights.write();
        mesh_.write();
    }

    buildPatchAndTree();
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
        mesh_.cells()[cellI].points
        (
            mesh_.faces(),
            mesh_.points()
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
    labelList& volumeStatus,
    volScalarField& cellWeights
) const
{
    labelHashSet cellsToRefine;

    // Determine/update the status of each cell
    forAll(volumeStatus, cellI)
    {
        if (volumeStatus[cellI] == searchableSurface::MIXED)
        {
            if (meshCutter_.cellLevel()[cellI] < minLevels_)
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


void Foam::backgroundMeshDecomposition::buildPatchAndTree()
{
    primitivePatch tmpBoundaryFaces
    (
        SubList<face>
        (
            mesh_.faces(),
            mesh_.nFaces() - mesh_.nInternalFaces(),
            mesh_.nInternalFaces()
        ),
        mesh_.points()
    );

    boundaryFacesPtr_.reset
    (
        new bPatch
        (
            tmpBoundaryFaces.localFaces(),
            tmpBoundaryFaces.localPoints()
        )
    );

    // Overall bb
    treeBoundBox overallBb(boundaryFacesPtr_().localPoints());

    Random& rnd = cvMesh_.rndGen();

    bFTreePtr_.reset
    (
        new indexedOctree<treeDataBPatch>
        (
            treeDataBPatch(false, boundaryFacesPtr_()),
            overallBb.extend(rnd, 1e-4),
            10, // maxLevel
            10, // leafSize
            3.0 // duplicity
        )
    );

    // Give the bounds of every processor to every other processor
    allBackgroundMeshBounds_[Pstream::myProcNo()] = overallBb;

    Pstream::gatherList(allBackgroundMeshBounds_);
    Pstream::scatterList(allBackgroundMeshBounds_);

    if (debug)
    {
        OFstream fStr
        (
            cvMesh_.time().path()
           /"backgroundMeshDecomposition_proc_"
          + name(Pstream::myProcNo())
          + "_boundaryFaces.obj"
        );

        const faceList& faces = boundaryFacesPtr_().localFaces();
        const List<point>& points = boundaryFacesPtr_().localPoints();

        Map<label> foamToObj(points.size());

        label vertI = 0;

        forAll(faces, i)
        {
            const face& f = faces[i];

            forAll(f, fPI)
            {
                if (foamToObj.insert(f[fPI], vertI))
                {
                    meshTools::writeOBJ(fStr, points[f[fPI]]);
                    vertI++;
                }
            }

            fStr<< 'f';

            forAll(f, fPI)
            {
                fStr<< ' ' << foamToObj[f[fPI]] + 1;
            }

            fStr<< nl;
        }
    }
}


Foam::autoPtr<Foam::mapDistribute> Foam::backgroundMeshDecomposition::buildMap
(
    const List<label>& toProc
) const
{
    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // 1. Count
    labelList nSend(Pstream::nProcs(), 0);

    forAll(toProc, i)
    {
        label procI = toProc[i];

        nSend[procI]++;
    }

    // Send over how many I need to receive
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList sendSizes(Pstream::nProcs());

    sendSizes[Pstream::myProcNo()] = nSend;

    combineReduce(sendSizes, UPstream::listEq());

    // 2. Size sendMap
    labelListList sendMap(Pstream::nProcs());

    forAll(nSend, procI)
    {
        sendMap[procI].setSize(nSend[procI]);

        nSend[procI] = 0;
    }

    // 3. Fill sendMap
    forAll(toProc, i)
    {
        label procI = toProc[i];

        sendMap[procI][nSend[procI]++] = i;
    }

    // Determine receive map
    // ~~~~~~~~~~~~~~~~~~~~~

    labelListList constructMap(Pstream::nProcs());

    // Local transfers first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label constructSize = constructMap[Pstream::myProcNo()].size();

    forAll(constructMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            label nRecv = sendSizes[procI][Pstream::myProcNo()];

            constructMap[procI].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[procI][i] = constructSize++;
            }
        }
    }

    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            constructSize,
            sendMap.xfer(),
            constructMap.xfer()
        )
    );
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
    meshCutter_
    (
        mesh_,
        labelList(mesh_.nCells(), 0),
        labelList(mesh_.nPoints(), 0)
    ),
    boundaryFacesPtr_(),
    bFTreePtr_(),
    allBackgroundMeshBounds_(Pstream::nProcs()),
    decomposeDict_
    (
        IOobject
        (
            "decomposeParDict",
            cvMesh_.time().system(),
            cvMesh_.time(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    decomposerPtr_(decompositionMethod::New(decomposeDict_)),
    mergeDist_(1e-6*mesh_.bounds().mag()),
    spanScale_(readScalar(coeffsDict_.lookup("spanScale"))),
    minCellSizeLimit_
    (
        coeffsDict_.lookupOrDefault<scalar>("minCellSizeLimit", 0.0)
    ),
    minLevels_(readLabel(coeffsDict_.lookup("minLevels"))),
    volRes_(readLabel(coeffsDict_.lookup("sampleResolution"))),
    maxCellWeightCoeff_(readScalar(coeffsDict_.lookup("maxCellWeightCoeff")))
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

    if (!decomposerPtr_().parallelAware())
    {
        FatalErrorIn
        (
            "void Foam::backgroundMeshDecomposition::initialRefinement() const"
        )
            << "You have selected decomposition method "
            << decomposerPtr_().typeName
            << " which is not parallel aware." << endl
            << exit(FatalError);
    }

    Info<< nl << "Building initial background mesh decomposition" << endl;

    initialRefinement();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::backgroundMeshDecomposition::~backgroundMeshDecomposition()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::backgroundMeshDecomposition::distribute
(
    volScalarField& cellWeights,
    List<DynamicList<point> >& cellVertices
)
{
    if (debug)
    {
        const_cast<Time&>(mesh_.time())++;
        Info<< "Time " << mesh_.time().timeName() << endl;
        cellWeights.write();
        mesh_.write();
    }

    while (true)
    {
        // Refine large cells if necessary

        scalar cellWeightLimit =
            maxCellWeightCoeff_
           *sum(cellWeights).value()
           /mesh_.globalData().nTotalCells();

        if (debug)
        {
            Info<< "    cellWeightLimit " << cellWeightLimit << endl;

            Pout<< "    sum(cellWeights) "
                << sum(cellWeights.internalField())
                << " max(cellWeights) "
                << max(cellWeights.internalField())
                << endl;
        }

        labelHashSet cellsToRefine;

        forAll(cellWeights, cWI)
        {
            if (cellWeights.internalField()[cWI] > cellWeightLimit)
            {
                cellsToRefine.insert(cWI);
            }
        }

        if (returnReduce(cellsToRefine.size(), sumOp<label>()) == 0)
        {
            break;
        }

        // Maintain 2:1 ratio
        labelList newCellsToRefine
        (
            meshCutter_.consistentRefinement
            (
                cellsToRefine.toc(),
                true                  // extend set
            )
        );

        if (debug && !cellsToRefine.empty())
        {
            Pout<< "    cellWeights too large in " << cellsToRefine.size()
                << " cells" << endl;
        }

        forAll(newCellsToRefine, nCTRI)
        {
            label cellI = newCellsToRefine[nCTRI];

            cellWeights.internalField()[cellI] /= 8.0;
        }

        // Mesh changing engine.
        polyTopoChange meshMod(mesh_);

        // Play refinement commands into mesh changer.
        meshCutter_.setRefinement(newCellsToRefine, meshMod);

        // Create mesh, return map from old to new mesh.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh
        (
            mesh_,
            false,  // inflate
            true,   // syncParallel
            true,   // orderCells (to reduce cell motion)
            false   // orderPoints
        );

        // Update fields
        mesh_.updateMesh(map);

        // Update numbering of cells/vertices.
        meshCutter_.updateMesh(map);

        {
            // Map cellVertices

            meshSearch cellSearch(mesh_);

            const labelList& reverseCellMap = map().reverseCellMap();

            List<DynamicList<point> > newCellVertices(mesh_.nCells());

            forAll(cellVertices, oldCellI)
            {
                DynamicList<point>& oldCellVertices =
                    cellVertices[oldCellI];

                if (findIndex(newCellsToRefine, oldCellI) >= 0)
                {
                    // This old cell was refined so the cell for the vertices
                    // in the new mesh needs to be searched for.

                    forAll (oldCellVertices, oPI)
                    {
                        const point& v = oldCellVertices[oPI];

                        label newCellI = cellSearch.findCell(v);

                        if (newCellI == -1)
                        {
                            Pout<< "findCell backgroundMeshDecomposition "
                                << v << " "
                                << oldCellI
                                << newCellI
                                << endl;
                        }

                        newCellVertices[newCellI].append(v);
                    }
                }
                else
                {
                    label newCellI = reverseCellMap[oldCellI];

                    forAll(oldCellVertices, oPI)
                    {
                        const point& v = oldCellVertices[oPI];

                        newCellVertices[newCellI].append(v);
                    }
                }
            }

            cellVertices.transfer(newCellVertices);
        }

        if (debug)
        {
            Info<< "    Background mesh refined from "
                << returnReduce(map().nOldCells(), sumOp<label>())
                << " to " << mesh_.globalData().nTotalCells()
                << " cells." << endl;

            const_cast<Time&>(mesh_.time())++;
            Info<< "Time " << mesh_.time().timeName() << endl;
            cellWeights.write();
            mesh_.write();
        }
    }

    if (debug)
    {
        printMeshData(mesh_);

        Pout<< "    Pre distribute sum(cellWeights) "
            << sum(cellWeights.internalField())
            << " max(cellWeights) "
            << max(cellWeights.internalField())
            << endl;
    }

    labelList newDecomp = decomposerPtr_().decompose
    (
        mesh_,
        mesh_.cellCentres(),
        cellWeights
    );

    Info<< "    Redistributing background mesh cells" << endl;

    fvMeshDistribute distributor(mesh_, mergeDist_);

    autoPtr<mapDistributePolyMesh> mapDist = distributor.distribute(newDecomp);

    meshCutter_.distribute(mapDist);

    if (debug)
    {
        printMeshData(mesh_);

        Pout<< "    Post distribute sum(cellWeights) "
            << sum(cellWeights.internalField())
            << " max(cellWeights) "
            << max(cellWeights.internalField())
            << endl;

        const_cast<Time&>(mesh_.time())++;
        Info<< "Time " << mesh_.time().timeName() << endl;
        mesh_.write();
        cellWeights.write();
    }

    mapDist().distributeCellData(cellVertices);

    buildPatchAndTree();

    return mapDist;
}

void Foam::backgroundMeshDecomposition::distributePoints
(
    List<point>& points
) const
{
    labelList toProc(processorPosition(points));

    autoPtr<mapDistribute> map(buildMap(toProc));

    map().distribute(points);
}


bool Foam::backgroundMeshDecomposition::positionOnThisProcessor
(
    const point& pt
) const
{
    return
        bFTreePtr_().getVolumeType(pt)
     == indexedOctree<treeDataBPatch>::INSIDE;
}


Foam::boolList Foam::backgroundMeshDecomposition::positionOnThisProcessor
(
    const List<point>& pts
) const
{
    boolList posProc(pts.size(), true);

    forAll(pts, pI)
    {
        posProc[pI] = positionOnThisProcessor(pts[pI]);
    }

    return posProc;
}


Foam::pointIndexHit Foam::backgroundMeshDecomposition::findLine
(
    const point& start,
    const point& end
) const
{
    return bFTreePtr_().findLine(start, end);
}


Foam::labelList Foam::backgroundMeshDecomposition::processorPosition
(
    const List<point>& pts
) const
{
    DynamicList<label> toCandidateProc;
    DynamicList<point> testPoints;
    labelList ptBlockStart(pts.size(), -1);
    labelList ptBlockSize(pts.size(), -1);

    label nTotalCandidates = 0;

    forAll(pts, pI)
    {
        const point& pt = pts[pI];

        label nCandidates = 0;

        forAll(allBackgroundMeshBounds_, procI)
        {
            if (allBackgroundMeshBounds_[procI].contains(pt))
            {
                toCandidateProc.append(procI);
                testPoints.append(pt);

                nCandidates++;
            }
        }

        ptBlockStart[pI] = nTotalCandidates;
        ptBlockSize[pI] = nCandidates;

        nTotalCandidates += nCandidates;
    }

    // Needed for reverseDistribute
    label preDistributionToCandidateProcSize = toCandidateProc.size();

    autoPtr<mapDistribute> map(buildMap(toCandidateProc));

    map().distribute(testPoints);

    List<bool> pointOnCandidate(testPoints.size(), false);

    // Test candidate points on candidate processors
    forAll(testPoints, tPI)
    {
        pointOnCandidate[tPI] = positionOnThisProcessor(testPoints[tPI]);
    }

    map().reverseDistribute
    (
        preDistributionToCandidateProcSize,
        pointOnCandidate
    );

    labelList ptProc(pts.size(), -1);

    forAll(pts, pI)
    {
        // Extract the sub list of results for this point

        SubList<bool> ptProcResults
        (
            pointOnCandidate,
            ptBlockSize[pI],
            ptBlockStart[pI]
        );

        forAll(ptProcResults, pPRI)
        {
            if (ptProcResults[pPRI])
            {
                ptProc[pI] = toCandidateProc[ptBlockStart[pI] + pPRI];

                break;
            }
        }

        if (ptProc[pI] < 0)
        {
            // No processor was found
            FatalErrorIn
            (
                "Foam::labelList"
                "Foam::backgroundMeshDecomposition::processorPosition"
                "("
                    "const List<point>& pts"
                ") const"
            )
                << "The position " << pts[pI]
                << " was not found in any part of the background mesh."
                << exit(FatalError);
        }
    }

    return ptProc;
}


Foam::List<Foam::List<Foam::pointIndexHit> >
Foam::backgroundMeshDecomposition::intersectsProcessors
(
    const List<point>& starts,
    const List<point>& ends,
    bool includeOwnProcessor
) const
{
    DynamicList<label> toCandidateProc;
    DynamicList<point> testStarts;
    DynamicList<point> testEnds;
    labelList segmentBlockStart(starts.size(), -1);
    labelList segmentBlockSize(starts.size(), -1);

    label nTotalCandidates = 0;

    forAll(starts, sI)
    {
        const point& s = starts[sI];
        const point& e = ends[sI];

        // Dummy point for treeBoundBox::intersects
        point p(vector::zero);

        label nCandidates = 0;

        forAll(allBackgroundMeshBounds_, procI)
        {
            // It is assumed that the sphere in question overlaps the source
            // processor, so don't test it, unless includeOwnProcessor is true
            if
            (
                (includeOwnProcessor || procI != Pstream::myProcNo())
              && allBackgroundMeshBounds_[procI].intersects(s, e, p)
            )
            {
                toCandidateProc.append(procI);
                testStarts.append(s);
                testEnds.append(e);

                nCandidates++;
            }
        }

        segmentBlockStart[sI] = nTotalCandidates;
        segmentBlockSize[sI] = nCandidates;

        nTotalCandidates += nCandidates;
    }

    // Needed for reverseDistribute
    label preDistributionToCandidateProcSize = toCandidateProc.size();

    autoPtr<mapDistribute> map(buildMap(toCandidateProc));

    map().distribute(testStarts);
    map().distribute(testEnds);

    List<pointIndexHit> segmentIntersectsCandidate(testStarts.size());

    // Test candidate segments on candidate processors
    forAll(testStarts, sI)
    {
        const point& s = testStarts[sI];
        const point& e = testEnds[sI];

        // If the sphere finds some elements of the patch, then it overlaps
        segmentIntersectsCandidate[sI] = bFTreePtr_().findLine(s, e);
    }

    map().reverseDistribute
    (
        preDistributionToCandidateProcSize,
        segmentIntersectsCandidate
    );

    List<List<pointIndexHit> > segmentHitProcs(starts.size());

    // Working storage for assessing processors
    DynamicList<pointIndexHit> tmpProcHits;

    forAll(starts, sI)
    {
        tmpProcHits.clear();

        // Extract the sub list of results for this point

        SubList<pointIndexHit> segmentProcResults
        (
            segmentIntersectsCandidate,
            segmentBlockSize[sI],
            segmentBlockStart[sI]
        );

        forAll(segmentProcResults, sPRI)
        {
            if (segmentProcResults[sPRI].hit())
            {
                tmpProcHits.append(segmentProcResults[sPRI]);

                tmpProcHits.last().setIndex
                (
                    toCandidateProc[segmentBlockStart[sI] + sPRI]
                );
            }
        }

        segmentHitProcs[sI] = tmpProcHits;
    }

    return segmentHitProcs;
}


Foam::labelListList Foam::backgroundMeshDecomposition::overlapsProcessors
(
    const List<point>& centres,
    const List<scalar>& radiusSqrs,
    bool includeOwnProcessor
) const
{
    DynamicList<label> toCandidateProc;
    DynamicList<point> testCentres;
    DynamicList<scalar> testRadiusSqrs;
    labelList sphereBlockStart(centres.size(), -1);
    labelList sphereBlockSize(centres.size(), -1);

    label nTotalCandidates = 0;

    forAll(centres, sI)
    {
        const point& c = centres[sI];
        scalar rSqr = radiusSqrs[sI];

        label nCandidates = 0;

        forAll(allBackgroundMeshBounds_, procI)
        {
            // It is assumed that the sphere in question overlaps the source
            // processor, so don't test it, unless includeOwnProcessor is true
            if
            (
                (includeOwnProcessor || procI != Pstream::myProcNo())
             && allBackgroundMeshBounds_[procI].overlaps(c, rSqr)
            )
            {
                toCandidateProc.append(procI);
                testCentres.append(c);
                testRadiusSqrs.append(rSqr);

                nCandidates++;
            }
        }

        sphereBlockStart[sI] = nTotalCandidates;
        sphereBlockSize[sI] = nCandidates;

        nTotalCandidates += nCandidates;
    }

    // Needed for reverseDistribute
    label preDistributionToCandidateProcSize = toCandidateProc.size();

    autoPtr<mapDistribute> map(buildMap(toCandidateProc));

    map().distribute(testCentres);
    map().distribute(testRadiusSqrs);

    List<bool> sphereOverlapsCandidate(testCentres.size(), false);

    // Test candidate spheres on candidate processors
    forAll(testCentres, sI)
    {
        const point& c = testCentres[sI];
        scalar rSqr = testRadiusSqrs[sI];

        // If the sphere finds some elements of the patch, then it overlaps
        sphereOverlapsCandidate[sI] = !bFTreePtr_().findSphere(c, rSqr).empty();
    }

    map().reverseDistribute
    (
        preDistributionToCandidateProcSize,
        sphereOverlapsCandidate
    );

    labelListList sphereProcs(centres.size());

    // Working storage for assessing processors
    DynamicList<label> tmpProcs;

    forAll(centres, sI)
    {
        tmpProcs.clear();

        // Extract the sub list of results for this point

        SubList<bool> sphereProcResults
        (
            sphereOverlapsCandidate,
            sphereBlockSize[sI],
            sphereBlockStart[sI]
        );

        forAll(sphereProcResults, sPRI)
        {
            if (sphereProcResults[sPRI])
            {
                tmpProcs.append(toCandidateProc[sphereBlockStart[sI] + sPRI]);
            }
        }

        sphereProcs[sI] = tmpProcs;
    }

    return sphereProcs;
}


// ************************************************************************* //
