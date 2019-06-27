/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
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
#include "fvMesh.H"
#include "volFields.H"
#include "OBJstream.H"
#include "PatchTools.H"
#include "foamVtkSurfaceWriter.H"

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
    const PackedBoolList& isLeakPoint,
    const bool distanceMode,
    const point& origin
)
{
    const cell& cFaces2 = mesh.cells()[cellI];

    // 1. Get topologically nearest face

    label minDist = labelMax;
    label minFaceI = -1;
    label nMin = 0;
    forAll(cFaces2, i)
    {
        label faceI = cFaces2[i];
        const topoDistanceData& info = allFaceInfo[faceI];
        if (info.distance() < minDist)
        {
            minDist = info.distance();
            minFaceI = faceI;
            nMin = 1;
        }
        else if (info.distance() == minDist)
        {
            nMin++;
        }
    }

    if (nMin > 1)
    {
        // 2. Check all faces with minDist for minimum (or maximum)
        //    distance to origin
        if (distanceMode)
        {
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
        }
        else
        {
            // Avoid leak points
            label minLeakPoints = labelMax;
            forAll(cFaces2, i)
            {
                label faceI = cFaces2[i];
                if (allFaceInfo[faceI].distance() == minDist)
                {
                    // Count number of leak points
                    label nLeak = 0;
                    {
                        const face& f = mesh.faces()[faceI];
                        forAll(f, fp)
                        {
                            if (isLeakPoint[f[fp]])
                            {
                                nLeak++;
                            }
                        }
                    }

                    if (nLeak < minLeakPoints)
                    {
                        minLeakPoints = nLeak;
                        minFaceI  = faceI;
                    }
                }
            }
        }
    }

    return minFaceI;
}


void Foam::shortestPathSet::calculateDistance
(
    const label iter,
    const polyMesh& mesh,
    const label cellI,

    List<topoDistanceData>& allFaceInfo,
    List<topoDistanceData>& allCellInfo
) const
{
    int dummyTrackData = 0;

    // Seed faces on cell1
    DynamicList<topoDistanceData> faceDist;
    DynamicList<label> cFaces1;

    if (cellI != -1)
    {
        const labelList& cFaces = mesh.cells()[cellI];
        faceDist.reserve(cFaces.size());
        cFaces1.reserve(cFaces.size());

        for (label facei : cFaces)
        {
            if (!allFaceInfo[facei].valid(dummyTrackData))
            {
                cFaces1.append(facei);
                faceDist.append(topoDistanceData(123, 0));
            }
        }
    }



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


    if (debug)
    {
        const fvMesh& fm = refCast<const fvMesh>(mesh);

        const_cast<fvMesh&>(fm).setInstance(fm.time().timeName());
        fm.write();
        volScalarField fld
        (
            IOobject
            (
                "allCellInfo" + Foam::name(iter),
                fm.time().timeName(),
                fm,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fm,
            dimensionedScalar(dimless, Zero)
        );
        forAll(allCellInfo, celli)
        {
            fld[celli] = allCellInfo[celli].distance();
        }
        forAll(fld.boundaryField(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];
            SubList<topoDistanceData> p(pp.patchSlice(allFaceInfo));
            scalarField pfld(fld.boundaryField()[patchi].size());
            forAll(pfld, i)
            {
                pfld[i] = 1.0*p[i].distance();
            }
            fld.boundaryFieldRef()[patchi] == pfld;
        }
        //Note: do not swap cell values so do not do
        //fld.correctBoundaryConditions();
        fld.write();
    }
}


void Foam::shortestPathSet::sync
(
    const polyMesh& mesh,
    PackedBoolList& isLeakFace,
    PackedBoolList& isLeakPoint,
    const label celli,
    point& origin,
    bool& findMinDistance
) const
{
    syncTools::syncPointList
    (
        mesh,
        isLeakPoint,
        orEqOp<unsigned int>(),
        0u
    );
    syncTools::syncFaceList
    (
        mesh,
        isLeakFace,
        orEqOp<unsigned int>()
    );
    // Sync origin, findMinDistance
    {
        typedef Tuple2<label, Tuple2<point, bool>> ProcData;

        ProcData searchData
        (
            celli,
            Tuple2<point, bool>(origin, findMinDistance)
        );
        Foam::combineReduce
        (
            searchData,
            [](ProcData& x, const ProcData& y)
            {
                if (y.first() != -1)
                {
                    x = y;
                }
            }
        );
        origin = searchData.second().first();
        findMinDistance = searchData.second().second();
    }
}


bool Foam::shortestPathSet::touchesWall
(
    const polyMesh& mesh,
    const label facei,

    PackedBoolList& isLeakFace,
    const PackedBoolList& isLeakPoint
) const
{
    // Check if facei touches leakPoint
    const face& f = mesh.faces()[facei];
    forAll(f, fp)
    {
        if (isLeakPoint[f[fp]])
        {
            if (isLeakFace.set(facei))
            {
                return true;
            }
        }
    }

    return false;
}


void Foam::shortestPathSet::genSamples
(
    const bool markLeakPath,
    const label maxIter,
    const polyMesh& mesh,
    const boolList& isBlockedFace,
    const point& insidePoint,
    const label insideCelli,
    const point& outsidePoint,

    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist,
    PackedBoolList& isLeakCell,
    PackedBoolList& isLeakFace,
    PackedBoolList& isLeakPoint
) const
{
    const topoDistanceData maxData(labelMax, labelMax);

    // Get the target point
    label outsideCelli = mesh.findCell(outsidePoint);

    // Maintain overall track length. Used to make curveDist continuous.
    label trackLength = 0;

    for (label iter = 0; iter < maxIter; iter++)
    {
        List<topoDistanceData> allFaceInfo(mesh.nFaces());
        List<topoDistanceData> allCellInfo(mesh.nCells());

        // Mark blocked faces with high distance
        forAll(isBlockedFace, facei)
        {
            if (isBlockedFace[facei])
            {
                allFaceInfo[facei] = maxData;
            }
        }

        if (markLeakPath)
        {
            // Mark any previously found leak path. This marks all
            // faces of all cells on the path. This will make sure that
            // ultimately the path will go through another gap.
            forAll(isLeakCell, celli)
            {
                if (celli != insideCelli && celli != outsideCelli)
                {
                    if (isLeakCell[celli])
                    {
                        allCellInfo[celli] = maxData;
                        const cell& cFaces = mesh.cells()[celli];
                        for (auto facei : cFaces)
                        {
                            allFaceInfo[facei] = maxData;
                        }
                    }
                }
            }
        }

        // Mark any previously found leak faces. These are faces that
        // are (point-)connected to an existing boundary.
        forAll(isLeakFace, facei)
        {
            if (isLeakFace[facei])
            {
                allFaceInfo[facei] = maxData;
            }
        }


        // Pass1: Get distance to insideCelli

        calculateDistance(iter, mesh, insideCelli, allFaceInfo, allCellInfo);



        // Pass2: walk from outside points backwards. Note: could be done
        //        using FaceCellWave as well but is overly complex since
        //        does not allow logic comparing all faces of a cell.

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();


        bool targetFound = false;
        if (outsideCelli != -1)
        {
            int dummyTrackData;
            targetFound = allCellInfo[outsideCelli].valid(dummyTrackData);
            if (!targetFound)
            {
                WarningInFunction
                    << "Point " << outsidePoint
                    << " not reachable by walk from " << insidePoint
                    << ". Probably mesh has island/regions."
                    << " Skipped route detection." << endl;
            }
        }
        reduce(targetFound, orOp<bool>());
        if (!targetFound)
        {
            break;
        }


        // Start with given target cell and walk back
        // If point happens to be on multiple processors, random pick
        label frontCellI = outsideCelli;
        point origin(outsidePoint);
        bool findMinDistance = true;

        while (true)
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
                // e.g. if we start from cell with value = 4, we have
                // neighbour faces 4, 4, 5, 5. Choose 4 (least distance
                // to seed) and continue...

                frontFaceI = findMinFace
                (
                    mesh,
                    frontCellI,
                    allFaceInfo,
                    isLeakPoint,
                    findMinDistance,    // mode : find min or find max
                    origin
                );


                // Loop until we hit a boundary face
                PackedBoolList isNewLeakPoint(isLeakPoint);
                while (mesh.isInternalFace(frontFaceI))
                {
                    if (isBlockedFace.size() && isBlockedFace[frontFaceI])
                    {
                        // Pout<< " Found blocked face" << endl;
                        frontCellI = -1;
                        break;
                    }

                    // Step to neighbouring cell
                    label nbrCellI = mesh.faceOwner()[frontFaceI];
                    if (nbrCellI == frontCellI)
                    {
                        nbrCellI = mesh.faceNeighbour()[frontFaceI];
                    }

                    if (nbrCellI == insideCelli)
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
                        isLeakPoint,
                        findMinDistance,
                        origin
                    );

                    const topoDistanceData& cInfo = allCellInfo[frontCellI];
                    const topoDistanceData& fInfo = allFaceInfo[frontFaceI];

                    if (fInfo.distance() <= cInfo.distance())
                    {
                        samplingPts.append(mesh.cellCentres()[frontCellI]);
                        samplingCells.append(frontCellI);
                        samplingFaces.append(-1);
                        samplingSegments.append(iter);
                        samplingCurveDist.append
                        (
                            trackLength+cInfo.distance()
                        );
                        isLeakCell.set(frontCellI);

                        if
                        (
                            touchesWall
                            (
                                mesh,
                                frontFaceI,
                                isLeakFace,
                                isLeakPoint
                            )
                        )
                        {
                            isNewLeakPoint.set(mesh.faces()[frontFaceI]);
                            origin = mesh.faceCentres()[frontFaceI];
                            findMinDistance = false;
                        }

                    }
                }
                isLeakPoint.transfer(isNewLeakPoint);
            }

            // Situation 1: we found the destination cell (do nothing),
            // frontCellI is -1 on all processors
            if (returnReduce(frontCellI == -1, andOp<bool>()))
            {
                break;
            }

            // Situation 2: we're on a coupled patch and might need to
            //              switch processor/cell. We need to transfer:
            //              -frontface -frontdistance -leak points/faces
            boolList isFront(mesh.nBoundaryFaces(), false);

            if (frontFaceI != -1)
            {
                isFront[frontFaceI-mesh.nInternalFaces()] = true;
            }
            syncTools::swapBoundaryFaceList(mesh, isFront);

            label minCellDistance = labelMax;
            if (frontCellI != -1)
            {
                minCellDistance = allCellInfo[frontCellI].distance();
            }
            reduce(minCellDistance, minOp<label>());

            // Sync all leak data
            sync
            (
                mesh,
                isLeakFace,
                isLeakPoint,
                frontCellI,
                origin,
                findMinDistance
            );



            const label oldFrontFaceI = frontFaceI;
            frontCellI = -1;
            frontFaceI = -1;
            forAll(pbm, patchI)
            {
                const polyPatch& pp = pbm[patchI];
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;
                    if
                    (
                        oldFrontFaceI == -1
                     && isFront[faceI-mesh.nInternalFaces()]
                     && (isBlockedFace.empty() || !isBlockedFace[faceI])
                    )
                    {
                        frontFaceI = faceI;
                        frontCellI = pp.faceCells()[i];
                        break;
                    }
                }

                if
                (
                    frontCellI != -1
                 && allCellInfo[frontCellI].distance() < minCellDistance
                )
                {
                    const topoDistanceData& cInfo = allCellInfo[frontCellI];

                    samplingPts.append(mesh.cellCentres()[frontCellI]);
                    samplingCells.append(frontCellI);
                    samplingFaces.append(-1);
                    samplingSegments.append(iter);
                    samplingCurveDist.append
                    (
                        trackLength+cInfo.distance()
                    );
                    isLeakCell.set(frontCellI);

                    // Check if frontFacei touches leakPoint
                    if
                    (
                        touchesWall
                        (
                            mesh,
                            frontFaceI,
                            isLeakFace,
                            isLeakPoint
                        )
                    )
                    {
                        isLeakPoint.set(mesh.faces()[frontFaceI]);
                        origin = mesh.faceCentres()[frontFaceI];
                        findMinDistance = false;
                    }

                    break;
                }

                // When seed cell is isolated by processor boundaries
                if (insideCelli != -1 && frontCellI == insideCelli)
                {
                    // Pout<< " Found connection seed cell!" << endl;
                    frontCellI = -1;
                    break;
                }
            }

            // Sync all leak data
            sync
            (
                mesh,
                isLeakFace,
                isLeakPoint,
                frontCellI,
                origin,
                findMinDistance
            );
        }


        // Recalculate the overall trackLength
        trackLength = returnReduce
        (
            (
                samplingCurveDist.size()
              ? samplingCurveDist.last()
              : 0
            ),
            maxOp<scalar>()
        );


        if (debug)
        {
            const fvMesh& fm = refCast<const fvMesh>(mesh);

            const_cast<fvMesh&>(fm).setInstance(fm.time().timeName());
            volScalarField fld
            (
                IOobject
                (
                    "isLeakCell",
                    fm.time().timeName(),
                    fm,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fm,
                dimensionedScalar(dimless, Zero)
            );
            forAll(isLeakCell, celli)
            {
                fld[celli] = isLeakCell[celli];
            }
            fld.correctBoundaryConditions();
            fld.write();
        }
    }
}


void Foam::shortestPathSet::genSamples
(
    const bool markLeakPath,
    const label maxIter,
    const polyMesh& mesh,
    const labelUList& wallPatches,
    const boolList& isBlockedFace
)
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    // Seed faces and points on 'real' boundary
    PackedBoolList isLeakFace(mesh.nFaces());
    PackedBoolList isLeakPoint(mesh.nPoints());
    {
        // Real boundaries
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        for (label patchi : wallPatches)
        {
            const polyPatch& pp = pbm[patchi];
            forAll(pp, i)
            {
                isLeakPoint.set(pp[i]);
            }
        }

        // Meshed boundaries
        forAll(isBlockedFace, facei)
        {
            if (isBlockedFace[facei])
            {
                isLeakPoint.set(mesh.faces()[facei]);
            }
        }

        syncTools::syncPointList
        (
            mesh,
            isLeakPoint,
            orEqOp<unsigned int>(),
            0u
        );
    }


    // All cells along leak paths
    PackedBoolList isLeakCell(mesh.nCells());

    label prevSegmenti = 0;
    scalar prevDistance = 0.0;

    for (auto insidePoint : insidePoints_)
    {
        const label insideCelli = mesh.findCell(insidePoint);

        for (auto outsidePoint : outsidePoints_)
        {
            const label nOldSamples = samplingSegments.size();

            // Append any leak path to sampling*
            genSamples
            (
                markLeakPath,
                maxIter,
                mesh,
                isBlockedFace,
                insidePoint,
                insideCelli,
                outsidePoint,

                samplingPts,
                samplingCells,
                samplingFaces,
                samplingSegments,
                samplingCurveDist,
                isLeakCell,
                isLeakFace,
                isLeakPoint
            );

            // Make segment, distance consecutive
            label maxSegment = 0;
            scalar maxDistance = 0.0;
            for (label i = nOldSamples; i < samplingSegments.size(); ++i)
            {
                samplingSegments[i] += prevSegmenti;
                maxSegment = max(maxSegment, samplingSegments[i]);
                samplingCurveDist[i] += prevDistance;
                maxDistance = max(maxDistance, samplingCurveDist[i]);
            }
            prevSegmenti = returnReduce(maxSegment, maxOp<label>());
            prevDistance = returnReduce(maxDistance, maxOp<scalar>());
        }
    }

    if (debug)
    {
        const faceList leakFaces(mesh.faces(), isLeakFace.sortedToc());

        OBJstream str(mesh.time().path()/"isLeakFace.obj");
        str.write(leakFaces, mesh.points(), false);
        Pout<< "Writing " << leakFaces.size() << " faces to " << str.name()
            << endl;
    }


    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    // Move into *this
    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
    );

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shortestPathSet::shortestPathSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const bool markLeakPath,
    const label maxIter,
    const labelUList& wallPatches,
    const pointField& insidePoints,
    const pointField& outsidePoints,
    const boolList& isBlockedFace
)
:
    sampledSet(name, mesh, searchEngine, axis),
    insidePoints_(insidePoints),
    outsidePoints_(outsidePoints)
{
    if (debug)
    {
        fileName outputDir =
        (
            mesh.time().globalPath()
          / functionObject::outputPrefix
          / mesh.pointsInstance()
        );
        outputDir.clean();

        Info<< "shortestPathSet : Writing blocked faces to "
            << outputDir << endl;

        const indirectPrimitivePatch setPatch
        (
            IndirectList<face>
            (
                mesh.faces(),
                findIndices(isBlockedFace, true)
            ),
            mesh.points()
        );

        if (Pstream::parRun())
        {
            // Topological merge
            labelList pointToGlobal;
            labelList uniqueMeshPointLabels;
            autoPtr<globalIndex> globalPoints;
            autoPtr<globalIndex> globalFaces;
            faceList mergedFaces;
            pointField mergedPoints;
            Foam::PatchTools::gatherAndMerge
            (
                mesh,
                setPatch.localFaces(),
                setPatch.meshPoints(),
                setPatch.meshPointMap(),

                pointToGlobal,
                uniqueMeshPointLabels,
                globalPoints,
                globalFaces,

                mergedFaces,
                mergedPoints
            );

            // Write
            if (Pstream::master())
            {
                vtk::surfaceWriter writer
                (
                    mergedPoints,
                    mergedFaces,
                    (outputDir / "blockedFace"),
                    false  // serial only - already merged
                );

                writer.writeGeometry();
            }
        }
        else
        {
            vtk::surfaceWriter writer
            (
                setPatch.localPoints(),
                setPatch.localFaces(),
                (outputDir / "blockedFace"),
                false  // serial only - redundant
            );

            writer.writeGeometry();
        }
    }

    genSamples(markLeakPath, maxIter, mesh, wallPatches, isBlockedFace);
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
    insidePoints_(dict.get<pointField>("insidePoints")),
    outsidePoints_(dict.get<pointField>("outsidePoints"))
{
    const label maxIter(dict.lookupOrDefault<label>("maxIter", 50));
    const bool markLeakPath(dict.lookupOrDefault<bool>("markLeakPath", true));

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    DynamicList<label> wallPatches(pbm.size());
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (!pp.coupled() && !isA<emptyPolyPatch>(pp))
        {
            wallPatches.append(patchi);
        }
    }

    genSamples(markLeakPath, maxIter, mesh, wallPatches, boolList());
}


// ************************************************************************* //
