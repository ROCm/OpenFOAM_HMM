/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    All to do with snapping to the surface

\*----------------------------------------------------------------------------*/

#include "autoHexMeshDriver.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "mapPolyMesh.H"
#include "motionSmoother.H"
#include "pointEdgePoint.H"
#include "PointEdgeWave.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::autoHexMeshDriver::getZonedSurfaces
(
    labelList& zonedSurfaces,
    labelList& unzonedSurfaces
) const
{    // Surfaces with zone information
    const wordList& faceZoneNames = surfaces().faceZoneNames();

    zonedSurfaces.setSize(faceZoneNames.size());
    label zonedI = 0;
    unzonedSurfaces.setSize(faceZoneNames.size());
    label unzonedI = 0;

    forAll(faceZoneNames, surfI)
    {
        if (faceZoneNames[surfI].size() > 0)
        {
            zonedSurfaces[zonedI++] = surfI;
        }
        else
        {
            unzonedSurfaces[unzonedI++] = surfI;
        }
    }
    zonedSurfaces.setSize(zonedI);
    unzonedSurfaces.setSize(unzonedI);
}


// Get faces to repatch. Returns map from face to patch.
Foam::Map<Foam::label> Foam::autoHexMeshDriver::getZoneBafflePatches
(
    const bool allowBoundary
) const
{
    Map<label> bafflePatch(mesh_.nFaces()/1000);

    const wordList& faceZoneNames = surfaces().faceZoneNames();
    const faceZoneMesh& fZones = mesh_.faceZones();

    forAll(faceZoneNames, surfI)
    {
        if (faceZoneNames[surfI].size() > 0)
        {
            // Get zone
            label zoneI = fZones.findZoneID(faceZoneNames[surfI]);

            const faceZone& fZone = fZones[zoneI];

            //// Get patch allocated for zone
            //label patchI = surfaceToCyclicPatch_[surfI];
            // Get patch of (first region) of surface
            label patchI = globalToPatch_[surfaces().globalRegion(surfI, 0)];

            Info<< "For surface "
                << surfaces()[surfI].IOobject::name()
                //<< surfaces().names()[surfI]
                << " found faceZone " << fZone.name()
                << " and patch " << mesh_.boundaryMesh()[patchI].name()
                << endl;


            forAll(fZone, i)
            {
                label faceI = fZone[i];

                if (allowBoundary || mesh_.isInternalFace(faceI))
                {
                    if (!bafflePatch.insert(faceI, patchI))
                    {
                        label oldPatchI = bafflePatch[faceI];

                        if (oldPatchI != patchI)
                        {
                            FatalErrorIn("getZoneBafflePatches(const bool)")
                                << "Face " << faceI
                                << " fc:" << mesh_.faceCentres()[faceI]
                                << " is in faceZone "
                                << mesh_.boundaryMesh()[oldPatchI].name()
                                << " and in faceZone "
                                << mesh_.boundaryMesh()[patchI].name()
                                << abort(FatalError);
                        }
                    }
                }
            }
        }
    }
    return bafflePatch;
}


// Calculate geometrically collocated points, Requires PackedList to be
// sizes and initalised!
Foam::label Foam::autoHexMeshDriver::getCollocatedPoints
(
    const scalar tol,
    const pointField& points,
    PackedList<1>& isCollocatedPoint
)
{
    labelList pointMap;
    pointField newPoints;
    bool hasMerged = mergePoints
    (
        points,                         // points
        tol,                            // mergeTol
        false,                          // verbose
        pointMap,
        newPoints
    );

    if (!returnReduce(hasMerged, orOp<bool>()))
    {
        return 0;
    }

    // Determine which newPoints are referenced more than once
    label nCollocated = 0;

    // Per old point the newPoint. Or -1 (not set yet) or -2 (already seen
    // twice)
    labelList firstOldPoint(newPoints.size(), -1);
    forAll(pointMap, oldPointI)
    {
        label newPointI = pointMap[oldPointI];

        if (firstOldPoint[newPointI] == -1)
        {
            // First use of oldPointI. Store.
            firstOldPoint[newPointI] = oldPointI;
        }
        else if (firstOldPoint[newPointI] == -2)
        {
            // Third or more reference of oldPointI -> non-manifold
            isCollocatedPoint.set(oldPointI, 1u);
            nCollocated++;
        }
        else
        {
            // Second reference of oldPointI -> non-manifold
            isCollocatedPoint.set(firstOldPoint[newPointI], 1u);
            nCollocated++;

            isCollocatedPoint.set(oldPointI, 1u);
            nCollocated++;

            // Mark with special value to save checking next time round
            firstOldPoint[newPointI] = -2;
        }
    }
    return returnReduce(nCollocated, sumOp<label>());
}


// Calculate displacement as average of patch points.
Foam::pointField Foam::autoHexMeshDriver::smoothPatchDisplacement
(
    const motionSmoother& meshMover
) const
{
    const indirectPrimitivePatch& pp = meshMover.patch();

    // Calculate geometrically non-manifold points on the patch to be moved.
    PackedList<1> nonManifoldPoint(pp.nPoints());
    label nNonManifoldPoints = getCollocatedPoints
    (
        SMALL,
        pp.localPoints(),
        nonManifoldPoint
    );
    Info<< "Found " << nNonManifoldPoints << " non-mainfold point(s)."
        << endl;


    // Average points
    // ~~~~~~~~~~~~~~

    // We determine three points:
    // - average of (centres of) connected patch faces
    // - average of (centres of) connected internal mesh faces
    // - as fallback: centre of any connected cell
    // so we can do something moderately sensible for non/manifold points.

    // Note: the averages are calculated properly parallel. This is
    // necessary to get the points shared by processors correct.


    const labelListList& pointFaces = pp.pointFaces();
    const labelList& meshPoints = pp.meshPoints();
    const pointField& points = pp.points();
    const polyMesh& mesh = meshMover.mesh();


    // Get average position of boundary face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vectorField avgBoundary(pointFaces.size(), vector::zero);
    labelList nBoundary(pointFaces.size(), 0);

    forAll(pointFaces, patchPointI)
    {
        const labelList& pFaces = pointFaces[patchPointI];

        forAll(pFaces, pfI)
        {
            avgBoundary[patchPointI] += pp[pFaces[pfI]].centre(points);
        }
        nBoundary[patchPointI] = pFaces.size();
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        avgBoundary,
        plusEqOp<point>(),  // combine op
        vector::zero,       // null value
        false               // no separation
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nBoundary,
        plusEqOp<label>(),  // combine op
        0,                  // null value
        false               // no separation
    );

    forAll(avgBoundary, i)
    {
        avgBoundary[i] /= nBoundary[i];
    }


    // Get average position of internal face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vectorField avgInternal;
    labelList nInternal;
    {
        vectorField globalSum(mesh.nPoints(), vector::zero);
        labelList globalNum(mesh.nPoints(), 0);

        // Note: no use of pointFaces
        const faceList& faces = mesh.faces();

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            const face& f = faces[faceI];
            const point& fc = mesh.faceCentres()[faceI];

            forAll(f, fp)
            {
                globalSum[f[fp]] += fc;
                globalNum[f[fp]]++;
            }
        }

        // Count coupled faces as internal ones (but only once)
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(patches, patchI)
        {
            if (Pstream::parRun() && isA<processorPolyPatch>(patches[patchI]))
            {
                const processorPolyPatch& pp =
                    refCast<const processorPolyPatch>(patches[patchI]);

                if (pp.myProcNo() < pp.neighbProcNo())
                {
                    const vectorField::subField faceCentres = pp.faceCentres();

                    forAll(pp, i)
                    {
                        const face& f = pp[i];
                        const point& fc = faceCentres[i];

                        forAll(f, fp)
                        {
                            globalSum[f[fp]] += fc;
                            globalNum[f[fp]]++;
                        }
                    }
                }
            }
            else if (isA<cyclicPolyPatch>(patches[patchI]))
            {
                const cyclicPolyPatch& pp =
                    refCast<const cyclicPolyPatch>(patches[patchI]);

                const vectorField::subField faceCentres = pp.faceCentres();

                for (label i = 0; i < pp.size()/2; i++)
                {
                    const face& f = pp[i];
                    const point& fc = faceCentres[i];

                    forAll(f, fp)
                    {
                        globalSum[f[fp]] += fc;
                        globalNum[f[fp]]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            globalSum,
            plusEqOp<vector>(), // combine op
            vector::zero,       // null value
            false               // no separation
        );
        syncTools::syncPointList
        (
            mesh,
            globalNum,
            plusEqOp<label>(),  // combine op
            0,                  // null value
            false               // no separation
        );

        avgInternal.setSize(meshPoints.size());
        nInternal.setSize(meshPoints.size());

        forAll(avgInternal, patchPointI)
        {
            label meshPointI = meshPoints[patchPointI];

            nInternal[patchPointI] = globalNum[meshPointI];

            if (nInternal[patchPointI] == 0)
            {
                avgInternal[patchPointI] = globalSum[meshPointI];
            }
            else
            {
                avgInternal[patchPointI] =
                    globalSum[meshPointI]
                  / nInternal[patchPointI];
            }
        }
    }


    // Precalculate any cell using mesh point (replacement of pointCells()[])
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList anyCell(mesh.nPoints(), -1);
    forAll(mesh.faceNeighbour(), faceI)
    {
        label own = mesh.faceOwner()[faceI];
        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            anyCell[f[fp]] = own;
        }
    }
    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        label own = mesh.faceOwner()[faceI];

        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            anyCell[f[fp]] = own;
        }
    }


    // Displacement to calculate.
    pointField patchDisp(meshPoints.size(), vector::zero);

    forAll(pointFaces, i)
    {
        label meshPointI = meshPoints[i];
        const point& currentPos = pp.points()[meshPointI];

        // Now we have the two average points: avgBoundary and avgInternal
        // and how many boundary/internal faces connect to the point
        // (nBoundary, nInternal)
        // Do some blending between the two.
        // Note: the following section has some reasoning behind it but the
        // blending factors can be experimented with.

        point newPos;

        if (nonManifoldPoint.get(i) == 0u)
        {
            // Points that are manifold. Weight the internal and boundary
            // by their number of faces and blend with
            scalar internalBlend = 0.1;
            scalar blend = 0.1;

            point avgPos =
                (
                   internalBlend*nInternal[i]*avgInternal[i]
                  +(1-internalBlend)*nBoundary[i]*avgBoundary[i]
                )
              / (internalBlend*nInternal[i]+(1-internalBlend)*nBoundary[i]);

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else if (nInternal[i] == 0)
        {
            // Non-manifold without internal faces. Use any connected cell
            // as internal point instead. Use precalculated any cell to avoid
            // e.g. pointCells()[meshPointI][0]

            const point& cc = mesh.cellCentres()[anyCell[meshPointI]];

            scalar cellCBlend = 0.8;
            scalar blend = 0.1;

            point avgPos = (1-cellCBlend)*avgBoundary[i] + cellCBlend*cc;

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else
        {
            // Non-manifold point with internal faces connected to them
            scalar internalBlend = 0.9;
            scalar blend = 0.1;

            point avgPos =
                internalBlend*avgInternal[i]
              + (1-internalBlend)*avgBoundary[i];

            newPos = (1-blend)*avgPos + blend*currentPos;
        }

        patchDisp[i] = newPos - currentPos;
    }

    return patchDisp;
}


Foam::tmp<Foam::scalarField> Foam::autoHexMeshDriver::edgePatchDist
(
    const pointMesh& pMesh,
    const indirectPrimitivePatch& pp
)
{
    const polyMesh& mesh = pMesh();

    // Set initial changed points to all the patch points
    List<pointEdgePoint> wallInfo(pp.nPoints());

    forAll(pp.localPoints(), ppI)
    {
        wallInfo[ppI] = pointEdgePoint(pp.localPoints()[ppI], 0.0);
    }

    // Current info on points
    List<pointEdgePoint> allPointInfo(mesh.nPoints());

    // Current info on edges
    List<pointEdgePoint> allEdgeInfo(mesh.nEdges());

    PointEdgeWave<pointEdgePoint> wallCalc
    (
        pMesh,
        pp.meshPoints(),
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh.globalData().nTotalPoints()  // max iterations
    );

    // Copy edge values into scalarField
    tmp<scalarField> tedgeDist(new scalarField(mesh.nEdges()));
    scalarField& edgeDist = tedgeDist();

    forAll(allEdgeInfo, edgeI)
    {
        edgeDist[edgeI] = Foam::sqrt(allEdgeInfo[edgeI].distSqr());
    }


    //{
    //    // For debugging: dump to file
    //    pointScalarField pointDist
    //    (
    //        IOobject
    //        (
    //            "pointDist",
    //            mesh.DB().timeName(),
    //            mesh.DB(),
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE
    //        ),
    //        pMesh,
    //        dimensionedScalar("pointDist", dimless, 0.0)
    //    );
    //
    //    forAll(allEdgeInfo, edgeI)
    //    {
    //        scalar d = Foam::sqrt(allEdgeInfo[edgeI].distSqr());
    //
    //        const edge& e = mesh.edges()[edgeI];
    //
    //        pointDist[e[0]] += d;
    //        pointDist[e[1]] += d;
    //    }
    //    forAll(pointDist, pointI)
    //    {
    //        pointDist[pointI] /= mesh.pointEdges()[pointI].size();
    //    }
    //    Info<< "Writing patch distance to " << pointDist.name()
    //        << " at time " << mesh.DB().timeName() << endl;
    //
    //    pointDist.write();
    //}

    return tedgeDist;
}


void Foam::autoHexMeshDriver::dumpMove
(
    const fileName& fName,
    const pointField& meshPts,
    const pointField& surfPts
)
{
    // Dump direction of growth into file
    Pout<< nl << "Dumping move direction to " << fName << nl
        << "View this Lightwave-OBJ file with e.g. javaview" << nl
        << endl;

    OFstream nearestStream(fName);

    label vertI = 0;

    forAll(meshPts, ptI)
    {
        meshTools::writeOBJ(nearestStream, meshPts[ptI]);
        vertI++;

        meshTools::writeOBJ(nearestStream, surfPts[ptI]);
        vertI++;

        nearestStream<< "l " << vertI-1 << ' ' << vertI << nl;
    }
}


// Check whether all displacement vectors point outwards of patch. Return true
// if so.
bool Foam::autoHexMeshDriver::outwardsDisplacement
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp
)
{
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();

    forAll(pointFaces, pointI)
    {
        const labelList& pFaces = pointFaces[pointI];

        vector disp(patchDisp[pointI]);

        scalar magDisp = mag(disp);

        if (magDisp > SMALL)
        {
            disp /= magDisp;

            bool outwards = meshTools::visNormal(disp, faceNormals, pFaces);

            if (!outwards)
            {
                Warning<< "Displacement " << patchDisp[pointI]
                    << " at mesh point " << pp.meshPoints()[pointI]
                    << " coord " << pp.points()[pp.meshPoints()[pointI]]
                    << " points through the surrounding patch faces" << endl;
                return false;
            }
        }
        else
        {
            //? Displacement small but in wrong direction. Would probably be ok.
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapPolyMesh> Foam::autoHexMeshDriver::createZoneBaffles
(
    List<labelPair>& baffles
)
{
    labelList zonedSurfaces;
    labelList unzonedSurfaces;
    getZonedSurfaces(zonedSurfaces, unzonedSurfaces);

    autoPtr<mapPolyMesh> map;

    // No need to sync; all processors will have all same zonedSurfaces.
    if (zonedSurfaces.size() > 0)
    {
        // Split internal faces on interface surfaces
        Info<< "Converting zoned faces into baffles ..." << endl;

        // Get faces (internal only) to be baffled. Map from face to patch
        // label.
        Map<label> faceToPatch(getZoneBafflePatches(false));

        label nZoneFaces = returnReduce(faceToPatch.size(), sumOp<label>());
        if (nZoneFaces > 0)
        {
            // Convert into labelLists
            labelList ownPatch(mesh_.nFaces(), -1);
            forAllConstIter(Map<label>, faceToPatch, iter)
            {
                ownPatch[iter.key()] = iter();
            }

            // Create baffles. both sides same patch.
            map = meshRefinerPtr_().createBaffles(ownPatch, ownPatch);

            // Get pairs of faces created.
            // Just loop over faceMap and store baffle if we encounter a slave
            // face.

            baffles.setSize(faceToPatch.size());
            label baffleI = 0;

            const labelList& faceMap = map().faceMap();
            const labelList& reverseFaceMap = map().reverseFaceMap();

            forAll(faceMap, faceI)
            {
                label oldFaceI = faceMap[faceI];

                // Does face originate from face-to-patch
                Map<label>::const_iterator iter = faceToPatch.find(oldFaceI);

                if (iter != faceToPatch.end())
                {
                    label masterFaceI = reverseFaceMap[oldFaceI];
                    if (faceI != masterFaceI)
                    {
                        baffles[baffleI++] = labelPair(masterFaceI, faceI);
                    }
                }
            }

            if (baffleI != faceToPatch.size())
            {
                FatalErrorIn("autoHexMeshDriver::createZoneBaffles(..)")
                    << "Had " << faceToPatch.size() << " patches to create "
                    << " but encountered " << baffleI
                    << " slave faces originating from patcheable faces."
                    << abort(FatalError);
            }

            if (debug_)
            {
                const_cast<Time&>(mesh_.time())++;
                Pout<< "Writing baffled mesh to time " << mesh_.time().timeName()
                    << endl;
                mesh_.write();
            }
        }
        Info<< "Created " << nZoneFaces << " baffles in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::autoHexMeshDriver::mergeZoneBaffles
(
    const List<labelPair>& baffles
)
{
    labelList zonedSurfaces;
    labelList unzonedSurfaces;
    getZonedSurfaces(zonedSurfaces, unzonedSurfaces);

    autoPtr<mapPolyMesh> map;

    // No need to sync; all processors will have all same zonedSurfaces.
    label nBaffles = returnReduce(baffles.size(), sumOp<label>());
    if (zonedSurfaces.size() > 0 && nBaffles > 0)
    {
        // Merge any baffles
        Info<< "Converting " << nBaffles << " baffles back into zoned faces ..."
            << endl;

        map = meshRefinerPtr_().mergeBaffles(baffles);

        Info<< "Converted baffles in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }

    return map;
}


Foam::scalarField Foam::autoHexMeshDriver::calcSnapDistance
(
    const indirectPrimitivePatch& pp
) const
{
    const dictionary& snapDict = dict_.subDict("snapDict");
    // When to snap
    scalar snapTol(readScalar(snapDict.lookup("snapTol")));


    const edgeList& edges = pp.edges();
    const labelListList& pointEdges = pp.pointEdges();
    const pointField& localPoints = pp.localPoints();

    scalarField maxEdgeLen(localPoints.size(), -GREAT);

    forAll(pointEdges, pointI)
    {
        const labelList& pEdges = pointEdges[pointI];

        forAll(pEdges, pEdgeI)
        {
            const edge& e = edges[pEdges[pEdgeI]];

            scalar len = e.mag(localPoints);

            maxEdgeLen[pointI] = max(maxEdgeLen[pointI], len);
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pp.meshPoints(),
        maxEdgeLen,
        maxEqOp<scalar>(),  // combine op
        -GREAT,             // null value
        false               // no separation
    );

    return snapTol*maxEdgeLen;
}


// Invert globalToPatch_ to get the patches related to surfaces.
Foam::labelList Foam::autoHexMeshDriver::getSurfacePatches() const
{
    // Set of patches originating from surface
    labelHashSet surfacePatchSet(surfaces().size());

    forAll(globalToPatch_, i)
    {
        if (globalToPatch_[i] != -1)
        {
            surfacePatchSet.insert(globalToPatch_[i]);
        }
    }

    DynamicList<label> surfacePatches(surfacePatchSet.size());

    for (label patchI = 0; patchI < mesh_.boundaryMesh().size(); patchI++)
    {
        if (surfacePatchSet.found(patchI))
        {
            surfacePatches.append(patchI);
        }
    }
    return surfacePatches.shrink();
}


void Foam::autoHexMeshDriver::preSmoothPatch
(
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
) const
{
    const dictionary& snapDict = dict_.subDict("snapDict");
    // Smoothing iterations
    label nSmoothPatch(readLabel(snapDict.lookup("nSmoothPatch")));
    // Snapping iterations
    label nSnap(readLabel(snapDict.lookup("nSnap")));

    labelList checkFaces;

    Info<< "Smoothing patch points ..." << endl;
    for (label smoothIter = 0; smoothIter < nSmoothPatch; smoothIter++)
    {
        Info<< "Smoothing iteration " << smoothIter << endl;
        checkFaces.setSize(mesh_.nFaces());
        forAll(checkFaces, faceI)
        {
            checkFaces[faceI] = faceI;
        }

        pointField patchDisp(smoothPatchDisplacement(meshMover));

        // The current mesh is the starting mesh to smooth from.
        meshMover.setDisplacement(patchDisp);
        meshMover.correct();

        scalar oldErrorReduction = -1;

        for (label snapIter = 0; snapIter < 2*nSnap; snapIter++)
        {
            Info<< nl << "Scaling iteration " << snapIter << endl;

            if (snapIter == nSnap)
            {
                Info<< "Displacement scaling for error reduction set to 0."
                    << endl;
                oldErrorReduction = meshMover.setErrorReduction(0.0);
            }

            // Try to adapt mesh to obtain displacement by smoothly
            // decreasing displacement at error locations.
            if (meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors))
            {
                Info<< "Successfully moved mesh" << endl;
                break;
            }
        }

        if (oldErrorReduction >= 0)
        {
            meshMover.setErrorReduction(oldErrorReduction);
        }
        Info<< endl;
    }


    // The current mesh is the starting mesh to smooth from.
    meshMover.correct();

    if (debug_)
    {
        Pout<< "Writing patch smoothed mesh to time " << mesh_.time().timeName()
            << endl;
        mesh_.write();
    }

    Info<< "Patch points smoothed in = "
        << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
}


Foam::vectorField Foam::autoHexMeshDriver::calcNearestSurface
(
    const scalarField& snapDist,
    motionSmoother& meshMover
) const
{
    Info<< "Calculating patchDisplacement as distance to nearest surface"
        << " point ..." << endl;

    const indirectPrimitivePatch& pp = meshMover.patch();
    const pointField& localPoints = pp.localPoints();

    // Divide surfaces into zoned and unzoned
    labelList zonedSurfaces;
    labelList unzonedSurfaces;
    getZonedSurfaces(zonedSurfaces, unzonedSurfaces);

    // Displacement per patch point
    vectorField patchDisp(localPoints.size(), vector::zero);


    // 1. All points to non-interface surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(localPoints, pointI)
    {
        pointIndexHit pHit;

        label surfI = surfaces().findNearest
        (
            unzonedSurfaces,
            localPoints[pointI],
            sqr(4*snapDist[pointI]),            // sqr of attract distance
            pHit
        );

        if (surfI != -1)
        {
            patchDisp[pointI] = pHit.hitPoint() - localPoints[pointI];
        }
        //else
        //{
        //   WarningIn("autoHexMeshDriver::calcNearestSurface(..)")
        //        << "For point:" << pointI
        //        << " coordinate:" << localPoints[pointI]
        //        << " did not find any surface within:" << 4*snapDist[pointI]
        //        << " meter." << endl;
        //}
    }


    // 2. All points on zones to their respective surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Surfaces with zone information
    const wordList& faceZoneNames = surfaces().faceZoneNames();

    forAll(zonedSurfaces, i)
    {
        label zoneSurfI = zonedSurfaces[i];

        const labelList surfacesToTest(1, zoneSurfI);

        label zoneI = mesh_.faceZones().findZoneID(faceZoneNames[zoneSurfI]);

        const faceZone& fZone = mesh_.faceZones()[zoneI];

        forAll(fZone, i)
        {
            const face& f = mesh_.faces()[fZone[i]];

            forAll(f, fp)
            {
                label meshPointI = f[fp];

                Map<label>::const_iterator iter =
                    pp.meshPointMap().find(meshPointI);

                if (iter != pp.meshPointMap().end())
                {
                    label pointI = iter();

                    pointIndexHit pHit;

                    label surfI = surfaces().findNearest
                    (
                        surfacesToTest,
                        localPoints[pointI],
                        sqr(4*snapDist[pointI]),    // sqr of attract distance
                        pHit
                    );

                    if (surfI != -1)
                    {
                        patchDisp[pointI] =
                            pHit.hitPoint() - localPoints[pointI];
                    }
                    else
                    {
                        WarningIn("autoHexMeshDriver::calcNearestSurface(..)")
                            << "For point:" << pointI
                            << " coordinate:" << localPoints[pointI]
                            << " did not find any surface within:"
                            << 4*snapDist[pointI]
                            << " meter." << endl;
                    }
                }
            }
        }
    }


    {
        scalarField magDisp(mag(patchDisp));

        Info<< "Wanted displacement : average:"
            << gSum(magDisp)/returnReduce(patchDisp.size(), sumOp<label>())
            << " min:" << gMin(magDisp)
            << " max:" << gMax(magDisp) << endl;
    }

    Info<< "Calculated surface displacement in = "
        << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;


    // Limit amount of movement.
    forAll(patchDisp, patchPointI)
    {
        scalar magDisp = mag(patchDisp[patchPointI]);

        if (magDisp > snapDist[patchPointI])
        {
            patchDisp[patchPointI] *= snapDist[patchPointI] / magDisp;

            Pout<< "Limiting displacement for " << patchPointI
                << " from " << magDisp << " to " << snapDist[patchPointI]
                << endl;
        }
    }

    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh_,
        pp.meshPoints(),
        patchDisp,
        minMagEqOp(),                   // combine op
        vector(GREAT, GREAT, GREAT),    // null value
        false                           // no separation
    );


    // Check for displacement being outwards.
    outwardsDisplacement(pp, patchDisp);

    // Set initial distribution of displacement field (on patches) from
    // patchDisp and make displacement consistent with b.c. on displacement
    // pointVectorField.
    meshMover.setDisplacement(patchDisp);

    if (debug_)
    {
        dumpMove
        (
            mesh_.time().path()/"patchDisplacement.obj",
            pp.localPoints(),
            pp.localPoints() + patchDisp
        );
    }

    return patchDisp;
}


void Foam::autoHexMeshDriver::smoothDisplacement(motionSmoother& meshMover)
 const
{
    const pointMesh& pMesh = meshMover.pMesh();
    const indirectPrimitivePatch& pp = meshMover.patch();

    Info<< "Smoothing displacement ..." << endl;

    const dictionary& snapDict = dict_.subDict("snapDict");
    // Smoothing iterations
    label nSmoothDisp(readLabel(snapDict.lookup("nSmoothDispl")));

    // Set edge diffusivity as inverse of distance to patch
    scalarField edgeGamma(1.0/(edgePatchDist(pMesh, pp) + SMALL));
    //scalarField edgeGamma(mesh_.nEdges(), 1.0);
    //scalarField edgeGamma(wallGamma(mesh, pp, 10, 1));

    // Get displacement field
    pointVectorField& disp = meshMover.displacement();

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        if ((iter % 10) == 0)
        {
            Info<< "Iteration " << iter << endl;
        }
        pointVectorField oldDisp(disp);

        meshMover.smooth(oldDisp, edgeGamma, false, disp);
    }
    Info<< "Displacement smoothed in = "
        << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug_)
    {
        Pout<< "Writing smoothed mesh to time " << mesh_.time().timeName()
            << endl;
        mesh_.write();

        Pout<< "Writing displacement field ..." << endl;
        disp.write();
        tmp<pointScalarField> magDisp(mag(disp));
        magDisp().write();

        Pout<< "Writing actual patch displacement ..." << endl;
        vectorField actualPatchDisp
        (
            IndirectList<point>(disp, pp.meshPoints())()
        );
        dumpMove
        (
            mesh_.time().path()/"actualPatchDisplacement.obj",
            pp.localPoints(),
            pp.localPoints() + actualPatchDisp
        );
    }
}


void Foam::autoHexMeshDriver::scaleMesh
(
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
    const dictionary& snapDict = dict_.subDict("snapDict");
    // Snapping iterations
    label nSnap(readLabel(snapDict.lookup("nSnap")));


    // Relax displacement until correct mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList checkFaces(identity(mesh_.nFaces()));

    scalar oldErrorReduction = -1;

    Info<< "Moving mesh ..." << endl;
    for (label iter = 0; iter < 2*nSnap; iter++)
    {
        Info<< nl << "Iteration " << iter << endl;

        if (iter == nSnap)
        {
            Info<< "Displacement scaling for error reduction set to 0." << endl;
            oldErrorReduction = meshMover.setErrorReduction(0.0);
        }

        if (meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors))
        {
            Info<< "Successfully moved mesh" << endl;

            break;
        }
        if (debug_)
        {
            Pout<< "Writing scaled mesh to time " << mesh_.time().timeName()
                << endl;
            mesh_.write();

            Pout<< "Writing displacement field ..." << endl;
            meshMover.displacement().write();
            tmp<pointScalarField> magDisp(mag(meshMover.displacement()));
            magDisp().write();

            const_cast<Time&>(mesh_.time())++;
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover.setErrorReduction(oldErrorReduction);
    }
    Info<< "Moved mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
}


// ************************************************************************* //
