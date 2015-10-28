/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "meshRefinement.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "OBJstream.H"
#include "triSurfaceMesh.H"
#include "treeDataCell.H"
#include "searchableSurfaces.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshRefinement::markSurfaceGapRefinement
(
    const scalar planarCos,

    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    // Get the gap level for the shells
    const labelList maxLevel(shells_.maxGapLevel());

    label oldNRefine = nRefine;

    if (max(maxLevel) > 0)
    {
        // Use cached surfaceIndex_ to detect if any intersection. If so
        // re-intersect to determine level wanted.

        // Collect candidate faces
        // ~~~~~~~~~~~~~~~~~~~~~~~

        labelList testFaces(getRefineCandidateFaces(refineCell));

        // Collect segments
        // ~~~~~~~~~~~~~~~~

        pointField start(testFaces.size());
        pointField end(testFaces.size());
        labelList minLevel(testFaces.size());

        calcCellCellRays
        (
            neiCc,
            neiLevel,
            testFaces,
            start,
            end,
            minLevel
        );


        // Collect cells to test for inside/outside in shell
        labelList cellToCompact(mesh_.nCells(), -1);
        labelList bFaceToCompact(mesh_.nFaces()-mesh_.nInternalFaces(), -1);
        List<FixedList<label, 3> > shellGapInfo;
        List<volumeType> shellGapMode;
        {
            DynamicField<point> compactToCc(mesh_.nCells()/10);
            DynamicList<label> compactToLevel(compactToCc.capacity());
            forAll(testFaces, i)
            {
                label faceI = testFaces[i];
                label own = mesh_.faceOwner()[faceI];
                if (cellToCompact[own] == -1)
                {
                    cellToCompact[own] = compactToCc.size();
                    compactToCc.append(cellCentres[own]);
                    compactToLevel.append(cellLevel[own]);
                }
                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    if (cellToCompact[nei] == -1)
                    {
                        cellToCompact[nei] = compactToCc.size();
                        compactToCc.append(cellCentres[nei]);
                        compactToLevel.append(cellLevel[nei]);
                    }
                }
                else
                {
                    label bFaceI = faceI - mesh_.nInternalFaces();
                    if (bFaceToCompact[bFaceI] == -1)
                    {
                        bFaceToCompact[bFaceI] = compactToCc.size();
                        compactToCc.append(neiCc[bFaceI]);
                        compactToLevel.append(neiLevel[bFaceI]);
                    }
                }
            }

            shells_.findHigherGapLevel
            (
                compactToCc,
                compactToLevel,
                shellGapInfo,
                shellGapMode
            );
        }


        //OBJstream str(mesh_.time().timePath()/"markSurfaceRefinement.obj");
        //Info<< "Dumping rays to " << str.name( ) << endl;

        const List<FixedList<label, 3> >& extendedGapLevel =
            surfaces_.extendedGapLevel();
        const List<volumeType>& extendedGapMode =
            surfaces_.extendedGapMode();

        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        {
            labelList surface2;
            List<pointIndexHit> hit2;
            labelList region2;
            vectorField normal2;

            surfaces_.findNearestIntersection
            (
                identity(surfaces_.surfaces().size()),
                start,
                end,

                surface1,
                hit1,
                region1,
                normal1,

                surface2,
                hit2,
                region2,
                normal2
            );
        }


        start.setSize(2*surface1.size());
        end.setSize(2*surface1.size());
        labelList map(2*surface1.size());
        pointField start2(2*surface1.size());
        pointField end2(2*surface1.size());
        labelList cellMap(2*surface1.size());

        label compactI = 0;

        forAll(surface1, i)
        {
            label surfI = surface1[i];

            if (surfI != -1)
            {
                label globalRegionI = surfaces_.globalRegion(surfI, region1[i]);

                label faceI = testFaces[i];
                const point& surfPt = hit1[i].hitPoint();

                label own = mesh_.faceOwner()[faceI];
                if (cellToCompact[own] != -1)
                {
                    // Combine info from shell and surface
                    label compactI = cellToCompact[own];
                    FixedList<label, 3> gapInfo;
                    volumeType gapMode;
                    mergeGapInfo
                    (
                        shellGapInfo[compactI],
                        shellGapMode[compactI],
                        extendedGapLevel[globalRegionI],
                        extendedGapMode[globalRegionI],
                        gapInfo,
                        gapMode
                    );

                    const point& cc = cellCentres[own];
                    bool okRay = generateRay
                    (
                        false,
                        surfPt,
                        normal1[i],
                        gapInfo,
                        gapMode,
                        surfPt+((cc-surfPt)&normal1[i])*normal1[i],
                        cellLevel[own],

                        start[compactI],
                        end[compactI],
                        start2[compactI],
                        end2[compactI]
                    );
                    cellMap[compactI] = own;

                    if (okRay)
                    {
                        map[compactI++] = i;
                    }
                }
                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    if (cellToCompact[nei] != -1)
                    {
                        // Combine info from shell and surface
                        label compactI = cellToCompact[nei];
                        FixedList<label, 3> gapInfo;
                        volumeType gapMode;
                        mergeGapInfo
                        (
                            shellGapInfo[compactI],
                            shellGapMode[compactI],
                            extendedGapLevel[globalRegionI],
                            extendedGapMode[globalRegionI],
                            gapInfo,
                            gapMode
                        );

                        const point& cc = cellCentres[nei];
                        bool okRay = generateRay
                        (
                            false,
                            surfPt,
                            normal1[i],
                            gapInfo,
                            gapMode,
                            surfPt+((cc-surfPt)&normal1[i])*normal1[i],
                            cellLevel[nei],

                            start[compactI],
                            end[compactI],
                            start2[compactI],
                            end2[compactI]
                        );
                        cellMap[compactI] = nei;

                        if (okRay)
                        {
                            map[compactI++] = i;
                        }
                    }
                }
                else
                {
                    label bFaceI = faceI - mesh_.nInternalFaces();

                    if (bFaceToCompact[bFaceI] != -1)
                    {
                        // Combine info from shell and surface
                        label compactI = bFaceToCompact[bFaceI];
                        FixedList<label, 3> gapInfo;
                        volumeType gapMode;
                        mergeGapInfo
                        (
                            shellGapInfo[compactI],
                            shellGapMode[compactI],
                            extendedGapLevel[globalRegionI],
                            extendedGapMode[globalRegionI],
                            gapInfo,
                            gapMode
                        );

                        const point& cc = neiCc[bFaceI];
                        bool okRay = generateRay
                        (
                            false,
                            surfPt,
                            normal1[i],
                            gapInfo,
                            gapMode,
                            surfPt+((cc-surfPt)&normal1[i])*normal1[i],
                            neiLevel[bFaceI],

                            start[compactI],
                            end[compactI],
                            start2[compactI],
                            end2[compactI]
                        );
                        cellMap[compactI] = -1;

                        if (okRay)
                        {
                            map[compactI++] = i;
                        }
                    }
                }
            }
        }

        //Info<< "Retesting " << returnReduce(compactI, sumOp<label>())
        //    << " out of " << returnReduce(start.size(), sumOp<label>())
        //    << endl;

        start.setSize(compactI);
        end.setSize(compactI);
        start2.setSize(compactI);
        end2.setSize(compactI);
        map.setSize(compactI);
        cellMap.setSize(compactI);

        testFaces = UIndirectList<label>(testFaces, map)();
        minLevel = UIndirectList<label>(minLevel, map)();
        surface1.clear();
        hit1.clear();
        region1.clear();
        normal1 = UIndirectList<vector>(normal1, map)();

        // Do intersections in first direction
        labelList surfaceHit;
        vectorField surfaceNormal;
        surfaces_.findNearestIntersection
        (
            start,
            end,
            surfaceHit,
            surfaceNormal
        );
        {
            // Do intersections in second direction and merge
            labelList surfaceHit2;
            vectorField surfaceNormal2;
            surfaces_.findNearestIntersection
            (
                start2,
                end2,
                surfaceHit2,
                surfaceNormal2
            );
            forAll(surfaceHit, i)
            {
                if (surfaceHit[i] == -1 && surfaceHit2[i] != -1)
                {
                    surfaceHit[i] = surfaceHit2[i];
                    surfaceNormal[i] = surfaceNormal2[i];
                }
            }
        }


        forAll(surfaceHit, i)
        {
            label surfI = surfaceHit[i];

            if (surfI != -1)
            {
                // Found intersection with surface. Check opposite normal.

                label cellI = cellMap[i];

                if (cellI != -1 && mag(normal1[i]&surfaceNormal[i]) > planarCos)
                {
                    if
                    (
                       !markForRefine
                        (
                            surfI,
                            nAllowRefine,
                            refineCell[cellI],
                            nRefine
                        )
                    )
                    {
                        break;
                    }
                }
            }
        }

        if
        (
            returnReduce(nRefine, sumOp<label>())
          > returnReduce(nAllowRefine, sumOp<label>())
        )
        {
            Info<< "Reached refinement limit." << endl;
        }
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


//Foam::labelList Foam::meshRefinement::extractHits
//(
//    const List<pointIndexHit>& info
//)
//{
//    labelList compactMap(info.size());
//    label compactI = 0;
//    forAll(info, i)
//    {
//        if (info[i].hit())
//        {
//            compactMap[compactI++] = i;
//        }
//    }
//    compactMap.setSize(compactI);
//    return compactMap;
//}
//
//
//void Foam::meshRefinement::mergeHits
//(
//    const pointField& samples,
//    const List<pointIndexHit>& extraInfo,
//    List<pointIndexHit>& info
//)
//{
//    if (extraInfo.size() != info.size())
//    {
//        FatalErrorIn
//        (
//            "meshRefinement::mergeHits(const List<pointIndexHit>&"
//            ", List<pointIndexHit>& info)"
//        )   << "Sizes differ " << extraInfo.size() << ' ' << info.size()
//            << exit(FatalError);
//    }
//
//    // Merge oppositeInfo2 into oppositeInfo
//    forAll(extraInfo, i)
//    {
//        if (!info[i].hit())
//        {
//            if (extraInfo[i].hit())
//            {
//                info[i] = extraInfo[i];
//            }
//        }
//        else if (extraInfo[i].hit())
//        {
//            const point& nearPt = samples[i];   //nearInfo[i].hitPoint();
//
//            scalar d = magSqr(info[i].hitPoint()-nearPt);
//            scalar d2 = magSqr(extraInfo[i].hitPoint()-nearPt);
//
//            if (d2 < d)
//            {
//                info[i] = extraInfo[i];
//            }
//        }
//    }
//}
//
//
//Foam::tmp<Foam::pointField> Foam::meshRefinement::extractPoints
//(
//    const List<pointIndexHit>& info
//)
//{
//    tmp<pointField> tfld(new pointField(info.size()));
//    pointField& fld = tfld();
//
//    forAll(info, i)
//    {
//        fld[i] = info[i].rawPoint();
//    }
//    return tfld;
//}
//
//
//Foam::meshRefinement::findNearestOppositeOp::findNearestOppositeOp
//(
//    const indexedOctree<treeDataTriSurface>& tree,
//    const point& oppositePoint,
//    const vector& oppositeNormal,
//    const scalar minCos
//)
//:
//    tree_(tree),
//    oppositePoint_(oppositePoint),
//    oppositeNormal_(oppositeNormal),
//    minCos_(minCos)
//{}
//
//
//void Foam::meshRefinement::findNearestOppositeOp::operator()
//(
//    const labelUList& indices,
//    const point& sample,
//    scalar& nearestDistSqr,
//    label& minIndex,
//    point& nearestPoint
//) const
//{
//    const treeDataTriSurface& shape = tree_.shapes();
//    const triSurface& patch = shape.patch();
//    const pointField& points = patch.points();
//
//    forAll(indices, i)
//    {
//        const label index = indices[i];
//        const labelledTri& f = patch[index];
//
//        pointHit nearHit = f.nearestPoint(sample, points);
//        scalar distSqr = sqr(nearHit.distance());
//
//        if (distSqr < nearestDistSqr)
//        {
//            // Nearer. Check if
//            // - a bit way from other hit
//            // - in correct search cone
//            vector d(nearHit.rawPoint()-oppositePoint_);
//            scalar normalDist(d&oppositeNormal_);
//
//            if (normalDist > Foam::sqr(SMALL) && normalDist/mag(d) > minCos_)
//            {
//                nearestDistSqr = distSqr;
//                minIndex = index;
//                nearestPoint = nearHit.rawPoint();
//            }
//        }
//    }
//}
//
//
//void Foam::meshRefinement::searchCone
//(
//    const label surfI,
//    labelList& nearMap,                 // cells
//    scalarField& nearGap,               // gap size
//    List<pointIndexHit>& nearInfo,      // nearest point on surface
//    List<pointIndexHit>& oppositeInfo   // detected point on gap (or miss)
//) const
//{
//    const labelList& cellLevel = meshCutter_.cellLevel();
//    const pointField& cellCentres = mesh_.cellCentres();
//    const scalar edge0Len = meshCutter_.level0EdgeLength();
//
//    const labelList& surfaceIndices = surfaces_.surfaces();
//    const List<FixedList<label, 3> >& extendedGapLevel =
//        surfaces_.extendedGapLevel();
//    const List<volumeType>& extendedGapMode = surfaces_.extendedGapMode();
//
//
//    label geomI = surfaceIndices[surfI];
//    const searchableSurface& geom = surfaces_.geometry()[geomI];
//
//    const triSurfaceMesh& s = refCast<const triSurfaceMesh>(geom);
//    const indexedOctree<treeDataTriSurface>& tree = s.tree();
//
//
//    const scalar searchCos(Foam::cos(degToRad(30)));
//
//    // Normals for ray shooting and inside/outside detection
//    vectorField nearNormal;
//    geom.getNormal(nearInfo, nearNormal);
//    // Regions
//    labelList nearRegion;
//    geom.getRegion(nearInfo, nearRegion);
//
//
//    // Now loop over all near points and search in the half cone
//    labelList map(nearInfo.size());
//    label compactI = 0;
//
//    oppositeInfo.setSize(nearInfo.size());
//
//    forAll(nearInfo, i)
//    {
//        label globalRegionI =
//            surfaces_.globalRegion(surfI, nearRegion[i]);
//
//        // Get updated gap information now we have the region
//        label nGapCells = extendedGapLevel[globalRegionI][0];
//        label minLevel = extendedGapLevel[globalRegionI][1];
//        label maxLevel = extendedGapLevel[globalRegionI][2];
//        volumeType mode = extendedGapMode[globalRegionI];
//
//        label cellI = nearMap[i];
//        label cLevel = cellLevel[cellI];
//
//        if (cLevel >= minLevel && cLevel < maxLevel)
//        {
//            scalar cellSize = edge0Len/pow(2, cLevel);
//
//            // Update gap size
//            nearGap[i] = nGapCells*cellSize;
//
//            const point& nearPt = nearInfo[i].hitPoint();
//            vector v(cellCentres[cellI]-nearPt);
//            scalar magV = mag(v);
//
//            // Like with ray shooting we want to
//            // - find triangles up to nearGap away on the wanted side of the
//            //   surface
//            // - find triangles up to 0.5*cellSize away on the unwanted side
//            //   of the surface. This is for cells straddling the surface
//            //   where
//            //   the cell centre might be on the wrong side of the surface
//
//            // Tbd: check that cell centre is inbetween the gap hits
//            // (only if the cell is far enough away)
//
//            scalar posNormalSize = 0.0;
//            scalar negNormalSize = 0.0;
//
//            if (mode == volumeType::OUTSIDE)
//            {
//                posNormalSize = nearGap[i];
//                if (magV < 0.5*cellSize)
//                {
//                    negNormalSize = 0.5*cellSize;
//                }
//            }
//            else if (mode == volumeType::INSIDE)
//            {
//                if (magV < 0.5*cellSize)
//                {
//                    posNormalSize = 0.5*cellSize;
//                }
//                negNormalSize = nearGap[i];
//            }
//            else
//            {
//                posNormalSize = nearGap[i];
//                negNormalSize = nearGap[i];
//            }
//
//            // Test with positive normal
//            oppositeInfo[compactI] = tree.findNearest
//            (
//                nearPt,
//                sqr(posNormalSize),
//                findNearestOppositeOp
//                (
//                    tree,
//                    nearPt,
//                    nearNormal[i],
//                    searchCos
//                )
//            );
//
//            if (oppositeInfo[compactI].hit())
//            {
//                map[compactI++] = i;
//            }
//            else
//            {
//                // Test with negative normal
//                oppositeInfo[compactI] = tree.findNearest
//                (
//                    nearPt,
//                    sqr(negNormalSize),
//                    findNearestOppositeOp
//                    (
//                        tree,
//                        nearPt,
//                        -nearNormal[i],
//                        searchCos
//                    )
//                );
//
//                if (oppositeInfo[compactI].hit())
//                {
//                    map[compactI++] = i;
//                }
//            }
//        }
//    }
//
//    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
//        << " hits on the correct side out of "
//        << returnReduce(map.size(), sumOp<label>()) << endl;
//    map.setSize(compactI);
//    oppositeInfo.setSize(compactI);
//
//    nearMap = UIndirectList<label>(nearMap, map)();
//    nearGap = UIndirectList<scalar>(nearGap, map)();
//    nearInfo = UIndirectList<pointIndexHit>(nearInfo, map)();
//    nearNormal = UIndirectList<vector>(nearNormal, map)();
//
//    // Exclude hits which aren't opposite enough. E.g. you might find
//    // a point on a perpendicular wall - but this does not consistute a gap.
//    vectorField oppositeNormal;
//    geom.getNormal(oppositeInfo, oppositeNormal);
//
//    compactI = 0;
//    forAll(oppositeInfo, i)
//    {
//        if ((nearNormal[i] & oppositeNormal[i]) < -0.707)
//        {
//            map[compactI++] = i;
//        }
//    }
//
//    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
//        << " hits opposite the nearest out of "
//        << returnReduce(map.size(), sumOp<label>()) << endl;
//    map.setSize(compactI);
//
//    nearMap = UIndirectList<label>(nearMap, map)();
//    nearGap = UIndirectList<scalar>(nearGap, map)();
//    nearInfo = UIndirectList<pointIndexHit>(nearInfo, map)();
//    oppositeInfo = UIndirectList<pointIndexHit>(oppositeInfo, map)();
//}


bool Foam::meshRefinement::generateRay
(
    const bool useSurfaceNormal,

    const point& nearPoint,
    const vector& nearNormal,
    const FixedList<label, 3>& gapInfo,
    const volumeType& mode,

    const point& cc,
    const label cLevel,

    point& start,
    point& end,
    point& start2,
    point& end2
) const
{
    // We want to handle the following cases:
    // - surface: small gap (marked with 'surface'). gap might be
    //            on inside or outside of surface.
    // - A: cell well inside the gap.
    // - B: cell well outside the gap.
    // - C: cell straddling the gap. cell centre might be inside
    //      or outside
    //
    //       +---+
    //       | B |
    //       +---+
    //
    //            +------+
    //            |      |
    //            |   C  |
    //    --------|------|----surface
    //            +------+
    //
    //        +---+
    //        | A |
    //        +---+
    //
    //
    //    --------------------surface
    //
    // So:
    // - find nearest point on surface
    // - in situation A,B decide if on wanted side of surface
    // - detect if locally a gap (and the cell inside the gap) by
    //   shooting a ray from the point on the surface in the direction
    //   of
    //   - A,B: the cell centre
    //   - C: the surface normal and/or negative surface normal
    //   and see we hit anything
    //
    // Variations of this scheme:
    // - always shoot in the direction of the surface normal. This needs
    //   then an additional check to make sure the cell centre is
    //   somewhere inside the gap
    // - instead of ray shooting use a 'constrained' nearest search
    //   by e.g. looking inside a search cone (implemented in searchCone).
    //   The problem with this constrained nearest is that it still uses
    //   the absolute nearest point on each triangle and only afterwards
    //   checks if it is inside the search cone.


    // Decide which near points are good:
    // - with updated minLevel and maxLevel and nearGap make sure
    //   the cell is still a candidate
    //   NOTE: inside the gap the nearest point on the surface will
    //         be HALF the gap size - otherwise we would have found
    //         a point on the opposite side
    // - if the mode is both sides
    // - or if the hit is inside the current cell (situation 'C',
    //   magV < 0.5cellSize)
    // - or otherwise if on the correct side

    bool okRay = false;

    if (cLevel >= gapInfo[1] && cLevel < gapInfo[2])
    {
        scalar cellSize = meshCutter_.level0EdgeLength()/pow(2, cLevel);

        // Calculate gap size
        scalar nearGap = gapInfo[0]*cellSize;

        // Distance to nearest
        vector v(cc-nearPoint);
        scalar magV = mag(v);

        if (useSurfaceNormal || magV < 0.5*cellSize)
        {
            const vector& n = nearNormal;

            // Situation 'C' above: cell too close. Use surface
            // normal to shoot rays

            if (mode == volumeType::OUTSIDE)
            {
                start = nearPoint+1e-6*n;
                end = nearPoint+nearGap*n;
                // Small second vector just to make sure to
                // refine cell
                start2 = nearPoint-1e-6*n;
                end2 = nearPoint-0.5*cellSize*n;
            }
            else if (mode == volumeType::INSIDE)
            {
                start = nearPoint-1e-6*n;
                end = nearPoint-nearGap*n;
                // Small second vector just to make sure to
                // refine cell
                start2 = nearPoint+1e-6*n;
                end2 = nearPoint+0.5*cellSize*n;
            }
            else if (mode == volumeType::MIXED)
            {
                start = nearPoint+1e-6*n;
                end = nearPoint+nearGap*n;

                start2 = nearPoint-1e-6*n;
                end2 = nearPoint-nearGap*n;
            }

            okRay = true;
        }
        else
        {
            // Situation 'A' or 'B' above: cell well away. Test if
            // cell on correct side of surface and shoot ray through
            // cell centre. Note: no need to shoot ray in other
            // direction since we're trying to detect cell inside
            // the gap.

            scalar s = (v&nearNormal);

            if
            (
                (mode == volumeType::MIXED)
             || (mode == volumeType::OUTSIDE && s > SMALL)
             || (mode == volumeType::INSIDE && s < -SMALL)
            )
            {
                // Use vector through cell centre
                vector n(v/(magV+ROOTVSMALL));

                start = nearPoint+1e-6*n;
                end = nearPoint+nearGap*n;
                // Dummy second vector
                start2 = start;
                end2 = start2;

                okRay = true;
            }
        }
    }

    return okRay;
}


//void Foam::meshRefinement::shootRays
//(
//    const label surfI,
//    labelList& nearMap,
//    scalarField& nearGap,
//    List<pointIndexHit>& nearInfo,
//    List<pointIndexHit>& oppositeInfo
//) const
//{
//    const labelList& cellLevel = meshCutter_.cellLevel();
//    const scalar edge0Len = meshCutter_.level0EdgeLength();
//
//    const labelList& surfaceIndices = surfaces_.surfaces();
//    const List<FixedList<label, 3> >& extendedGapLevel =
//        surfaces_.extendedGapLevel();
//    const List<volumeType>& extendedGapMode = surfaces_.extendedGapMode();
//
//
//    label geomI = surfaceIndices[surfI];
//    const searchableSurface& geom = surfaces_.geometry()[geomI];
//
//
//    // Normals for ray shooting and inside/outside detection
//    vectorField nearNormal;
//    geom.getNormal(nearInfo, nearNormal);
//    // Regions
//    labelList nearRegion;
//    geom.getRegion(nearInfo, nearRegion);
//
//    labelList map(nearInfo.size());
//
//    pointField start(nearInfo.size());
//    pointField end(nearInfo.size());
//    pointField start2(nearInfo.size());
//    pointField end2(nearInfo.size());
//
//    label compactI = 0;
//    forAll(nearInfo, i)
//    {
//        label globalRegionI =
//            surfaces_.globalRegion(surfI, nearRegion[i]);
//
//        label cellI = nearMap[i];
//        label cLevel = cellLevel[cellI];
//        scalar cellSize = edge0Len/pow(2, cLevel);
//
//        // Update gap size
//        nearGap[i] = extendedGapLevel[globalRegionI][0]*cellSize;
//
//        // Construct one or two rays to test for oppositeness
//        bool okRay = generateRay
//        (
//            false,
//            nearInfo[i].hitPoint(),
//            nearNormal[i],
//            extendedGapLevel[globalRegionI],
//            extendedGapMode[globalRegionI],
//
//            cellCentres[cellI],
//            cellLevel[cellI],
//
//            start[compactI],
//            end[compactI],
//            start2[compactI],
//            end2[compactI]
//        );
//        if (okRay)
//        {
//            map[compactI++] = i;
//        }
//    }
//
//    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
//        << " hits on the correct side out of "
//        << returnReduce(map.size(), sumOp<label>()) << endl;
//    map.setSize(compactI);
//    start.setSize(compactI);
//    end.setSize(compactI);
//    start2.setSize(compactI);
//    end2.setSize(compactI);
//
//    nearMap = UIndirectList<label>(nearMap, map)();
//    nearGap = UIndirectList<scalar>(nearGap, map)();
//    nearInfo = UIndirectList<pointIndexHit>(nearInfo, map)();
//
//    geom.findLineAny(start, end, oppositeInfo);
//
//    List<pointIndexHit> oppositeInfo2;
//    geom.findLineAny(start2, end2, oppositeInfo2);
//    mergeHits(extractPoints(nearInfo), oppositeInfo2, oppositeInfo);
//}


void Foam::meshRefinement::selectGapCandidates
(
    const labelList& refineCell,
    const label nRefine,

    labelList& cellMap,
    List<FixedList<label, 3> >& shellGapInfo,
    List<volumeType>& shellGapMode
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    // Collect cells to test
    cellMap.setSize(cellLevel.size()-nRefine);
    label compactI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            cellMap[compactI++] = cellI;
        }
    }
    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
        << " unmarked cells out of "
        << mesh_.globalData().nTotalCells() << endl;
    cellMap.setSize(compactI);

    // Do test to see whether cells are inside/outside shell with
    // applicable specification (minLevel <= celllevel < maxLevel)
    shells_.findHigherGapLevel
    (
        pointField(cellCentres, cellMap),
        UIndirectList<label>(cellLevel, cellMap)(),
        shellGapInfo,
        shellGapMode
    );

    // Compact out hits

    labelList map(shellGapInfo.size());
    compactI = 0;
    forAll(shellGapInfo, i)
    {
        if (shellGapInfo[i][2] > 0)
        {
            map[compactI++] = i;
        }
    }

    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
        << " cells inside gap shells out of "
        << mesh_.globalData().nTotalCells() << endl;

    map.setSize(compactI);
    cellMap = UIndirectList<label>(cellMap, map)();
    shellGapInfo = UIndirectList<FixedList<label, 3> >(shellGapInfo, map)();
    shellGapMode = UIndirectList<volumeType>(shellGapMode, map)();
}


void Foam::meshRefinement::mergeGapInfo
(
    const FixedList<label, 3>& shellGapInfo,
    const volumeType shellGapMode,
    const FixedList<label, 3>& surfGapInfo,
    const volumeType surfGapMode,
    FixedList<label, 3>& gapInfo,
    volumeType& gapMode
) const
{
    if (surfGapInfo[0] == 0)
    {
        gapInfo = shellGapInfo;
        gapMode = shellGapMode;
    }
    else if (shellGapInfo[0] == 0)
    {
        gapInfo = surfGapInfo;
        gapMode = surfGapMode;
    }
    else
    {
        // Both specify a level. Does surface level win? Or does information
        // need to be merged?

        //gapInfo[0] = max(surfGapInfo[0], shellGapInfo[0]);
        //gapInfo[1] = min(surfGapInfo[1], shellGapInfo[1]);
        //gapInfo[2] = max(surfGapInfo[2], shellGapInfo[2]);
        gapInfo = surfGapInfo;
        gapMode = surfGapMode;
    }
}


Foam::label Foam::meshRefinement::markInternalGapRefinement
(
    const scalar planarCos,
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();
    const scalar edge0Len = meshCutter_.level0EdgeLength();

    const List<FixedList<label, 3> >& extendedGapLevel =
        surfaces_.extendedGapLevel();
    const List<volumeType>& extendedGapMode = surfaces_.extendedGapMode();

    // Get the gap level for the shells
    const labelList maxLevel(shells_.maxGapLevel());

    label oldNRefine = nRefine;

    if (max(maxLevel) > 0)
    {
        // Collect cells to test
        labelList cellMap;
        List<FixedList<label, 3> > shellGapInfo;
        List<volumeType> shellGapMode;
        selectGapCandidates
        (
            refineCell,
            nRefine,

            cellMap,
            shellGapInfo,
            shellGapMode
        );


        // Find nearest point and normal on the surfaces
        List<pointIndexHit> nearInfo;
        vectorField nearNormal;
        labelList nearSurface;
        labelList nearRegion;
        {
            // Now we have both the cell-level and the gap size information. Use
            // this to calculate the gap size
            scalarField gapSize(cellMap.size());
            forAll(cellMap, i)
            {
                label cellI = cellMap[i];
                scalar cellSize = edge0Len/pow(2, cellLevel[cellI]);
                gapSize[i] = shellGapInfo[i][0]*cellSize;
            }

            surfaces_.findNearestRegion
            (
                identity(surfaces_.surfaces().size()),
                pointField(cellCentres, cellMap),
                sqr(gapSize),
                nearSurface,
                nearInfo,
                nearRegion,
                nearNormal
            );
        }



        labelList map(nearInfo.size());
        pointField start(nearInfo.size());
        pointField end(nearInfo.size());
        pointField start2(nearInfo.size());
        pointField end2(nearInfo.size());

        label compactI = 0;

        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                label globalRegionI = surfaces_.globalRegion
                (
                    nearSurface[i],
                    nearRegion[i]
                );

                // Combine info from shell and surface
                FixedList<label, 3> gapInfo;
                volumeType gapMode;
                mergeGapInfo
                (
                    shellGapInfo[i],
                    shellGapMode[i],
                    extendedGapLevel[globalRegionI],
                    extendedGapMode[globalRegionI],
                    gapInfo,
                    gapMode
                );


                // Construct one or two rays to test for oppositeness
                bool okRay = generateRay
                (
                    false,
                    nearInfo[i].hitPoint(),
                    nearNormal[i],
                    gapInfo,
                    gapMode,

                    cellCentres[cellMap[i]],
                    cellLevel[cellMap[i]],

                    start[compactI],
                    end[compactI],
                    start2[compactI],
                    end2[compactI]
                );
                if (okRay)
                {
                    map[compactI++] = i;
                }
            }
        }

        Info<< "Selected " << returnReduce(compactI, sumOp<label>())
            << " cells for testing out of "
            << mesh_.globalData().nTotalCells() << endl;
        map.setSize(compactI);
        start.setSize(compactI);
        end.setSize(compactI);
        start2.setSize(compactI);
        end2.setSize(compactI);

        cellMap = UIndirectList<label>(cellMap, map)();
        nearNormal = UIndirectList<vector>(nearNormal, map)();
        shellGapInfo.clear();
        shellGapMode.clear();
        nearInfo.clear();
        nearSurface.clear();
        nearRegion.clear();


        // Do intersections in first direction
        labelList surfaceHit;
        vectorField surfaceNormal;
        surfaces_.findNearestIntersection
        (
            start,
            end,
            surfaceHit,
            surfaceNormal
        );
        {
            // Do intersections in second direction and merge
            labelList surfaceHit2;
            vectorField surfaceNormal2;
            surfaces_.findNearestIntersection
            (
                start2,
                end2,
                surfaceHit2,
                surfaceNormal2
            );
            forAll(surfaceHit, i)
            {
                if (surfaceHit[i] == -1 && surfaceHit2[i] != -1)
                {
                    surfaceHit[i] = surfaceHit2[i];
                    surfaceNormal[i] = surfaceNormal2[i];
                }
            }
        }


        forAll(surfaceHit, i)
        {
            label surfI = surfaceHit[i];

            if (surfI != -1 && mag(nearNormal[i]&surfaceNormal[i]) > planarCos)
            {
                // Found intersection with surface. Check opposite normal.

                label cellI = cellMap[i];

                if
                (
                   !markForRefine
                    (
                        surfI,
                        nAllowRefine,
                        refineCell[cellI],
                        nRefine
                    )
                )
                {
                    break;
                }
            }
        }

        if
        (
            returnReduce(nRefine, sumOp<label>())
          > returnReduce(nAllowRefine, sumOp<label>())
        )
        {
            Info<< "Reached refinement limit." << endl;
        }
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


Foam::label Foam::meshRefinement::markSmallFeatureRefinement
(
    const scalar planarCos,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    const labelList& surfaceIndices = surfaces_.surfaces();
    const List<FixedList<label, 3> >& extendedGapLevel =
        surfaces_.extendedGapLevel();
    const List<volumeType>& extendedGapMode = surfaces_.extendedGapMode();

    label oldNRefine = nRefine;

    // Check that we're using any gap refinement
    labelList shellMaxLevel(shells_.maxGapLevel());

    if (max(shellMaxLevel) == 0)
    {
        return 0;
    }

    //- Force calculation of tetBasePt
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellTree();


    forAll(surfaceIndices, surfI)
    {
        label geomI = surfaceIndices[surfI];
        const searchableSurface& geom = surfaces_.geometry()[geomI];


        // Get the element index in a roundabout way. Problem is e.g.
        // distributed surface where local indices differ from global
        // ones (needed for getRegion call)

        pointField ctrs;
        labelList region;
        vectorField normal;
        {
            // Representative local coordinates and bounding sphere
            scalarField radiusSqr;
            geom.boundingSpheres(ctrs, radiusSqr);

            List<pointIndexHit> info;
            geom.findNearest(ctrs, radiusSqr, info);

            forAll(info, i)
            {
                if (!info[i].hit())
                {
                    FatalErrorIn("meshRefinement::markSmallFeatureRefinement")
                        << "fc:" << ctrs[i]
                        << " radius:" << radiusSqr[i]
                        << exit(FatalError);
                }
            }

            geom.getRegion(info, region);
            geom.getNormal(info, normal);
        }

        // Do test to see whether triangles are inside/outside shell with
        // applicable specification (minLevel <= celllevel < maxLevel)
        List<FixedList<label, 3> > shellGapInfo;
        List<volumeType> shellGapMode;
        shells_.findHigherGapLevel
        (
            ctrs,
            labelList(ctrs.size(), 0),
            shellGapInfo,
            shellGapMode
        );


        labelList map(ctrs.size());
        labelList cellMap(ctrs.size());
        label compactI = 0;

        pointField start(ctrs.size());
        pointField end(ctrs.size());
        pointField start2(ctrs.size());
        pointField end2(ctrs.size());

        forAll(ctrs, i)
        {
            if (shellGapInfo[i][2] > 0)
            {
                label globalRegionI = surfaces_.globalRegion(surfI, region[i]);

                // Combine info from shell and surface
                FixedList<label, 3> gapInfo;
                volumeType gapMode;
                mergeGapInfo
                (
                    shellGapInfo[i],
                    shellGapMode[i],
                    extendedGapLevel[globalRegionI],
                    extendedGapMode[globalRegionI],
                    gapInfo,
                    gapMode
                );

                //- Option 1: use octree nearest searching inside polyMesh
                //label cellI = mesh_.findCell(pt);

                //- Option 2: use octree 'inside' searching inside polyMesh. Is
                //            much faster.
                label cellI = -1;
                const indexedOctree<treeDataCell>& tree = mesh_.cellTree();
                if (tree.nodes().size() && tree.bb().contains(ctrs[i]))
                {
                    cellI = tree.findInside(ctrs[i]);
                }

                if (cellI != -1 && refineCell[cellI] == -1)
                {
                    // Construct one or two rays to test for oppositeness
                    // Note that we always want to use the surface normal
                    // and not the vector from cell centre to surface point

                    bool okRay = generateRay
                    (
                        true,               // always use surface normal
                        ctrs[i],
                        normal[i],
                        gapInfo,
                        gapMode,

                        cellCentres[cellI],
                        cellLevel[cellI],

                        start[compactI],
                        end[compactI],
                        start2[compactI],
                        end2[compactI]
                    );
                    if (okRay)
                    {
                        cellMap[compactI] = cellI;
                        map[compactI++] = i;
                    }
                }
            }
        }

        Info<< "Selected " << returnReduce(compactI, sumOp<label>())
            << " cells containing triangle centres out of "
            << mesh_.globalData().nTotalCells() << endl;
        map.setSize(compactI);
        cellMap.setSize(compactI);
        start.setSize(compactI);
        end.setSize(compactI);
        start2.setSize(compactI);
        end2.setSize(compactI);

        ctrs.clear();
        region.clear();
        shellGapInfo.clear();
        shellGapMode.clear();
        normal = UIndirectList<vector>(normal, map)();


        // Do intersections in first direction.
        // Note passing in of -1 for cell level such that it does not
        // discard surfaces with level 0. (since can be overridden with
        // gap level)
        labelList surfaceHit;
        vectorField surfaceNormal;
        surfaces_.findNearestIntersection
        (
            start,
            end,
            surfaceHit,
            surfaceNormal
        );
        {
            // Do intersections in second direction and merge
            labelList surfaceHit2;
            vectorField surfaceNormal2;
            surfaces_.findNearestIntersection
            (
                start2,
                end2,
                surfaceHit2,
                surfaceNormal2
            );
            forAll(surfaceHit, i)
            {
                if (surfaceHit[i] == -1 && surfaceHit2[i] != -1)
                {
                    surfaceHit[i] = surfaceHit2[i];
                    surfaceNormal[i] = surfaceNormal2[i];
                }
            }
        }


        label nGaps = 0;

        forAll(surfaceHit, i)
        {
            label surfI = surfaceHit[i];

            if (surfI != -1 && mag(normal[i]&surfaceNormal[i]) > planarCos)
            {
                nGaps++;

                if
                (
                   !markForRefine
                    (
                        surfI,
                        nAllowRefine,
                        refineCell[cellMap[i]],
                        nRefine
                    )
                )
                {
                    break;
                }
            }
        }

        Info<< "For surface " << geom.name() << " found "
            << returnReduce(nGaps, sumOp<label>())
            << " cells in small gaps" << endl;

        if
        (
            returnReduce(nRefine, sumOp<label>())
          > returnReduce(nAllowRefine, sumOp<label>())
        )
        {
            Info<< "Reached refinement limit." << endl;
        }
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// ************************************************************************* //
