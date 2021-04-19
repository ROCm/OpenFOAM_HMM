/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
#include "fvMesh.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "removeCells.H"
#include "unitConversion.H"
#include "bitSet.H"

#include "FaceCellWave.H"
#include "volFields.H"
#include "wallPoints.H"
#include "searchableSurfaces.H"
#include "distributedTriSurfaceMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::label Foam::meshRefinement::markFakeGapRefinement
//(
//    const scalar planarCos,
//
//    const label nAllowRefine,
//    const labelList& neiLevel,
//    const pointField& neiCc,
//
//    labelList& refineCell,
//    label& nRefine
//) const
//{
//    label oldNRefine = nRefine;
//
//
//    // Collect candidate faces (i.e. intersecting any surface and
//    // owner/neighbour not yet refined.
//    const labelList testFaces(getRefineCandidateFaces(refineCell));
//
//    // Collect segments
//    pointField start(testFaces.size());
//    pointField end(testFaces.size());
//    labelList minLevel(testFaces.size());
//
//    calcCellCellRays
//    (
//        neiCc,
//        neiLevel,
//        testFaces,
//        start,
//        end,
//        minLevel
//    );
//
//
//    // Re-use the gap shooting methods. This needs:
//    //  - shell gapLevel : faked. Set to 0,labelMax
//    //  - surface gapLevel : faked by overwriting
//
//
//    List<FixedList<label, 3>>& surfGapLevel = const_cast
//    <
//        List<FixedList<label, 3>>&
//    >(surfaces_.extendedGapLevel());
//
//    List<volumeType>& surfGapMode = const_cast
//    <
//        List<volumeType>&
//    >(surfaces_.extendedGapMode());
//
//    const List<FixedList<label, 3>> surfOldLevel(surfGapLevel);
//    const List<volumeType> surfOldMode(surfGapMode);
//
//    // Set the extended gap levels
//    forAll(surfaces_.gapLevel(), regioni)
//    {
//        surfGapLevel[regioni] = FixedList<label, 3>
//        ({
//            3,
//           -1,
//            surfaces_.gapLevel()[regioni]+1
//        });
//    }
//    surfGapMode = volumeType::MIXED;
//
//Pout<< "gapLevel was:" << surfOldLevel << endl;
//Pout<< "gapLevel now:" << surfGapLevel << endl;
//Pout<< "gapMode was:" << surfOldMode << endl;
//Pout<< "gapMode now:" << surfGapMode << endl;
//Pout<< "nRefine was:" << oldNRefine << endl;
//
//
//
//    List<List<FixedList<label, 3>>>& shellGapLevel = const_cast
//    <
//        List<List<FixedList<label, 3>>>&
//    >(shells_.extendedGapLevel());
//
//    List<List<volumeType>>& shellGapMode = const_cast
//    <
//        List<List<volumeType>>&
//    >(shells_.extendedGapMode());
//
//    const List<List<FixedList<label, 3>>> shellOldLevel(shellGapLevel);
//    const List<List<volumeType>> shellOldMode(shellGapMode);
//
//    // Set the extended gap levels
//    forAll(shellGapLevel, shelli)
//    {
//        shellGapLevel[shelli] =  FixedList<label, 3>({3, -1, labelMax});
//        shellGapMode[shelli] = volumeType::MIXED;
//    }
//Pout<< "shellLevel was:" << shellOldLevel << endl;
//Pout<< "shellLevel now:" << shellGapLevel << endl;
//
//    const label nAdditionalRefined = markSurfaceGapRefinement
//    (
//        planarCos,
//
//        nAllowRefine,
//        neiLevel,
//        neiCc,
//
//        refineCell,
//        nRefine
//    );
//
//Pout<< "nRefine now:" << nRefine << endl;
//
//    // Restore
//    surfGapLevel = surfOldLevel;
//    surfGapMode = surfOldMode;
//    shellGapLevel = shellOldLevel;
//    shellGapMode = shellOldMode;
//
//    return nAdditionalRefined;
//}


void Foam::meshRefinement::markOutsideFaces
(
    const labelList& cellLevel,
    const labelList& neiLevel,
    const labelList& refineCell,
    bitSet& isOutsideFace
) const
{
    // Get faces:
    // - on outside of cell set
    // - inbetween same cell level (i.e. quads)

    isOutsideFace.setSize(mesh_.nFaces());
    isOutsideFace = Zero;

    forAll(mesh_.faceNeighbour(), facei)
    {
        label own = mesh_.faceOwner()[facei];
        label nei = mesh_.faceNeighbour()[facei];
        if
        (
            (cellLevel[own] == cellLevel[nei])
         && (
                (refineCell[own] != -1)
             != (refineCell[nei] != -1)
            )
        )
        {
            isOutsideFace.set(facei);
        }
    }
    {

        const label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();

        labelList neiRefineCell(nBnd);
        syncTools::swapBoundaryCellList(mesh_, refineCell, neiRefineCell);
        for (label bFacei = 0; bFacei < nBnd; ++bFacei)
        {
            label facei = mesh_.nInternalFaces()+bFacei;
            label own = mesh_.faceOwner()[facei];

            if
            (
                (cellLevel[own] == neiLevel[bFacei])
             && (
                    (refineCell[own] != -1)
                 != (neiRefineCell[bFacei] != -1)
                )
            )
            {
                isOutsideFace.set(facei);
            }
        }
    }
}


Foam::label Foam::meshRefinement::countFaceDirs
(
    const bitSet& isOutsideFace,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];
    const vectorField& faceAreas = mesh_.faceAreas();

    Vector<bool> haveDirs(vector::uniform(false));

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];

        if (isOutsideFace[facei])
        {
            const vector& n = faceAreas[facei];
            scalar magSqrN = magSqr(n);

            if (magSqrN > ROOTVSMALL)
            {
                for
                (
                    direction dir = 0;
                    dir < pTraits<vector>::nComponents;
                    dir++
                )
                {
                    if (Foam::sqr(n[dir]) > 0.99*magSqrN)
                    {
                        haveDirs[dir] = true;
                        break;
                    }
                }
            }
        }
    }

    label nDirs = 0;
    forAll(haveDirs, dir)
    {
        if (haveDirs[dir])
        {
            nDirs++;
        }
    }
    return nDirs;
}


void Foam::meshRefinement::growSet
(
    const labelList& neiLevel,
    const bitSet& isOutsideFace,
    labelList& refineCell,
    label& nRefine
) const
{
    // Get cells with three or more outside faces
    const cellList& cells = mesh_.cells();
    forAll(cells, celli)
    {
        if (refineCell[celli] == -1)
        {
            if (countFaceDirs(isOutsideFace, celli) == 3)
            {
                // Mark cell with any value
                refineCell[celli] = 0;
                nRefine++;
            }
        }
    }
}


void Foam::meshRefinement::markMultiRegionCell
(
    const label celli,
    const FixedList<label, 3>& surface,

    Map<FixedList<label, 3>>& cellToRegions,
    bitSet& isMultiRegion
) const
{
    if (!isMultiRegion[celli])
    {
        Map<FixedList<label, 3>>::iterator fnd = cellToRegions.find(celli);
        if (!fnd.found())
        {
            cellToRegions.insert(celli, surface);
        }
        else if (fnd() != surface)
        {
            // Found multiple intersections on cell
            isMultiRegion.set(celli);
        }
    }
}


void Foam::meshRefinement::detectMultiRegionCells
(
    const labelListList& faceZones,
    const labelList& testFaces,

    const labelList& surface1,
    const List<pointIndexHit>& hit1,
    const labelList& region1,

    const labelList& surface2,
    const List<pointIndexHit>& hit2,
    const labelList& region2,

    bitSet& isMultiRegion
) const
{
    isMultiRegion.clear();
    isMultiRegion.setSize(mesh_.nCells());

    Map<FixedList<label, 3>> cellToRegions(testFaces.size());

    forAll(testFaces, i)
    {
        const pointIndexHit& info1 = hit1[i];
        if (info1.hit())
        {
            const label facei = testFaces[i];
            const labelList& fz1 = faceZones[surface1[i]];
            const FixedList<label, 3> surfaceInfo1
            ({
                surface1[i],
                region1[i],
                (fz1.size() ? fz1[info1.index()] : region1[i])
            });

            markMultiRegionCell
            (
                mesh_.faceOwner()[facei],
                surfaceInfo1,
                cellToRegions,
                isMultiRegion
            );
            if (mesh_.isInternalFace(facei))
            {
                markMultiRegionCell
                (
                    mesh_.faceNeighbour()[facei],
                    surfaceInfo1,
                    cellToRegions,
                    isMultiRegion
                );
            }

            const pointIndexHit& info2 = hit2[i];

            if (info2.hit() && info1 != info2)
            {
                const labelList& fz2 = faceZones[surface2[i]];
                const FixedList<label, 3> surfaceInfo2
                ({
                    surface2[i],
                    region2[i],
                    (fz2.size() ? fz2[info2.index()] : region2[i])
                });

                markMultiRegionCell
                (
                    mesh_.faceOwner()[facei],
                    surfaceInfo2,
                    cellToRegions,
                    isMultiRegion
                );
                if (mesh_.isInternalFace(facei))
                {
                    markMultiRegionCell
                    (
                        mesh_.faceNeighbour()[facei],
                        surfaceInfo2,
                        cellToRegions,
                        isMultiRegion
                    );
                }
            }
        }
    }


    if (debug&meshRefinement::MESH)
    {
        volScalarField multiCell
        (
            IOobject
            (
                "multiCell",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, 1, 0, 0, 0),
                0.0
            )
        );
        forAll(isMultiRegion, celli)
        {
            if (isMultiRegion[celli])
            {
                multiCell[celli] = 1.0;
            }
        }

        Info<< "Writing all multi cells to " << multiCell.name() << endl;
        multiCell.write();
    }
}


Foam::label Foam::meshRefinement::markProximityRefinementWave
(
    const scalar planarCos,
    const labelList& blockedSurfaces,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    labelListList faceZones(surfaces_.surfaces().size());
    {
        // If triSurface do additional zoning based on connectivity
        for (const label surfi : blockedSurfaces)
        {
            const label geomi = surfaces_.surfaces()[surfi];
            const searchableSurface& s = surfaces_.geometry()[geomi];
            if (isA<triSurfaceMesh>(s) && !isA<distributedTriSurfaceMesh>(s))
            {
                const triSurfaceMesh& surf = refCast<const triSurfaceMesh>(s);
                const labelListList& edFaces = surf.edgeFaces();
                boolList isOpenEdge(edFaces.size(), false);
                forAll(edFaces, i)
                {
                    if (edFaces[i].size() == 1)
                    {
                        isOpenEdge[i] = true;
                    }
                }

                labelList faceZone;
                const label nZones = surf.markZones(isOpenEdge, faceZone);
                if (nZones > 1)
                {
                    faceZones[surfi].transfer(faceZone);
                }
            }
        }
    }


    // Re-work the blockLevel of the blockedSurfaces into a length scale
    // and a number of cells to cross
    scalarField regionToBlockSize(surfaces_.blockLevel().size(), 0);
    //label nIters = 2;

    for (const label surfi : blockedSurfaces)
    {
        const label geomi = surfaces_.surfaces()[surfi];
        const searchableSurface& s = surfaces_.geometry()[geomi];
        const label nRegions = s.regions().size();
        for (label regioni = 0; regioni < nRegions; regioni++)
        {
            const label globalRegioni = surfaces_.globalRegion(surfi, regioni);
            const label bLevel = surfaces_.blockLevel()[globalRegioni];
            regionToBlockSize[globalRegioni] =
                meshCutter_.level0EdgeLength()/pow(2.0, bLevel);

            //const label mLevel = surfaces_.maxLevel()[globalRegioni];
            //// TBD: check for higher cached level of surface due to vol
            ////      refinement. Problem: might still miss refinement bubble
            ////      fully inside thin channel
            //if (isA<triSurfaceMesh>(s) && !isA<distributedTriSurfaceMesh>(s))
            //{
            //    const triSurfaceMesh& surf = refCast<const triSurfaceMesh>(s);
            //}

            //nIters = max(nIters, (2<<(mLevel-bLevel)));
        }
    }
    // Clever limiting of the number of iterations (= max cells in the channel)
    // since it has too many problematic issues, e.g. with volume refinement
    // and the real check uses the proper distance anyway just disable.
    const label nIters = mesh_.globalData().nTotalCells();


    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined)
    const labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
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
    // TBD. correct nIters for actual cellLevel (since e.g. volume refinement
    //      might add to cell level). Unfortunately we only have minLevel,
    //      not maxLevel. Problem: what if volume refinement only in middle of
    //      channel? Workaround: have dummy surface with e.g. maxLevel 100 to
    //      force nIters to be high enough.


    // Test for all intersections (with surfaces of higher gap level than
    // minLevel) and cache per cell the max surface level and the local normal
    // on that surface.

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    vectorField normal1;

    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    vectorField normal2;

    surfaces_.findNearestIntersection
    (
        blockedSurfaces,
        start,
        end,

        surface1,
        hit1,
        region1,    // local region
        normal1,

        surface2,
        hit2,
        region2,    // local region
        normal2
    );


    // Detect cells that are using multiple surface regions
    bitSet isMultiRegion;
    detectMultiRegionCells
    (
        faceZones,
        testFaces,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2,

        isMultiRegion
    );


    label n = 0;
    forAll(testFaces, i)
    {
        if (hit1[i].hit())
        {
            n++;
        }
    }

    List<wallPoints> faceDist(n);
    labelList changedFaces(n);
    n = 0;

    DynamicList<point> originLocation(2);
    DynamicList<scalar> originDistSqr(2);
    DynamicList<FixedList<label, 3>> originSurface(2);
    //DynamicList<point> originNormal(2);


    //- To avoid walking through surfaces we mark all faces that have been
    //  intersected. We can either mark only those faces intersecting
    //  blockedSurfaces (i.e. with a 'blockLevel') or mark faces intersecting
    //  any (refinement) surface (this includes e.g. faceZones). This is
    //  much easier since that information is already cached
    //  (meshRefinement::intersectedFaces())

    //bitSet isBlockedFace(mesh_.nFaces());
    forAll(testFaces, i)
    {
        if (hit1[i].hit())
        {
            const label facei = testFaces[i];
            //isBlockedFace.set(facei);
            const point& fc = mesh_.faceCentres()[facei];
            const labelList& fz1 = faceZones[surface1[i]];

            originLocation.clear();
            originDistSqr.clear();
            //originNormal.clear();
            originSurface.clear();

            originLocation.append(hit1[i].hitPoint());
            originDistSqr.append(magSqr(fc-hit1[i].hitPoint()));
            //originNormal.append(normal1[i]);
            originSurface.append
            (
                FixedList<label, 3>
                ({
                    surface1[i],
                    region1[i],
                    (fz1.size() ? fz1[hit1[i].index()] : region1[i])
                })
            );

            if (hit2[i].hit() && hit1[i] != hit2[i])
            {
                const labelList& fz2 = faceZones[surface2[i]];
                originLocation.append(hit2[i].hitPoint());
                originDistSqr.append(magSqr(fc-hit2[i].hitPoint()));
                //originNormal.append(normal2[i]);
                originSurface.append
                (
                    FixedList<label, 3>
                    ({
                        surface2[i],
                        region2[i],
                        (fz2.size() ? fz2[hit2[i].index()] : region2[i])
                    })
                );
            }

            // Collect all seed data. Currently walking does not look at
            // surface direction - if so pass in surface normal as well
            faceDist[n] = wallPoints
            (
                originLocation,     // origin
                originDistSqr,      // distance to origin
                originSurface       // surface+region+zone
                //originNormal        // normal at origin
            );
            changedFaces[n] = facei;
            n++;
        }
    }


    // Clear intersection info
    surface1.clear();
    hit1.clear();
    region1.clear();
    normal1.clear();
    surface2.clear();
    hit2.clear();
    region2.clear();
    normal2.clear();


    List<wallPoints> allFaceInfo(mesh_.nFaces());
    List<wallPoints> allCellInfo(mesh_.nCells());

    // Any refinement surface (even a faceZone) should stop the gap walking.
    // This is exactly the information which is cached in the surfaceIndex_
    // field.
    const bitSet isBlockedFace(intersectedFaces());

    wallPoints::trackData td(isBlockedFace);
    FaceCellWave<wallPoints, wallPoints::trackData> wallDistCalc
    (
        mesh_,
        changedFaces,
        faceDist,
        allFaceInfo,
        allCellInfo,
        0,            // max iterations
        td
    );
    wallDistCalc.iterate(nIters);


    if (debug&meshRefinement::MESH)
    {
        // Dump current nearest opposite surfaces
        volScalarField distance
        (
            IOobject
            (
                "gapSize",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimLength,  //dimensionSet(0, 1, 0, 0, 0),
                -1
            )
        );

        forAll(allCellInfo, celli)
        {
            if (allCellInfo[celli].valid(wallDistCalc.data()))
            {
                const point& cc = mesh_.C()[celli];
                // Nearest surface points
                const List<point>& origin = allCellInfo[celli].origin();

                // Find 'opposite' pair with minimum distance
                for (label i = 0; i < origin.size(); i++)
                {
                    for (label j = i + 1; j < origin.size(); j++)
                    {
                        if (((cc-origin[i]) & (cc-origin[j])) < 0)
                        {
                            const scalar d(mag(origin[i]-origin[j]));
                            if (distance[celli] < 0)
                            {
                                distance[celli] = d;
                            }
                            else
                            {
                                distance[celli] = min(distance[celli], d);
                            }
                        }
                    }
                }
            }
        }
        distance.correctBoundaryConditions();

        Info<< "Writing measured gap distance to "
            << distance.name() << endl;
        distance.write();
    }



    // Detect tight gaps:
    // - cell is inbetween the two surfaces
    // - two surfaces are planarish
    // - two surfaces are not too far apart
    //   (number of walking iterations is a too-coarse measure)

    scalarField smallGapDistance(mesh_.nCells(), 0.0);
    label nMulti = 0;
    label nSmallGap = 0;

    //OBJstream str(mesh_.time().timePath()/"multiRegion.obj");


    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].valid(wallDistCalc.data()))
        {
            const point& cc = mesh_.C()[celli];

            const List<point>& origin = allCellInfo[celli].origin();
            const List<FixedList<label, 3>>& surface =
                allCellInfo[celli].surface();

            // Find pair with minimum distance
            for (label i = 0; i < origin.size(); i++)
            {
                for (label j = i + 1; j < origin.size(); j++)
                {
                    scalar maxDist
                    (
                        max
                        (
                            mag(cc-origin[i]),
                            mag(cc-origin[j])
                        )
                    );

                    //if (isMultiRegion[celli])
                    //{
                    //    // The intersection locations are too inaccurate
                    //    // (since not proper nearest, just a cell-cell ray
                    //    //  intersection) so include always
                    //
                    //    smallGapDistance[celli] =
                    //        max(smallGapDistance[celli], maxDist);
                    //
                    //
                    //    str.write(linePointRef(cc, origin[i]));
                    //    str.write(linePointRef(cc, origin[j]));
                    //
                    //    nMulti++;
                    //}
                    //else
                    if (((cc-origin[i]) & (cc-origin[j])) < 0)
                    {
                        const label globalRegioni = surfaces_.globalRegion
                        (
                            surface[i][0],
                            surface[i][1]
                        );
                        const label globalRegionj = surfaces_.globalRegion
                        (
                            surface[j][0],
                            surface[j][1]
                        );

                        const scalar maxSize = max
                        (
                            regionToBlockSize[globalRegioni],
                            regionToBlockSize[globalRegionj]
                        );

                        if
                        (
                            magSqr(origin[i]-origin[j])
                          < Foam::sqr(2*maxSize)
                        )
                        {
                            smallGapDistance[celli] =
                                max(smallGapDistance[celli], maxDist);
                            nSmallGap++;
                        }
                    }
                }
            }
        }
    }


    if (debug)
    {
        Info<< "Marked for blocking due to intersecting multiple surfaces  : "
            << returnReduce(nMulti, sumOp<label>()) << " cells." << endl;
        Info<< "Marked for blocking due to close opposite surfaces         : "
            << returnReduce(nSmallGap, sumOp<label>()) << " cells." << endl;
    }

    if (debug&meshRefinement::MESH)
    {
        volScalarField distance
        (
            IOobject
            (
                "smallGapDistance",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, 1, 0, 0, 0),
                0.0
            )
        );
        distance.field() = smallGapDistance;
        distance.correctBoundaryConditions();

        Info<< "Writing all small-gap cells to "
            << distance.name() << endl;
        distance.write();
    }


    // Mark refinement
    const label oldNRefine = nRefine;
    forAll(smallGapDistance, celli)
    {
        if (smallGapDistance[celli] > SMALL)
        {
            if
            (
                !markForRefine
                (
                    0,                      // mark level
                    nAllowRefine,
                    refineCell[celli],
                    nRefine
                )
            )
            {
                if (debug)
                {
                    Pout<< "Stopped refining since reaching my cell"
                        << " limit of " << mesh_.nCells()+7*nRefine
                        << endl;
                }
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

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::removeGapCells
(
    const scalar planarAngle,
    const labelList& minSurfaceLevel,
    const labelList& globalToMasterPatch,
    const label growIter
)
{
    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList refineCell(mesh_.nCells(), -1);
    label nRefine = 0;
    //markProximityRefinement
    //(
    //    Foam::cos(degToRad(planarAngle)),
    //
    //    minSurfaceLevel,                                // surface min level
    //    labelList(minSurfaceLevel.size(), labelMax),    // surfaceGapLevel
    //
    //    labelMax/Pstream::nProcs(), //nAllowRefine,
    //    neiLevel,
    //    neiCc,
    //
    //    refineCell,
    //    nRefine
    //);


    // Determine minimum blockLevel per surface
    Map<label> surfToBlockLevel;

    forAll(surfaces_.surfaces(), surfi)
    {
        const label geomi = surfaces_.surfaces()[surfi];
        const searchableSurface& s = surfaces_.geometry()[geomi];
        const label nRegions = s.regions().size();

        label minBlockLevel = labelMax;
        for (label regioni = 0; regioni < nRegions; regioni++)
        {
            const label globalRegioni = surfaces_.globalRegion(surfi, regioni);
            minBlockLevel = min
            (
                minBlockLevel,
                surfaces_.blockLevel()[globalRegioni]
            );
        }

        if (minBlockLevel < labelMax)
        {
            surfToBlockLevel.insert(surfi, minBlockLevel);
        }
    }


    markProximityRefinementWave
    (
        Foam::cos(degToRad(planarAngle)),
        surfToBlockLevel.sortedToc(),

        labelMax/Pstream::nProcs(), //nAllowRefine,
        neiLevel,
        neiCc,

        refineCell,
        nRefine
    );


    //// Mark big-gap refinement
    //markFakeGapRefinement
    //(
    //    Foam::cos(degToRad(planarAngle)),
    //
    //    labelMax/Pstream::nProcs(), //nAllowRefine,
    //    neiLevel,
    //    neiCc,
    //
    //    refineCell,
    //    nRefine
    //);


    Info<< "Marked for blocking due to close opposite surfaces         : "
        << returnReduce(nRefine, sumOp<label>()) << " cells." << endl;

    // Remove outliers, i.e. cells with all points exposed
    if (growIter)
    {
        labelList oldRefineCell(refineCell);

        // Pass1: extend the set to fill in gaps
        bitSet isOutsideFace;
        for (label iter = 0; iter < growIter; iter++)
        {
            // Get outside faces
            markOutsideFaces
            (
                meshCutter_.cellLevel(),
                neiLevel,
                refineCell,
                isOutsideFace
            );
            // Extend with cells with three outside faces
            growSet(neiLevel, isOutsideFace, refineCell, nRefine);
        }


        // Pass2: erode back to original set if pass1 didn't help
        for (label iter = 0; iter < growIter; iter++)
        {
            // Get outside faces. Ignore cell level.
            markOutsideFaces
            (
                labelList(mesh_.nCells(), 0),
                labelList(neiLevel.size(), 0),
                refineCell,
                isOutsideFace
            );

            // Unmark cells with three or more outside faces
            for (label celli = 0; celli < mesh_.nCells(); celli++)
            {
                if (refineCell[celli] != -1 && oldRefineCell[celli] == -1)
                {
                    if (countFaceDirs(isOutsideFace, celli) >= 3)
                    {
                        refineCell[celli] = -1;
                        --nRefine;
                    }
                }
            }
        }

        Info<< "Marked for blocking after filtering                        : "
            << returnReduce(nRefine, sumOp<label>()) << " cells." << endl;
    }


    // Determine patch for every mesh face
    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    labelList unnamedSurfaces(surfaceZonesInfo::getUnnamedSurfaces(surfZones));
    const label defaultRegion(surfaces_.globalRegion(unnamedSurfaces[0], 0));

    const labelList nearestRegion
    (
        nearestIntersection
        (
            unnamedSurfaces,
            defaultRegion
        )
    );

    // Pack
    labelList cellsToRemove(nRefine);
    nRefine = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell[cellI] != -1)
        {
            cellsToRemove[nRefine++] = cellI;
        }
    }

    // Remove cells
    removeCells cellRemover(mesh_);
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));

    labelList exposedPatches(exposedFaces.size());
    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        exposedPatches[i] = globalToMasterPatch[nearestRegion[facei]];
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


// ************************************************************************* //
