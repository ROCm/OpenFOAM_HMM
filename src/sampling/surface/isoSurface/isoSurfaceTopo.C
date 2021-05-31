/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "isoSurfaceTopo.H"
#include "polyMesh.H"
#include "volFields.H"
#include "tetCell.H"
#include "tetMatcher.H"
#include "tetPointRef.H"
#include "DynamicField.H"
#include "syncTools.H"
#include "uindirectPrimitivePatch.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "isoSurfaceBaseMethods.H"
defineIsoSurfaceInterpolateMethods(Foam::isoSurfaceTopo);


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoSurfaceTopo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::isoSurfaceTopo::minTetQ
(
    const label facei,
    const label faceBasePtI
) const
{
    const scalar ownQuality =
        polyMeshTetDecomposition::minQuality
        (
            mesh_,
            mesh_.cellCentres()[mesh_.faceOwner()[facei]],
            facei,
            true,
            faceBasePtI
        );

    if (mesh_.isInternalFace(facei))
    {
        const scalar neiQuality =
            polyMeshTetDecomposition::minQuality
            (
                mesh_,
                mesh_.cellCentres()[mesh_.faceNeighbour()[facei]],
                facei,
                false,
                faceBasePtI
            );

        if (neiQuality < ownQuality)
        {
            return neiQuality;
        }
    }

    return ownQuality;
}


void Foam::isoSurfaceTopo::fixTetBasePtIs()
{
    // Determine points used by two faces on the same cell
    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();
    const labelList& faceOwn = mesh_.faceOwner();
    const labelList& faceNei = mesh_.faceNeighbour();


    // Get face triangulation base point
    tetBasePtIs_ = mesh_.tetBasePtIs();


    // Pre-filter: mark all cells with illegal base points
    bitSet problemCells(cells.size());

    forAll(tetBasePtIs_, facei)
    {
        if (tetBasePtIs_[facei] == -1)
        {
            problemCells.set(faceOwn[facei]);

            if (mesh_.isInternalFace(facei))
            {
                problemCells.set(faceNei[facei]);
            }
        }
    }


    // Mark all points that are shared by just two faces within an adjacent
    // problem cell as problematic
    bitSet problemPoints(mesh_.points().size());

    {
        // Number of times a point occurs in a cell.
        // Used to detect dangling vertices (count = 2)
        Map<label> pointCount;

        // Analyse problem cells for points shared by two faces only
        for (const label celli : problemCells)
        {
            pointCount.clear();

            for (const label facei : cells[celli])
            {
                for (const label pointi : faces[facei])
                {
                    ++pointCount(pointi);
                }
            }

            forAllConstIters(pointCount, iter)
            {
                if (iter.val() == 1)
                {
                    FatalErrorInFunction
                        << "point:" << iter.key()
                        << " at:" << mesh_.points()[iter.key()]
                        << " only used by one face" << nl
                        << exit(FatalError);
                }
                else if (iter.val() == 2)
                {
                    problemPoints.set(iter.key());
                }
            }
        }
    }


    // For all faces which form a part of a problem-cell, check if the base
    // point is adjacent to any problem points. If it is, re-calculate the base
    // point so that it is not.
    label nAdapted = 0;
    forAll(tetBasePtIs_, facei)
    {
        if
        (
            problemCells.test(faceOwn[facei])
         || (mesh_.isInternalFace(facei) && problemCells.test(faceNei[facei]))
        )
        {
            const face& f = faces[facei];

            // Check if either of the points adjacent to the base point is a
            // problem point. If not, the existing base point can be retained.
            const label fp0 = tetBasePtIs_[facei] < 0 ? 0 : tetBasePtIs_[facei];

            if
            (
                !problemPoints.test(f.rcValue(fp0))
             && !problemPoints.test(f.fcValue(fp0))
            )
            {
                continue;
            }

            // A new base point is required. Pick the point that results in the
            // least-worst tet and which is not adjacent to any problem points.
            scalar maxQ = -GREAT;
            label maxFp = -1;
            forAll(f, fp)
            {
                if
                (
                    !problemPoints.test(f.rcValue(fp))
                 && !problemPoints.test(f.fcValue(fp))
                )
                {
                    const scalar q = minTetQ(facei, fp);
                    if (q > maxQ)
                    {
                        maxQ = q;
                        maxFp = fp;
                    }
                }
            }

            if (maxFp != -1)
            {
                // Success! Set the new base point
                tetBasePtIs_[facei] = maxFp;
            }
            else
            {
                // No point was found on face that would not result in some
                // duplicate triangle. Do what? Continue and hope? Spit an
                // error? Silently or noisily reduce the filtering level?

                tetBasePtIs_[facei] = 0;
            }

            ++nAdapted;
        }
    }


    if (debug)
    {
        Pout<< "isoSurface : adapted starting point of triangulation on "
            << nAdapted << " faces." << endl;
    }

    syncTools::syncFaceList(mesh_, tetBasePtIs_, maxEqOp<label>());
}


Foam::label Foam::isoSurfaceTopo::generatePoint
(
    const label facei,
    const bool edgeIsDiag,
    const edge& vertices,

    DynamicList<edge>& pointToVerts,
    DynamicList<label>& pointToFace,
    DynamicList<bool>& pointFromDiag,
    EdgeMap<label>& vertsToPoint
) const
{
    // Generate new point, unless it already exists for edge

    label pointi = vertsToPoint.lookup(vertices, -1);
    if (pointi == -1)
    {
        pointi = pointToVerts.size();

        pointToVerts.append(vertices);
        pointToFace.append(facei);
        pointFromDiag.append(edgeIsDiag);
        vertsToPoint.insert(vertices, pointi);
    }

    return pointi;
}


void Foam::isoSurfaceTopo::generateTriPoints
(
    const label facei,
    const int tetCutIndex,
    const tetCell& tetLabels,
    const FixedList<bool, 6>& edgeIsDiag,// Per tet edge whether is face diag

    DynamicList<edge>& pointToVerts,
    DynamicList<label>& pointToFace,
    DynamicList<bool>& pointFromDiag,

    EdgeMap<label>& vertsToPoint,
    DynamicList<label>& verts       // Every three verts is new triangle
) const
{
    // Form the vertices of the triangles for each case
    switch (tetCutIndex)
    {
        case 0x00:
        case 0x0F:
        break;

        case 0x01:
        case 0x0E:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    tetLabels.edge(0),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    tetLabels.edge(1),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    tetLabels.edge(2),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (tetCutIndex == 0x0E)
            {
                // Flip normals
                const label sz = verts.size();
                std::swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x02:
        case 0x0D:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    tetLabels.reverseEdge(0),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    tetLabels.reverseEdge(3),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    tetLabels.edge(4),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (tetCutIndex == 0x0D)
            {
                // Flip normals
                const label sz = verts.size();
                std::swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x03:
        case 0x0C:
        {
            const label p0p2
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    tetLabels.edge(1),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            const label p1p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    tetLabels.reverseEdge(3),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    tetLabels.edge(2),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p1p3);
            verts.append(p0p2);
            verts.append(p1p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    tetLabels.edge(4),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p0p2);

            if (tetCutIndex == 0x0C)
            {
                // Flip normals
                const label sz = verts.size();
                std::swap(verts[sz-5], verts[sz-4]);
                std::swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x04:
        case 0x0B:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    tetLabels.reverseEdge(1),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    tetLabels.reverseEdge(4),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    tetLabels.reverseEdge(5),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (tetCutIndex == 0x0B)
            {
                // Flip normals
                const label sz = verts.size();
                std::swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x05:
        case 0x0A:
        {
            const label p0p1
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    tetLabels.edge(0),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            const label p2p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    tetLabels.reverseEdge(5),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append(p0p1);
            verts.append(p2p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    tetLabels.edge(2),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p0p1);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    tetLabels.edge(4),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p2p3);

            if (tetCutIndex == 0x0A)
            {
                // Flip normals
                const label sz = verts.size();
                std::swap(verts[sz-5], verts[sz-4]);
                std::swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x06:
        case 0x09:
        {
            const label p0p1
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    tetLabels.edge(0),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            const label p2p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    tetLabels.reverseEdge(5),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append(p0p1);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    tetLabels.reverseEdge(3),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p2p3);
            verts.append(p0p1);
            verts.append(p2p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    tetLabels.edge(1),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (tetCutIndex == 0x09)
            {
                // Flip normals
                const label sz = verts.size();
                std::swap(verts[sz-5], verts[sz-4]);
                std::swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x08:
        case 0x07:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    tetLabels.reverseEdge(2),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    tetLabels.reverseEdge(5),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    tetLabels.edge(3),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            if (tetCutIndex == 0x07)
            {
                // Flip normals
                const label sz = verts.size();
                std::swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;
    }
}


void Foam::isoSurfaceTopo::generateTriPoints
(
    const label celli,
    const bool isTet,

    DynamicList<edge>& pointToVerts,
    DynamicList<label>& pointToFace,
    DynamicList<bool>& pointFromDiag,

    EdgeMap<label>& vertsToPoint,
    DynamicList<label>& verts,
    DynamicList<label>& faceLabels
) const
{
    const faceList& faces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();
    const cell& cFaces = mesh_.cells()[celli];

    if (isTet)
    {
        // For tets don't do cell-centre decomposition, just use the
        // tet points and values

        const label startTrii = verts.size();

        const label facei = cFaces[0];
        const face& f0 = faces[facei];

        // Get the other point from f1. Tbd: check if not duplicate face
        // (ACMI / ignoreBoundaryFaces_).
        const face& f1 = faces[cFaces[1]];
        label oppositeI = -1;
        forAll(f1, fp)
        {
            oppositeI = f1[fp];
            if (!f0.found(oppositeI))
            {
                break;
            }
        }

        label p0 = f0[0];
        label p1 = f0[1];
        label p2 = f0[2];

        if (faceOwner[facei] == celli)
        {
            std::swap(p1, p2);
        }

        generateTriPoints
        (
            facei,
            getTetCutIndex
            (
                pVals_[p0],
                pVals_[p1],
                pVals_[p2],
                pVals_[oppositeI],
                iso_
            ),
            tetCell(p0, p1, p2, oppositeI),
            FixedList<bool, 6>(false),  // edgeIsDiag = false

            pointToVerts,
            pointToFace,
            pointFromDiag,
            vertsToPoint,
            verts       // Every three verts is new triangle
        );

        for (label nTris = (verts.size()-startTrii)/3; nTris; --nTris)
        {
            faceLabels.append(facei);
        }
    }
    else
    {
        for (const label facei : cFaces)
        {
            if
            (
               !mesh_.isInternalFace(facei)
             && ignoreBoundaryFaces_.test(facei-mesh_.nInternalFaces())
            )
            {
                continue;
            }

            const label startTrii = verts.size();

            const face& f = faces[facei];

            label fp0 = tetBasePtIs_[facei];

            // Fallback
            if (fp0 < 0)
            {
                fp0 = 0;
            }

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); ++i)
            {
                const label nextFp = f.fcIndex(fp);

                FixedList<bool, 6> edgeIsDiag(false);

                label p0 = f[fp0];
                label p1 = f[fp];
                label p2 = f[nextFp];
                if (faceOwner[facei] == celli)
                {
                    std::swap(p1, p2);
                    if (i != 2) edgeIsDiag[1] = true;
                    if (i != f.size()-1) edgeIsDiag[0] = true;
                }
                else
                {
                    if (i != 2) edgeIsDiag[0] = true;
                    if (i != f.size()-1) edgeIsDiag[1] = true;
                }

                generateTriPoints
                (
                    facei,
                    getTetCutIndex
                    (
                        pVals_[p0],
                        pVals_[p1],
                        pVals_[p2],
                        cVals_[celli],
                        iso_
                    ),
                    tetCell(p0, p1, p2, mesh_.nPoints()+celli),
                    edgeIsDiag,

                    pointToVerts,
                    pointToFace,
                    pointFromDiag,
                    vertsToPoint,
                    verts       // Every three verts is new triangle
                );

                fp = nextFp;
            }

            for (label nTris = (verts.size()-startTrii)/3; nTris; --nTris)
            {
                faceLabels.append(facei);
            }
        }
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::isoSurfaceTopo::triangulateOutside
(
    const bool filterDiag,
    const primitivePatch& pp,
    const boolUList& pointFromDiag,
    const labelUList& pointToFace,
    const label cellID,

    // outputs
    DynamicList<face>& compactFaces,
    DynamicList<label>& compactCellIDs
)
{
    // We can form pockets:
    // - 1. triangle on face
    // - 2. multiple triangles on interior (from diag edges)
    // - the edge loop will be pocket since it is only the diag
    //   edges that give it volume?

    // Retriangulate the exterior loops
    const labelListList& edgeLoops = pp.edgeLoops();
    const labelList& mp = pp.meshPoints();

    for (const labelList& loop : edgeLoops)
    {
        if (loop.size() > 2)
        {
            compactFaces.append(face(loop.size()));
            face& f = compactFaces.last();

            label fpi = 0;
            forAll(f, i)
            {
                const label pointi = mp[loop[i]];
                if (filterDiag && pointFromDiag[pointi])
                {
                    const label prevPointi = mp[loop[loop.fcIndex(i)]];
                    if
                    (
                        pointFromDiag[prevPointi]
                     && (pointToFace[pointi] != pointToFace[prevPointi])
                    )
                    {
                        f[fpi++] = pointi;
                    }
                    else
                    {
                        // Filter out diagonal point
                    }
                }
                else
                {
                    f[fpi++] = pointi;
                }
            }

            if (fpi > 2)
            {
                f.resize(fpi);
            }
            else
            {
                // Keep original face
                forAll(f, i)
                {
                    const label pointi = mp[loop[i]];
                    f[i] = pointi;
                }
            }
            compactCellIDs.append(cellID);
        }
    }
}


void Foam::isoSurfaceTopo::removeInsidePoints
(
    Mesh& s,
    const bool filterDiag,

    // inputs
    const boolUList& pointFromDiag,
    const labelUList& pointToFace,
    const labelUList& start,                // Per cell the starting triangle

    // outputs
    DynamicList<label>& pointCompactMap,    // Per returned point the original
    DynamicList<label>& compactCellIDs      // Per returned tri the cellID
)
{
    pointCompactMap.clear();
    compactCellIDs.clear();

    const pointField& points = s.points();

    DynamicList<face> compactFaces(s.size()/8);

    for (label celli = 0; celli < start.size()-1; ++celli)
    {
        // All triangles for the current cell

        const label nTris = start[celli+1]-start[celli];

        if (nTris)
        {
            const primitivePatch pp
            (
                SubList<face>(s, nTris, start[celli]),
                points
            );

            triangulateOutside
            (
                filterDiag,
                pp,
                pointFromDiag,
                pointToFace,
                //protectedFace,
                celli,

                compactFaces,
                compactCellIDs
            );
        }
    }


    // Compact out unused points
    labelList oldToCompact(points.size(), -1);
    pointCompactMap.clear();  // Extra safety (paranoid)

    for (face& f : compactFaces)
    {
        forAll(f, fp)
        {
            label pointi = f[fp];
            label compacti = oldToCompact[pointi];
            if (compacti == -1)
            {
                compacti = pointCompactMap.size();
                oldToCompact[pointi] = compacti;
                pointCompactMap.append(pointi);
            }
            f[fp] = compacti;
        }
    }

    pointField compactPoints
    (
        UIndirectList<point>(s.points(), pointCompactMap)
    );

    surfZoneList newZones(s.surfZones());

    s.clear();
    Mesh updated
    (
        std::move(compactPoints),
        std::move(compactFaces),
        s.surfZones()
    );

    s.transfer(updated);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceTopo::isoSurfaceTopo
(
    const polyMesh& mesh,
    const scalarField& cellValues,
    const scalarField& pointValues,
    const scalar iso,
    const isoSurfaceParams& params,
    const bitSet& ignoreCells
)
:
    isoSurfaceBase(mesh, cellValues, pointValues, iso, params),
    cellCutType_(mesh.nCells(), cutType::UNVISITED)
{
    if (debug)
    {
        Pout<< "isoSurfaceTopo:" << nl
            << "    cell min/max  : " << minMax(cVals_) << nl
            << "    point min/max : " << minMax(pVals_) << nl
            << "    isoValue      : " << iso << nl
            << "    filter        : "
            << isoSurfaceParams::filterNames[params.filter()] << nl
            << "    mesh span     : " << mesh.bounds().mag() << nl
            << "    ignoreCells   : " << ignoreCells.count()
            << " / " << cVals_.size() << nl
            << endl;
    }

    this->ignoreCyclics();

    label nBlockedCells = 0;

    // Mark ignoreCells as blocked
    nBlockedCells += blockCells(cellCutType_, ignoreCells);

    // Mark cells outside bounding box as blocked
    nBlockedCells +=
        blockCells(cellCutType_, params.getClipBounds(), volumeType::OUTSIDE);


    fixTetBasePtIs();

    // Determine cell cuts
    const label nCutCells = calcCellCuts(cellCutType_);

    if (debug)
    {
        Pout<< "isoSurfaceTopo : candidate cells cut "
            << nCutCells
            << " blocked " << nBlockedCells
            << " total " << mesh_.nCells() << endl;
    }

    if (debug && isA<fvMesh>(mesh))
    {
        const auto& fvmesh = dynamicCast<const fvMesh>(mesh);

        volScalarField debugField
        (
            IOobject
            (
                "isoSurfaceTopo.cutType",
                fvmesh.time().timeName(),
                fvmesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvmesh,
            dimensionedScalar(dimless, Zero)
        );

        auto& debugFld = debugField.primitiveFieldRef();

        forAll(cellCutType_, celli)
        {
            debugFld[celli] = cellCutType_[celli];
        }

        Pout<< "Writing cut types:"
            << debugField.objectPath() << endl;

        debugField.write();
    }


    // Per cell: 5 pyramids cut, each generating 2 triangles
    //  - pointToVerts : from generated iso point to originating mesh verts
    DynamicList<edge> pointToVerts(10*nCutCells);
    //  - pointToFace : from generated iso point to originating mesh face
    DynamicList<label> pointToFace(10*nCutCells);
    //  - pointFromDiag : if generated iso point is on face diagonal
    DynamicList<bool> pointFromDiag(10*nCutCells);

    // Per cell: number of intersected edges:
    //          - four faces cut so 4 mesh edges + 4 face-diagonal edges
    //          - 4 of the pyramid edges
    EdgeMap<label> vertsToPoint(12*nCutCells);
    DynamicList<label> verts(12*nCutCells);
    // Per cell: 5 pyramids cut (since only one pyramid not cut)
    DynamicList<label> faceLabels(5*nCutCells);
    DynamicList<label> cellLabels(5*nCutCells);


    labelList startTri(mesh_.nCells()+1, Zero);

    for (label celli = 0; celli < mesh_.nCells(); ++celli)
    {
        startTri[celli] = faceLabels.size();
        if ((cellCutType_[celli] & cutType::ANYCUT) != 0)
        {
            generateTriPoints
            (
                celli,

                // Same as tetMatcher::test(mesh_, celli),
                bool(cellCutType_[celli] & cutType::TETCUT),  // isTet

                pointToVerts,
                pointToFace,
                pointFromDiag,

                vertsToPoint,
                verts,
                faceLabels
            );

            for (label i = startTri[celli]; i < faceLabels.size(); ++i)
            {
                cellLabels.append(celli);
            }
        }
    }
    startTri[mesh_.nCells()] = faceLabels.size();


    pointToVerts_.transfer(pointToVerts);
    meshCells_.transfer(cellLabels);
    pointToFace_.transfer(pointToFace);

    pointField allPoints
    (
        this->interpolateTemplate
        (
            mesh_.cellCentres(),
            mesh_.points()
        )
    );


    // Assign to MeshedSurface
    faceList allTris(faceLabels.size());
    label verti = 0;
    for (face& allTri : allTris)
    {
        allTri.resize(3);
        allTri[0] = verts[verti++];
        allTri[1] = verts[verti++];
        allTri[2] = verts[verti++];
    }


    surfZoneList allZones(one{}, surfZone("allFaces", allTris.size()));

    Mesh::clear();
    Mesh updated
    (
        std::move(allPoints),
        std::move(allTris),
        std::move(allZones)
    );
    Mesh::transfer(updated);

    // Now:
    // - generated faces and points are assigned to *this
    // - per point we know:
    //  - pointOnDiag: whether it is on a face-diagonal edge
    //  - pointToFace_: from what pyramid (cell+face) it was produced
    //    (note that the pyramid faces are shared between multiple mesh faces)
    //  - pointToVerts_ : originating mesh vertex or cell centre


    if (debug)
    {
        Pout<< "isoSurfaceTopo : generated "
            << Mesh::size() << " faces "
            << Mesh::points().size() << " points" << endl;
    }


    if (params.filter() != filterType::NONE)
    {
        // Triangulate outside (filter edges to cell centres and optionally
        // face diagonals)
        DynamicList<label> pointCompactMap(size()); // Back to original point
        DynamicList<label> compactCellIDs(size());  // Per tri the cell

        removeInsidePoints
        (
            *this,
            (params.filter() == filterType::DIAGCELL ? true : false),

            pointFromDiag,
            pointToFace_,
            startTri,
            pointCompactMap,
            compactCellIDs
        );

        pointToVerts_ = UIndirectList<edge>(pointToVerts_, pointCompactMap)();
        pointToFace_ = UIndirectList<label>(pointToFace_, pointCompactMap)();
        meshCells_.transfer(compactCellIDs);

        pointCompactMap.clearStorage();
        compactCellIDs.clearStorage();

        if (debug)
        {
            Pout<< "isoSurfaceTopo :"
                " after removing cell centre and face-diag triangles "
                << Mesh::size() << " faces "
                << Mesh::points().size() << " points"
                << endl;
        }
    }

    // Not required after this stage
    pointFromDiag.clearStorage();


    if (params.filter() == filterType::DIAGCELL)
    {
        // We remove verts on face diagonals. This is in fact just
        // straightening the edges of the face through the cell. This can
        // close off 'pockets' of triangles and create open or
        // multiply-connected triangles

        // Solved by eroding open-edges
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Mark points on mesh outside.
        bitSet isBoundaryPoint(mesh.nPoints());
        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            ++facei
        )
        {
            isBoundaryPoint.set(mesh.faces()[facei]);
        }


        // Include faces that would be exposed from mesh subset
        if (nBlockedCells)
        {
            const labelList& faceOwn = mesh_.faceOwner();
            const labelList& faceNei = mesh_.faceNeighbour();

            for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
            {
                // If only one cell is blocked, the face corresponds
                // to an exposed subMesh face
                if
                (
                    (cellCutType_[faceOwn[facei]] == cutType::BLOCKED)
                 != (cellCutType_[faceNei[facei]] == cutType::BLOCKED)
                )
                {
                    isBoundaryPoint.set(mesh.faces()[facei]);
                }
            }
        }


        // The list of surface faces that should be retained after erosion
        Mesh& surf = *this;
        labelList faceAddr(identity(surf.size()));

        bitSet faceSelection;

        while (true)
        {
            // Shadow the surface for the purposes of erosion
            uindirectPrimitivePatch erosion
            (
                UIndirectList<face>(surf, faceAddr),
                surf.points()
            );

            faceSelection.clear();
            faceSelection.resize(erosion.size());

            const labelList& mp = erosion.meshPoints();
            const edgeList& surfEdges = erosion.edges();
            const labelListList& edgeFaces = erosion.edgeFaces();

            label nEdgeRemove = 0;

            forAll(edgeFaces, edgei)
            {
                const labelList& eFaces = edgeFaces[edgei];
                if (eFaces.size() == 1)
                {
                    // Open edge. Check that vertices do not originate
                    // from a boundary face
                    const edge& e = surfEdges[edgei];

                    const edge& verts0 = pointToVerts_[mp[e.first()]];
                    const edge& verts1 = pointToVerts_[mp[e.second()]];

                    if
                    (
                        isBoundaryPoint.test(verts0.first())
                     && isBoundaryPoint.test(verts0.second())
                     && isBoundaryPoint.test(verts1.first())
                     && isBoundaryPoint.test(verts1.second())
                    )
                    {
                        // Open edge on boundary face. Keep
                    }
                    else
                    {
                        // Open edge. Mark for erosion
                        faceSelection.set(eFaces[0]);
                        ++nEdgeRemove;
                    }
                }
            }

            if (debug)
            {
                Pout<< "isoSurfaceTopo :"
                    << " removing " << faceSelection.count()
                    << " / " << faceSelection.size()
                    << " faces on " << nEdgeRemove << " open edges" << endl;
            }

            if (returnReduce(faceSelection.none(), andOp<bool>()))
            {
                break;
            }

            // Remove the faces from the addressing
            inplaceSubset(faceSelection, faceAddr, true);  // True = remove
        }


        // Finished erosion (if any)
        // - retain the faces listed in the updated addressing

        if (surf.size() != faceAddr.size())
        {
            faceSelection.clear();
            faceSelection.resize(surf.size());
            faceSelection.set(faceAddr);

            inplaceSubsetMesh(faceSelection);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoSurfaceTopo::inplaceSubsetMesh(const bitSet& include)
{
    labelList pointMap;
    labelList faceMap;
    Mesh filtered
    (
        Mesh::subsetMesh(include, pointMap, faceMap)
    );
    Mesh::transfer(filtered);

    meshCells_ = UIndirectList<label>(meshCells_, faceMap)();

    pointToVerts_ = UIndirectList<edge>(pointToVerts_, pointMap)();
    pointToFace_ = UIndirectList<label>(pointToFace_, pointMap)();
}


// ************************************************************************* //
