/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Description
    Collapse short edges and combines edges that are in line.

    - collapse short edges. Length of edges to collapse provided as argument.
    - merge two edges if they are in line. Maximum angle provided as argument.
    - remove unused points.

    Optionally removes cells. Can remove faces and points but does not check
    for nonsense resulting topology.

    When collapsing an edge with one point on the boundary it will leave
    the boundary point intact. When both points inside it chooses random. When
    both points on boundary random again.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "fvMesh.H"
#include "mapPolyMesh.H"
#include "mathematicalConstants.H"
#include "PackedBoolList.H"
#include "unitConversion.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "PointEdgeWave.H"
#include "pointEdgeCollapse.H"
#include "motionSmoother.H"
#include "removePoints.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void filterFace
(
    const label faceI,
    face& f,
    const List<pointEdgeCollapse>& allPointInfo,
    const Map<DynamicList<label> >& collapseStrings
)
{
    label newFp = 0;

    face oldFace = f;

    forAll(f, fp)
    {
        label pointI = f[fp];

        label collapsePoint = allPointInfo[pointI].collapseIndex();

        if (collapseStrings.found(collapsePoint))
        {
            collapsePoint = collapseStrings[collapsePoint][0];
        }

        if (collapsePoint == -1)
        {
            WarningIn
                (
                    "filterFace"
                    "(const label, face&, const List<pointEdgeCollapse>&)"
                )
                << "Point " << pointI << " was not visited by PointEdgeWave"
                << endl;
        }
        else if (collapsePoint == -2)
        {
            f[newFp++] = pointI;
        }
        else
        {
            if (findIndex(SubList<label>(f, newFp), collapsePoint) == -1)
            {
                f[newFp++] = collapsePoint;
            }
        }
    }


    // Check for pinched face. Tries to correct
    // - consecutive duplicate vertex. Removes duplicate vertex.
    // - duplicate vertex with one other vertex in between (spike).
    // Both of these should not really occur! and should be checked before
    // collapsing edges.

    const label size = newFp;

    newFp = 2;

    for (label fp = 2; fp < size; fp++)
    {
        label fp1 = fp-1;
        label fp2 = fp-2;

        label pointI = f[fp];

        // Search for previous occurrence.
        label index = findIndex(SubList<label>(f, fp), pointI);

        if (index == fp1)
        {
            WarningIn
            (
                "Foam::edgeCollapser::filterFace(const label faceI, "
                "face& f) const"
            )   << "Removing consecutive duplicate vertex in face "
                << f << endl;
            // Don't store current pointI
        }
        else if (index == fp2)
        {
            WarningIn
            (
                "Foam::edgeCollapser::filterFace(const label faceI, "
                "face& f) const"
            )   << "Removing non-consecutive duplicate vertex in face "
                << f << endl;
            // Don't store current pointI and remove previous
            newFp--;
        }
        else if (index != -1)
        {
            WarningIn
            (
                "Foam::edgeCollapser::filterFace(const label faceI, "
                "face& f) const"
            )   << "Pinched face " << f << endl;
            f[newFp++] = pointI;
        }
        else
        {
            f[newFp++] = pointI;
        }
    }

    f.setSize(newFp);
}


bool setRefinement
(
    const polyMesh& mesh,
    polyTopoChange& meshMod,
    const List<pointEdgeCollapse>& allPointInfo
)
{
    const cellList& cells = mesh.cells();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const labelListList& pointFaces = mesh.pointFaces();
    const pointZoneMesh& pointZones = mesh.pointZones();

    globalIndex globalStrings(mesh.nPoints());

    boolList removedPoints(mesh.nPoints(), false);

    // Create strings of edges.
    // Map from collapseIndex(=global master point) to set of points
    Map<DynamicList<label> > collapseStrings;
    {
        // 1. Count elements per collapseIndex
        Map<label> nPerIndex(mesh.nPoints()/10);
        forAll(allPointInfo, pointI)
        {
            const label collapseIndex = allPointInfo[pointI].collapseIndex();

            if (collapseIndex != -1 && collapseIndex != -2)
            {
                Map<label>::iterator fnd = nPerIndex.find(collapseIndex);
                if (fnd != nPerIndex.end())
                {
                    fnd()++;
                }
                else
                {
                    nPerIndex.insert(collapseIndex, 1);
                }
            }
        }

        // 2. Size
        collapseStrings.resize(2*nPerIndex.size());
        forAllConstIter(Map<label>, nPerIndex, iter)
        {
            collapseStrings.insert(iter.key(), DynamicList<label>(iter()));
        }

        // 3. Fill
        forAll(allPointInfo, pointI)
        {
            const label collapseIndex = allPointInfo[pointI].collapseIndex();

            if (collapseIndex != -1 && collapseIndex != -2)
            {
                collapseStrings[collapseIndex].append(pointI);
            }
        }
    }


    bool meshChanged = false;

    // Current faces (is also collapseStatus: f.size() < 3)
    faceList newFaces(mesh.faces());

    // Current cellCollapse status
    boolList cellRemoved(mesh.nCells(), false);

    label nUnvisited = 0;
    label nUncollapsed = 0;
    label nCollapsed = 0;

    forAll(allPointInfo, pI)
    {
        const pointEdgeCollapse& pec = allPointInfo[pI];

        if (pec.collapseIndex() == -1)
        {
            nUnvisited++;
        }
        else if (pec.collapseIndex() == -2)
        {
            nUncollapsed++;
        }
        else if (pec.collapseIndex() >= 0)
        {
            nCollapsed++;
        }
    }

    label nPoints = allPointInfo.size();

    reduce(nPoints, sumOp<label>());
    reduce(nUnvisited, sumOp<label>());
    reduce(nUncollapsed, sumOp<label>());
    reduce(nCollapsed, sumOp<label>());

    Info<< incrIndent;
    Info<< indent << "Number of points : " << nPoints << nl
        << indent << "Not visited      : " << nUnvisited << nl
        << indent << "Not collapsed    : " << nUncollapsed << nl
        << indent << "Collapsed        : " << nCollapsed << nl
        << endl;
    Info<< decrIndent;

    do
    {
        // Update face collapse from edge collapses
        forAll(newFaces, faceI)
        {
            filterFace(faceI, newFaces[faceI], allPointInfo, collapseStrings);
        }

        // Check if faces to be collapsed cause cells to become collapsed.
        label nCellCollapsed = 0;

        forAll(cells, cellI)
        {
            if (!cellRemoved[cellI])
            {
                const cell& cFaces = cells[cellI];

                label nFaces = cFaces.size();

                forAll(cFaces, i)
                {
                    label faceI = cFaces[i];

                    if (newFaces[faceI].size() < 3)
                    {
                        --nFaces;

                        if (nFaces < 4)
                        {
                            Pout<< "Cell:" << cellI
                                << " uses faces:" << cFaces
                                << " of which too many are marked for removal:"
                                << endl
                                << "   ";
                            forAll(cFaces, j)
                            {
                                if (newFaces[cFaces[j]].size() < 3)
                                {
                                    Pout<< ' '<< cFaces[j];
                                }
                            }
                            Pout<< endl;

                            cellRemoved[cellI] = true;

                            // Collapse all edges of cell to nothing
                            //collapseEdges(cellEdges[cellI]);

                            nCellCollapsed++;

                            break;
                        }
                    }
                }
            }
        }

        if (nCellCollapsed == 0)
        {
            break;
        }
    } while (true);


    // Keep track of faces that have been done already.
    boolList doneFace(mesh.nFaces(), false);

    {
        // Mark points used.
        boolList usedPoint(mesh.nPoints(), false);

        forAll(cellRemoved, cellI)
        {
            if (cellRemoved[cellI])
            {
                meshMod.removeCell(cellI, -1);
            }
        }

        // Remove faces
        forAll(newFaces, faceI)
        {
            const face& f = newFaces[faceI];

            if (f.size() < 3)
            {
                meshMod.removeFace(faceI, -1);
                meshChanged = true;

                // Mark face as been done.
                doneFace[faceI] = true;
            }
            else
            {
                // Kept face. Mark vertices
                forAll(f, fp)
                {
                    usedPoint[f[fp]] = true;
                }
            }
        }

        // Remove unused vertices that have not been marked for removal already
        forAll(usedPoint, pointI)
        {
            if (!usedPoint[pointI])
            {
                removedPoints[pointI] = true;
                meshMod.removePoint(pointI, -1);
                meshChanged = true;
            }
        }
    }

    // Modify the point location of the remaining points
    forAll(allPointInfo, pointI)
    {
        const label collapseIndex = allPointInfo[pointI].collapseIndex();
        const point& collapsePoint = allPointInfo[pointI].collapsePoint();

        if
        (
            removedPoints[pointI] == false
         && collapseIndex != -1
         && collapseIndex != -2
        )
        {
            meshMod.modifyPoint
            (
                pointI,
                collapsePoint,
                pointZones.whichZone(pointI),
                false
            );
        }
    }


    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const faceZoneMesh& faceZones = mesh.faceZones();

    // Renumber faces that use points
    forAll(allPointInfo, pointI)
    {
        if (removedPoints[pointI] == true)
        {
            const labelList& changedFaces = pointFaces[pointI];

            forAll(changedFaces, changedFaceI)
            {
                label faceI = changedFaces[changedFaceI];

                if (!doneFace[faceI])
                {
                    doneFace[faceI] = true;

                    // Get current zone info
                    label zoneID = faceZones.whichZone(faceI);

                    bool zoneFlip = false;

                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = faceZones[zoneID];

                        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                    }

                    // Get current connectivity
                    label own = faceOwner[faceI];
                    label nei = -1;
                    label patchID = -1;

                    if (mesh.isInternalFace(faceI))
                    {
                        nei = faceNeighbour[faceI];
                    }
                    else
                    {
                        patchID = boundaryMesh.whichPatch(faceI);
                    }

                    meshMod.modifyFace
                    (
                        newFaces[faceI],            // face
                        faceI,                      // faceI to change
                        own,                        // owner
                        nei,                        // neighbour
                        false,                      // flipFaceFlux
                        patchID,                    // patch
                        zoneID,
                        zoneFlip
                    );
                    meshChanged = true;
                }
            }
        }
    }

    return meshChanged;
}


// Get faceEdges in order of face points, i.e. faceEdges[0] is between
// f[0] and f[1]
labelList getSortedEdges
(
    const edgeList& edges,
    const labelList& f,
    const labelList& edgeLabels
)
{
    labelList faceEdges(edgeLabels.size(), -1);

    // Find starting pos in f for every edgeLabels
    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];

        const edge& e = edges[edgeI];

        label fp = findIndex(f, e[0]);
        label fp1 = f.fcIndex(fp);

        if (f[fp1] == e[1])
        {
            // EdgeI between fp -> fp1
            faceEdges[fp] = edgeI;
        }
        else
        {
            // EdgeI between fp-1 -> fp
            faceEdges[f.rcIndex(fp)] = edgeI;
        }
    }

    return faceEdges;
}


// Return master point edge needs to be collapsed to (or -1)
label edgeMaster
(
    const labelList& boundaryPoint,
    const bool flipEdge,
    const edge& e
)
{
    label masterPoint = -1;

    label e0 = e[0];
    label e1 = e[1];

    if (flipEdge)
    {
        e0 = e[1];
        e1 = e[0];
    }

    // Check if one of the points is on a processor
//    if
//    (
//        boundaryPoint[e0] > 0
//     && boundaryPoint[e1] > 0
//    )
//    {
//        if (boundaryPoint[e0] != boundaryPoint[e1])
//        {
//            return -1;
//        }
//    }
//
//    if (boundaryPoint[e0] > 0)
//    {
//        return e0;
//    }
//    else if (boundaryPoint[e1] > 0)
//    {
//        return e1;
//    }

    // Collapse edge to boundary point.
    if (boundaryPoint[e0] == 0)
    {
        if (boundaryPoint[e1] == 0)
        {
            // Both points on boundary. Choose one to collapse to.
            // Note: should look at feature edges/points!
            masterPoint = e0;
        }
        else
        {
            masterPoint = e0;
        }
    }
    else
    {
        if (boundaryPoint[e1] == 0)
        {
            masterPoint = e1;
        }
        else
        {
            // None on boundary. Choose arbitrary.
            // Note: should look at geometry?
            masterPoint = e0;
        }
    }

    return masterPoint;
}


// Create consistent set of collapses.
//  collapseEdge : per edge:
//      -1 : do not collapse
//       0 : collapse to start
//       1 : collapse to end
//  Note: collapseEdge has to be parallel consistent (in orientation)
label syncCollapse
(
    const polyMesh& mesh,
    const globalIndex& globalStrings,
    const labelList& collapseEdge,
    List<pointEdgeCollapse>& allPointInfo
)
{
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    label nCollapsed = 0;

    DynamicList<label> initPoints(mesh.nPoints());
    DynamicList<pointEdgeCollapse> initPointInfo(mesh.nPoints());
    allPointInfo.resize(mesh.nPoints());

    // Initialise edges to no collapse
    List<pointEdgeCollapse> allEdgeInfo
    (
        mesh.nEdges(),
        pointEdgeCollapse(vector::zero, -1)
    );

    // Mark selected edges for collapse
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (collapseEdge[edgeI] != -1)
        {
            label masterPointI = e[collapseEdge[edgeI]];

            const pointEdgeCollapse pec
            (
                points[masterPointI],
                globalStrings.toGlobal(masterPointI)
            );

            // Mark as collapsable but with nonsense master so it gets
            // overwritten and starts an update wave
            allEdgeInfo[edgeI] = pointEdgeCollapse
            (
                points[masterPointI],
                labelMax
            );

            initPointInfo.append(pec);
            initPoints.append(e.start());

            initPointInfo.append(pec);
            initPoints.append(e.end());

            nCollapsed++;
        }
    }

    PointEdgeWave<pointEdgeCollapse> collapsePropagator
    (
        mesh,
        initPoints,
        initPointInfo,
        allPointInfo,
        allEdgeInfo,
        mesh.globalData().nTotalPoints()  // Maximum number of iterations
    );

    return nCollapsed;
}


void syncCollapseEdge(const polyMesh& mesh, labelList& collapseEdge)
{
    // Check whether edge point order is reversed from mesh to coupledPatch
    const globalMeshData& globalData = mesh.globalData();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();
    const labelList& coupledMeshEdges = globalData.coupledPatchMeshEdges();
    const indirectPrimitivePatch& coupledPatch = globalData.coupledPatch();
    const PackedBoolList& cppOrientation = globalData.globalEdgeOrientation();
    PackedBoolList meshToPatchSameOrientation(coupledMeshEdges.size(), true);

    forAll(coupledMeshEdges, eI)
    {
        const label meshEdgeIndex = coupledMeshEdges[eI];

        if (collapseEdge[meshEdgeIndex] != -1)
        {
            const edge& meshEdge = mesh.edges()[meshEdgeIndex];
            const edge& coupledPatchEdge = coupledPatch.edges()[eI];

            if
            (
                meshEdge[0] == coupledPatch.meshPoints()[coupledPatchEdge[1]]
             && meshEdge[1] == coupledPatch.meshPoints()[coupledPatchEdge[0]]
            )
            {
                meshToPatchSameOrientation[eI] = false;
            }
        }
    }


    labelList cppEdgeData(map.constructSize());

    forAll(coupledMeshEdges, eI)
    {
        const label meshEdgeIndex = coupledMeshEdges[eI];

        cppEdgeData[eI] = collapseEdge[meshEdgeIndex];

        if
        (
            (collapseEdge[meshEdgeIndex] != -1)
         && (meshToPatchSameOrientation[eI] != cppOrientation[eI])
        )
        {
            cppEdgeData[eI] = 1-cppEdgeData[eI];
        }
    }


    // Synchronise cppEdgeData
    // Use minEqOp reduction, so that edge will only be collapsed on processor
    // boundary if both processors agree to collapse it
    globalData.syncData
    (
        cppEdgeData,
        globalData.globalEdgeSlaves(),
        globalData.globalEdgeTransformedSlaves(),
        map,
        minEqOp<label>()
    );


    forAll(coupledMeshEdges, eI)
    {
        const label meshEdgeIndex = coupledMeshEdges[eI];

        collapseEdge[meshEdgeIndex] = cppEdgeData[eI];

        if
        (
            (cppEdgeData[eI] != -1)
         && (meshToPatchSameOrientation[eI] != cppOrientation[eI])
        )
        {
            collapseEdge[meshEdgeIndex] = 1-collapseEdge[meshEdgeIndex];
        }
    }
}


// Mark (in collapseEdge) any edges to collapse
label collapseSmallEdges
(
    const polyMesh& mesh,
    const labelList& boundaryPoint,
    const scalarField& minEdgeLen,
    labelList& collapseEdge
)
{
    // Work out which edges to collapse
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nCollapsed = 0;

    forAll(mesh.edges(), edgeI)
    {
        if (collapseEdge[edgeI] == -1)
        {
            const edge& e = mesh.edges()[edgeI];

            if (e.mag(mesh.points()) < minEdgeLen[edgeI])
            {
                label masterPointI = edgeMaster(boundaryPoint, false, e);
                if (masterPointI == e[0])
                {
                    collapseEdge[edgeI] = 0;
                }
                else
                {
                    collapseEdge[edgeI] = 1;
                }
                nCollapsed++;
            }
        }
    }
    return nCollapsed;
}


// Mark (in collapseEdge) any edges to merge
label mergeEdges
(
    const polyMesh& mesh,
    const scalar maxCos,
    const labelList& boundaryPoint,
    const scalarField& minEdgeLen,
    labelList& collapseEdge
)
{
    const edgeList& edges = mesh.edges();
    const labelListList& pointEdges = mesh.pointEdges();

    // Point removal engine
    removePoints pointRemover(mesh, false);

    // Find out points that can be deleted
    boolList pointCanBeDeleted;
    label nTotRemove = pointRemover.countPointUsage(maxCos, pointCanBeDeleted);


    // Rework point-to-remove into edge-to-collapse.

    label nCollapsed = 0;

    if (nTotRemove > 0)
    {
        forAll(pointEdges, pointI)
        {
            if (pointCanBeDeleted[pointI])
            {
                const labelList& pEdges = pointEdges[pointI];

                if (pEdges.size() == 2)
                {
                    // Always the case?

                    label e0 = pEdges[0];
                    label e1 = pEdges[1];

                    if
                    (
                        collapseEdge[e0] == -1
                     && minEdgeLen[e0] >= 0
                     && collapseEdge[e1] == -1
                     && minEdgeLen[e1] >= 0
                    )
                    {
                        // Get the two vertices on both sides of the point
                        label leftV = edges[e0].otherVertex(pointI);
                        label rightV = edges[e1].otherVertex(pointI);

                        // Can collapse pointI onto either leftV or rightV.
                        // Preferentially choose an internal point to hopefully
                        // give less distortion

                        if (boundaryPoint[leftV] == -1)
                        {
                            collapseEdge[e0] = findIndex(edges[e0], leftV);
                        }
                        else
                        {
                            collapseEdge[e1] = findIndex(edges[e1], rightV);
                        }
                    }
                }
            }
        }
    }
    return nCollapsed;
}


// Make consistent set of collapses that does not collapse any cells
label consistentCollapse
(
    const bool allowCellCollapse,
    const polyMesh& mesh,
    const globalIndex& globalStrings,
    labelList& collapseEdge,
    List<pointEdgeCollapse>& allPointInfo
)
{
    // Make sure we don't collapse cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    while (true)
    {
        // Sync collapseEdge
        syncCollapseEdge(mesh, collapseEdge);


        // Get collapsed faces

        label nAdditionalCollapsed = 0;

        PackedBoolList isCollapsedFace(mesh.nFaces());
        forAll(mesh.faceEdges(), faceI)
        {
            const labelList& fEdges = mesh.faceEdges()[faceI];

            // Count number of remaining edges
            label nEdges = fEdges.size();
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];
                if (collapseEdge[edgeI] != -1)
                {
                    nEdges--;
                }
            }

            if (nEdges < 3)
            {
                // Face is collapsed.
                isCollapsedFace[faceI] = 1;

                if (nEdges == 1)
                {
                    // Cannot collapse face down to single edge.

                    //- Option 1: collapse remaining edge as well. However
                    //  if this edge is on the coupled processor patch this
                    //  logic clashes with that of syncCollapseEdge
                    //  (do not collapse if any not collapse)
                    //forAll(fEdges, fEdgeI)
                    //{
                    //    label edgeI = fEdges[fEdgeI];
                    //    if (collapseEdge[edgeI] == -1)
                    //    {
                    //        collapseEdge[edgeI] = 0;
                    //        nAdditionalCollapsed++;
                    //    }
                    //}

                    //- Option 2: uncollapse this face.
                    forAll(fEdges, fEdgeI)
                    {
                        label edgeI = fEdges[fEdgeI];
                        if (collapseEdge[edgeI] != -1)
                        {
                            collapseEdge[edgeI] = -1;
                            nAdditionalCollapsed++;
                        }
                    }
                }
            }
        }

        //Pout<< "nAdditionalCollapsed : " << nAdditionalCollapsed << endl;


        label nUncollapsed = 0;

        if (!allowCellCollapse)
        {
            // Check collapsed cells

            forAll(mesh.cells(), cellI)
            {
                const cell& cFaces = mesh.cells()[cellI];
                label nFaces = cFaces.size();
                forAll(cFaces, i)
                {
                    label faceI = cFaces[i];
                    if (isCollapsedFace[faceI])
                    {
                        nFaces--;
                        if (nFaces < 4)
                        {
                            // Unmark this face for collapse
                            const labelList& fEdges = mesh.faceEdges()[faceI];

                            forAll(fEdges, fEdgeI)
                            {
                                label edgeI = fEdges[fEdgeI];
                                if (collapseEdge[edgeI] != -1)
                                {
                                    collapseEdge[edgeI] = -1;
                                    nUncollapsed++;
                                }
                            }

                            // Uncollapsed this face.
                            isCollapsedFace[faceI] = 0;
                            nFaces++;
                        }
                    }
                }
            }
            //Pout<< "** disallowing cells : nUncollapsed : "
            //    << nUncollapsed << endl;
        }


        if
        (
            returnReduce(nUncollapsed+nAdditionalCollapsed, sumOp<label>())
            == 0
        )
        {
            break;
        }
    }


    // Create consistent set of collapses
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: requires collapseEdge to be synchronised. (above loop makes sure
    //       of that)

    return syncCollapse(mesh, globalStrings, collapseEdge, allPointInfo);
}


// Check mesh and mark points on faces in error
// Returns boolList with points in error set
PackedBoolList checkMeshQuality(const polyMesh& mesh)
{
    //mesh.checkMesh(true);

    IOdictionary meshQualityDict
    (
        IOobject
        (
            "meshQualityDict",
            mesh.time().system(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    labelHashSet badFaces(mesh.nFaces()/100);
    DynamicList<label> checkFaces(mesh.nFaces());

    const vectorField& fAreas = mesh.faceAreas();

    scalar faceAreaLimit = SMALL;

    forAll(fAreas, fI)
    {
        if (mag(fAreas[fI]) > faceAreaLimit)
        {
            checkFaces.append(fI);
        }
    }

    motionSmoother::checkMesh
    (
        false,
        mesh,
        meshQualityDict,
        checkFaces,
        badFaces
    );

    label nBadFaces = badFaces.size();
    reduce(nBadFaces, sumOp<label>());

    Info<< nl << "Number of bad faces          : " << nBadFaces << endl;

    PackedBoolList isErrorPoint(mesh.nPoints());
    forAllConstIter(labelHashSet, badFaces, iter)
    {
        const face& f = mesh.faces()[iter.key()];

        forAll(f, pI)
        {
            isErrorPoint[f[pI]] = 1;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        isErrorPoint,
        orEqOp<unsigned int>(),
        0
    );

    return isErrorPoint;
}


// Check mesh with collapses (newMesh), updates minEdgeLen, nFrozenEdges
void checkMeshAndFreezeEdges
(
    const polyMesh& newMesh,
    const labelList& oldToNewMesh,
    const polyMesh& oldMesh,
    scalarField& minEdgeLen,
    label& nFrozenEdges
)
{
    PackedBoolList isErrorPoint = checkMeshQuality(newMesh);

    forAll(oldMesh.edges(), edgeI)
    {
        const edge& e = oldMesh.edges()[edgeI];
        label newStart = oldToNewMesh[e[0]];
        label newEnd = oldToNewMesh[e[1]];

        if
        (
            (newStart >= 0 && isErrorPoint[newStart])
         || (newEnd >= 0 && isErrorPoint[newEnd])
        )
        {
            // Gradual decrease? For now just hard disable
            if (minEdgeLen[edgeI] > -GREAT/2)
            {
                minEdgeLen[edgeI] = -GREAT;
                nFrozenEdges++;
            }
        }
    }
}


// Mark boundary points
// boundaryPoint:
// + -1 : point not on boundary
// +  0 : point on a real boundary
// + >0 : point on a processor patch with that ID
labelList findBoundaryPoints(const polyMesh& mesh)
{
    const faceList& faces = mesh.faces();
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();


    labelList boundaryPoint(mesh.nPoints(), -1);

    // Get all points on a boundary
    label nIntFaces = mesh.nInternalFaces();
    for (label faceI = nIntFaces; faceI < mesh.nFaces(); faceI++)
    {
        const face& f = faces[faceI];

        forAll(f, fp)
        {
           boundaryPoint[f[fp]] = 0;
        }
    }

    // Get all processor boundary points and the processor patch label
    // that they are on.
    forAll(bMesh, patchI)
    {
        const polyPatch& patch = bMesh[patchI];

        if (isA<coupledPolyPatch>(patch))
        {
            if (patchI == 0)
            {
                // We mark 'normal' boundary points with 0 so make sure this
                // coupled patch is not 0.
                FatalErrorIn("findBoundaryPoints(const polyMesh&)")
                    << "Your patches should have non-coupled ones before any"
                    << " coupled ones. Current patches " << bMesh.names()
                    << exit(FatalError);
            }

            forAll(patch, fI)
            {
                const face& f = patch[fI];

                forAll(f, fp)
                {
                    boundaryPoint[f[fp]] = patchI;
                }
            }
        }
    }

    return boundaryPoint;
}


// Main program:
int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Merges small and in-line edges.\n"
        "Collapses faces and optionally cells to a point."
    );

#   include "addOverwriteOption.H"

    argList::validArgs.append("edge length [m]");
    argList::validArgs.append("merge angle (degrees)");
    argList::addBoolOption
    (
        "allowCellCollapse",
        "Allow collapsing of cells to a single point"
    );
    argList::addBoolOption
    (
        "checkMeshQuality",
        "Only collapse if not exceeding given meshQualityDict limits"
    );



#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    scalar minLen = args.argRead<scalar>(1);
    const scalar angle = args.argRead<scalar>(2);

    const bool allowCellCollapse = args.optionFound("allowCellCollapse");
    const bool overwrite = args.optionFound("overwrite");
    const bool checkQuality = args.optionFound("checkMeshQuality");

    scalar maxCos = Foam::cos(degToRad(angle));

    Info<< "Merging:" << nl
        << "    edges with length less than " << minLen << " meters" << nl
        << "    edges split by a point with edges in line to within " << angle
        << " degrees" << nl
        << endl;

    if (allowCellCollapse)
    {
        Info<< "Allowing collapse of cells down to a point." << nl
            << endl;
    }
    else
    {
        Info<< "Disallowing collapse of cells down to a point." << nl
            << endl;
    }


    if (checkQuality)
    {
        Info<< "Selectively disabling wanted collapses until resulting quality"
            << " satisfies constraints in system/meshQualityDict" << nl
            << endl;
    }



    // To mark master of collapes
    globalIndex globalStrings(mesh.nPoints());


    // Local collapse length. Any edge below this length gets (attempted)
    // collapsed. Use either aa gradually decreasing value
    // (so from minLen to e.g. 0.5*minLen) or a hard limit (GREAT)
    scalarField minEdgeLen(mesh.nEdges(), minLen);
    label nFrozenEdges = 0;



    // Initial mesh check
    // ~~~~~~~~~~~~~~~~~~
    // Do not allow collapses in regions of error.
    // Updates minEdgeLen, nFrozenEdges
    if (checkQuality)
    {
        checkMeshAndFreezeEdges
        (
            mesh,
            identity(mesh.nPoints()),
            mesh,
            minEdgeLen,
            nFrozenEdges
        );
        Info<< "Initial frozen edges "
            << returnReduce(nFrozenEdges, sumOp<label>())
            << " out of " << returnReduce(mesh.nEdges(), sumOp<label>())
            << endl;
    }


    // Mark points on boundary
    const labelList boundaryPoint = findBoundaryPoints(mesh);

    // Copy of current set of topology changes. Used to generate final mesh.
    polyTopoChange savedMeshMod(mesh.boundaryMesh().size());

    // Keep track of whether mesh has changed at all
    bool meshChanged = false;



    // Main loop
    // ~~~~~~~~~
    // It tries and do some collapses, checks the resulting mesh and
    // 'freezes' some edges (by marking in minEdgeLen) and tries again.
    // This will iterate ultimately to the situation where every edge is
    // frozen and nothing gets collapsed.

    do
    {
        // Per edge collapse status:
        // -1 : not collapsed
        //  0 : collapse to start
        //  1 : collapse to end
        labelList collapseEdge(mesh.nEdges(), -1);


        // Work out which edges to collapse
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // This is by looking at minEdgeLen (to avoid frozen edges)
        // and marking in collapseEdge.

        // Small edges
        label nSmallCollapsed = collapseSmallEdges
        (
            mesh,
            boundaryPoint,
            minEdgeLen,
            collapseEdge
        );
        reduce(nSmallCollapsed, sumOp<label>());

        Info<< indent << "Collapsing " << nSmallCollapsed
            << " small edges" << endl;


        // Inline edges
        label nMerged = mergeEdges
        (
            mesh,
            maxCos,
            boundaryPoint,
            minEdgeLen,
            collapseEdge
        );

        reduce(nMerged, sumOp<label>());
        Info<< indent << "Collapsing " << nMerged << " in line edges"
            << endl;


        // Merge edge collapses into consistent collapse-network. Make sure
        // no cells get collapsed.
        List<pointEdgeCollapse> allPointInfo;
        label nLocalCollapse = consistentCollapse
        (
            allowCellCollapse,
            mesh,
            globalStrings,
            collapseEdge,
            allPointInfo
        );

        reduce(nLocalCollapse, sumOp<label>());
        Info<< "nLocalCollapse = " << nLocalCollapse << endl;

        if (nLocalCollapse == 0)
        {
            break;
        }


        // There are collapses so mesh will get changed
        meshChanged = true;


        // Apply collapses to current mesh
        polyTopoChange meshMod(mesh);

        // Insert mesh refinement into polyTopoChange.
        setRefinement(mesh, meshMod, allPointInfo);

        // Do all changes
        Info<< indent << "Applying changes to the mesh" << nl
            //<< decrIndent
            << endl;

        savedMeshMod = meshMod;

        autoPtr<fvMesh> newMeshPtr;
        autoPtr<mapPolyMesh> mapPtr = meshMod.makeMesh
        (
            newMeshPtr,
            IOobject
            (
                mesh.name(),
                mesh.instance(),
                mesh.time(),  // register with runTime
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            true        // parallel sync
        );

        fvMesh& newMesh = newMeshPtr();
        const mapPolyMesh& map = mapPtr();

        // Update fields
        newMesh.updateMesh(map);
        if (map.hasMotionPoints())
        {
            newMesh.movePoints(map.preMotionPoints());
        }




        // If no checks needed exit.
        if (!checkQuality)
        {
            break;
        }

        // Mesh check
        // ~~~~~~~~~~~~~~~~~~
        // Do not allow collapses in regions of error.
        // Updates minEdgeLen, nFrozenEdges
        label nOldFrozenEdges = returnReduce(nFrozenEdges, sumOp<label>());
        checkMeshAndFreezeEdges
        (
            newMesh,
            map.reversePointMap(),
            mesh,
            minEdgeLen,
            nFrozenEdges
        );
        label nNewFrozenEdges = returnReduce(nFrozenEdges, sumOp<label>());

        Info<< "Frozen edges "
            << returnReduce(nFrozenEdges, sumOp<label>())
            << " out of " << returnReduce(mesh.nEdges(), sumOp<label>())
            << endl;

        if (nNewFrozenEdges == nOldFrozenEdges)
        {
            break;
        }

    } while (true);


    if (meshChanged)
    {
        // Apply changes to current mesh
        autoPtr<mapPolyMesh> mapPtr = savedMeshMod.changeMesh(mesh, false);
        const mapPolyMesh& map = mapPtr();

        // Update fields
        mesh.updateMesh(map);
        if (map.hasMotionPoints())
        {
            mesh.movePoints(map.preMotionPoints());
        }


        // Write resulting mesh
        if (!overwrite)
        {
            runTime++;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        Info<< nl << "Writing collapsed mesh to time "
            << runTime.timeName() << nl << endl;

        mesh.write();
    }


    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
