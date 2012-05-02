/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

    Cannot remove cells. Can remove faces and points but does not check
    for nonsense resulting topology.

    When collapsing an edge with one point on the boundary it will leave
    the boundary point intact. When both points inside it chooses random. When
    both points on boundary random again. Note: it should in fact use features
    where if one point is on a feature it collapses to that one. Alas we don't
    have features on a polyMesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "mapPolyMesh.H"
#include "mathematicalConstants.H"
#include "PackedBoolList.H"
#include "SortableList.H"
#include "unitConversion.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "PointEdgeWave.H"
#include "pointEdgeCollapse.H"
#include "motionSmoother.H"

#include "OFstream.H"
#include "meshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label findIndex
(
    const labelList& elems,
    const label start,
    const label size,
    const label val
)
{
    for (label i = start; i < size; i++)
    {
        if (elems[i] == val)
        {
            return i;
        }
    }
    return -1;
}


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
            if (findIndex(f, 0, newFp, collapsePoint) == -1)
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
        label index = findIndex(f, 0, fp, pointI);

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

    // Create strings of edges
    Map<DynamicList<label> > collapseStrings;

    forAll(allPointInfo, pointI)
    {
        const label collapseIndex = allPointInfo[pointI].collapseIndex();

        if (collapseIndex != -1 && collapseIndex != -2)
        {
            collapseStrings(collapseIndex).append(pointI);
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
                            Info<< "Cell:" << cellI
                                << " uses faces:" << cFaces
                                << " of which too many are marked for removal:"
                                << endl
                                << "   ";
                            forAll(cFaces, j)
                            {
                                if (newFaces[cFaces[j]].size() < 3)
                                {
                                    Info<< ' '<< cFaces[j];
                                }
                            }
                            Info<< endl;

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

    // Print regions:
//    printRegions();

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


// Merges edges which are in straight line. I.e. edge split by point.
//label mergeEdges
//(
//    const polyMesh& mesh,
//    const scalar maxCos,
//    List<pointEdgeCollapse>& allPointInfo
//)
//{
//    const pointField& points = mesh.points();
//    const edgeList& edges = mesh.edges();
//    const labelListList& pointEdges = mesh.pointEdges();
//    const labelList& region = collapser.pointRegion();
//    const labelList& master = collapser.pointRegionMaster();
//
//    label nCollapsed = 0;
//
//    forAll(pointEdges, pointI)
//    {
//        const labelList& pEdges = pointEdges[pointI];
//
//        if (pEdges.size() == 2)
//        {
//            const edge& leftE = edges[pEdges[0]];
//            const edge& rightE = edges[pEdges[1]];
//
//            // Get the two vertices on both sides of the point
//            label leftV = leftE.otherVertex(pointI);
//            label rightV = rightE.otherVertex(pointI);
//
//            // Collapse only if none of the points part of merge network
//            // or all of networks with different masters.
//            label midMaster = -1;
//            if (region[pointI] != -1)
//            {
//                midMaster = master[region[pointI]];
//            }
//
//            label leftMaster = -2;
//            if (region[leftV] != -1)
//            {
//                leftMaster = master[region[leftV]];
//            }
//
//            label rightMaster = -3;
//            if (region[rightV] != -1)
//            {
//                rightMaster = master[region[rightV]];
//            }
//
//            if
//            (
//                midMaster != leftMaster
//             && midMaster != rightMaster
//             && leftMaster != rightMaster
//            )
//            {
//                // Check if the two edge are in line
//                vector leftVec = points[pointI] - points[leftV];
//                leftVec /= mag(leftVec) + VSMALL;
//
//                vector rightVec = points[rightV] - points[pointI];
//                rightVec /= mag(rightVec) + VSMALL;
//
//                if ((leftVec & rightVec) > maxCos)
//                {
//                    // Collapse one (left) side of the edge. Make left vertex
//                    // the master.
//                    //if (collapser.unaffectedEdge(pEdges[0]))
//                    const edge& e = mesh.edges()[pEdges[0]];
//
//                    if
//                    (
//                        allPointInfo[e[0]].collapseIndex() < 0
//                     && allPointInfo[e[1]].collapseIndex() < 0
//                    )
//                    {
//                        //pointEdgeCollapse pec(mesh.points()[leftV], leftV);
//
//                        allPointInfo[e[0]] = pec;
//                        allPointInfo[e[1]] = pec;
//
//                        //collapser.collapseEdge(pEdges[0], leftV);
//                        nCollapsed++;
//                    }
//                }
//            }
//        }
//    }
//
//    return nCollapsed;
//}


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


label collapseSmallEdges
(
    const polyMesh& mesh,
    const scalarList& freezeEdges,
    const labelList& boundaryPoint,
    const scalar minLen,
    List<pointEdgeCollapse>& allPointInfo
)
{
    const pointField& points = mesh.points();
    const edgeList& edges = mesh.edges();

    // Store collapse direction in collapseEdge
    //   -1 -> Do not collapse
    //    0 -> Collapse to start point
    //    1 -> Collapse to end point
    labelList collapseEdge(edges.size(), -1);

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (e.mag(points) < minLen*freezeEdges[edgeI])
        {
            collapseEdge[edgeI] = 0;
        }
    }

    // Check whether edge point order is reversed from mesh to coupledPatch
//    const globalMeshData& globalData = mesh.globalData();
//    const mapDistribute& map = globalData.globalEdgeSlavesMap();
//    const labelList& coupledMeshEdges = globalData.coupledPatchMeshEdges();
//    const indirectPrimitivePatch& coupledPatch = globalData.coupledPatch();
//    const PackedBoolList& cppOrientation = globalData.globalEdgeOrientation();
//    PackedBoolList meshToPatchSameOrientation(coupledMeshEdges.size(), true);
//
//    forAll(coupledMeshEdges, eI)
//    {
//        const label meshEdgeIndex = coupledMeshEdges[eI];
//
//        if (collapseEdge[meshEdgeIndex] != -1)
//        {
//            const edge& meshEdge = edges[meshEdgeIndex];
//            const edge& coupledPatchEdge = coupledPatch.edges()[eI];
//
//            if
//            (
//                meshEdge[0] == coupledPatch.meshPoints()[coupledPatchEdge[1]]
//             && meshEdge[1] == coupledPatch.meshPoints()[coupledPatchEdge[0]]
//            )
//            {
//                meshToPatchSameOrientation[eI] = false;
//            }
//        }
//    }
//
//
//    labelList cppEdgeData(coupledMeshEdges.size(), -1);
//
//    forAll(coupledMeshEdges, eI)
//    {
//        const label meshEdgeIndex = coupledMeshEdges[eI];
//
//        if (collapseEdge[meshEdgeIndex] != -1)
//        {
//            if (meshToPatchSameOrientation[eI] == cppOrientation[eI])
//            {
//                cppEdgeData[eI] = 0;
//            }
//            else
//            {
//                cppEdgeData[eI] = 1;
//            }
//        }
//    }
//
//
//    // Synchronise cppEdgeData
//    // Use minEqOp reduction, so that edge will only be collapsed on processor
//    // boundary if both processors agree to collapse it
//    globalData.syncData
//    (
//        cppEdgeData,
//        globalData.globalEdgeSlaves(),
//        globalData.globalEdgeTransformedSlaves(),
//        map,
//        minEqOp<label>()
//    );
//
//
//    forAll(coupledMeshEdges, eI)
//    {
//        const label meshEdgeIndex = coupledMeshEdges[eI];
//
//        if (collapseEdge[meshEdgeIndex] != -1)
//        {
//            if (meshToPatchSameOrientation[eI] == cppOrientation[eI])
//            {
//                collapseEdge[meshEdgeIndex] = 0;
//            }
//            else
//            {
//                collapseEdge[meshEdgeIndex] = 1;
//            }
//        }
//    }

    label nCollapsed = 0;

    DynamicList<label> initPoints(mesh.nPoints());
    DynamicList<pointEdgeCollapse> initPointInfo(mesh.nPoints());
    allPointInfo.resize(mesh.nPoints());

    globalIndex globalStrings(mesh.nPoints());

    List<pointEdgeCollapse> allEdgeInfo(mesh.nEdges());
    forAll(allEdgeInfo, edgeI)
    {
        allEdgeInfo[edgeI] = pointEdgeCollapse(vector::zero, -1);
    }

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (collapseEdge[edgeI] != -1)
        {
            const label master =
                edgeMaster
                (
                    boundaryPoint,
                    collapseEdge[edgeI],
                    e
                );

//            if (master != -1)
            {
                pointEdgeCollapse pec
                (
                    points[master],
                    globalStrings.toGlobal(master)
                );

                allEdgeInfo[edgeI] = pec;

                initPointInfo.append(pec);
                initPoints.append(e.start());

                initPointInfo.append(pec);
                initPoints.append(e.end());

                nCollapsed++;
            }
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


// Faces which have edges just larger than collapse length but faces which
// are very small. This one tries to collapse them if it can be done with
// edge collapse. For faces where a face gets replace by two edges use
// collapseFaces
//label collapseHighAspectFaces
//(
//    const polyMesh& mesh,
//    const PackedBoolList& boundaryPoint,
//    const Map<label>& processorPoints,
//    const scalar areaFac,
//    const scalar edgeRatio,
//    edgeCollapser& collapser
//)
//{
//    const pointField& points = mesh.points();
//    const edgeList& edges = mesh.edges();
//    const faceList& faces = mesh.faces();
//    const labelListList& faceEdges = mesh.faceEdges();
//
//    scalarField magArea(mag(mesh.faceAreas()));
//
//    label maxIndex = findMax(magArea);
//
//    scalar minArea = areaFac * magArea[maxIndex];
//
//    Info<< "Max face area:" << magArea[maxIndex] << endl
//        << "Collapse area factor:" << areaFac << endl
//        << "Collapse area:" << minArea << endl;
//
//    label nCollapsed = 0;
//
//    forAll(faces, faceI)
//    {
//        if (magArea[faceI] < minArea)
//        {
//            const face& f = faces[faceI];
//
//            // Get the edges in face point order
//            labelList fEdges(getSortedEdges(edges, f, faceEdges[faceI]));
//
//            SortableList<scalar> lengths(fEdges.size());
//            forAll(fEdges, i)
//            {
//                lengths[i] = edges[fEdges[i]].mag(points);
//            }
//            lengths.sort();
//
//
//            label edgeI = -1;
//
//            if (f.size() == 4)
//            {
//                // Compare second largest to smallest
//                if (lengths[2] > edgeRatio*lengths[0])
//                {
//                    // Collapse smallest only. Triangle should be cleared
//                    // next time around.
//                    edgeI = fEdges[lengths.indices()[0]];
//                }
//            }
//            else if (f.size() == 3)
//            {
//                // Compare second largest to smallest
//                if (lengths[1] > edgeRatio*lengths[0])
//                {
//                    edgeI = fEdges[lengths.indices()[0]];
//                }
//            }
//
//
//            if (edgeI != -1)
//            {
//                label master = edgeMaster
//                    (
//                        boundaryPoint,
//                        processorPoints,
//                        false,
//                        edges[edgeI]
//                    );
//
//                if (master != -1)// && collapser.unaffectedEdge(edgeI))
//                {
//                    collapser.collapseEdge(edgeI, master);
//                    nCollapsed++;
//                }
//            }
//        }
//    }
//
//    return nCollapsed;
//}


void set(const labelList& elems, const bool val, boolList& status)
{
    forAll(elems, i)
    {
        status[elems[i]] = val;
    }
}


// Tries to simplify polygons to face of minSize (4=quad, 3=triangle)
//label simplifyFaces
//(
//    const polyMesh& mesh,
//    const PackedBoolList& boundaryPoint,
//    const Map<label>& processorPoints,
//    const label minSize,
//    const scalar lenGap,
//    edgeCollapser& collapser
//)
//{
//    const pointField& points = mesh.points();
//    const edgeList& edges = mesh.edges();
//    const faceList& faces = mesh.faces();
//    const cellList& cells = mesh.cells();
//    const labelListList& faceEdges = mesh.faceEdges();
//    const labelList& faceOwner = mesh.faceOwner();
//    const labelList& faceNeighbour = mesh.faceNeighbour();
//    const labelListList& pointCells = mesh.pointCells();
//    const labelListList& cellEdges = mesh.cellEdges();
//
//    label nCollapsed = 0;
//
//    boolList protectedEdge(mesh.nEdges(), false);
//
//    forAll(faces, faceI)
//    {
//        const face& f = faces[faceI];
//
//        if
//        (
//            f.size() > minSize
//         && cells[faceOwner[faceI]].size() >= 6
//         && (
//                mesh.isInternalFace(faceI)
//             && cells[faceNeighbour[faceI]].size() >= 6
//            )
//        )
//        {
//            // Get the edges in face point order
//            labelList fEdges(getSortedEdges(edges, f, faceEdges[faceI]));
//
//            SortableList<scalar> lengths(fEdges.size());
//            forAll(fEdges, i)
//            {
//                lengths[i] = edges[fEdges[i]].mag(points);
//            }
//            lengths.sort();
//
//
//            // Now find a gap in length between consecutive elements greater
//            // than lenGap.
//
//            label gapPos = -1;
//
//            for (label i = f.size()-1-minSize; i >= 0; --i)
//            {
//                if (lengths[i+1] > lenGap*lengths[i])
//                {
//                    gapPos = i;
//
//                    break;
//                }
//            }
//
//            if (gapPos != -1)
//            {
//                //for (label i = gapPos; i >= 0; --i)
//                label i = 0;  // Hack: collapse smallest edge only.
//                {
//                    label edgeI = fEdges[lengths.indices()[i]];
//
//                    if (!protectedEdge[edgeI])
//                    {
//                        const edge& e = edges[edgeI];
//
//                        label master
//                       = edgeMaster(boundaryPoint, processorPoints, false, e);
//
//                        if (master != -1)
//                        {
//                            collapser.collapseEdge(edgeI, master);
//
//                            // Protect all other edges on all cells using edge
//                            // points.
//
//                            const labelList& pCells0 = pointCells[e[0]];
//
//                            forAll(pCells0, i)
//                            {
//                              set(cellEdges[pCells0[i]], true, protectedEdge);
//                            }
//                            const labelList& pCells1 = pointCells[e[1]];
//
//                            forAll(pCells1, i)
//                            {
//                              set(cellEdges[pCells1[i]], true, protectedEdge);
//                            }
//
//                            nCollapsed++;
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    return nCollapsed;
//}


labelHashSet checkMeshQuality(const polyMesh& mesh)
{
    //mesh.checkMesh(true);
    labelHashSet freezePoints;

    IOdictionary meshQualityDict
    (
        IOobject
        (
            "meshQualityControls",
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

    forAllConstIter(labelHashSet, badFaces, iter)
    {
        const face& f = mesh.faces()[iter.key()];

        forAll(f, pI)
        {
            freezePoints.insert(f[pI]);
        }
    }

//    const edgeList& edges = mesh.edges();
//
//    label nFrozenEdges = 0;
//
//    OFstream str("frozenEdges.obj");
//
//    label count = 0;
//    forAll(edges, eI)
//    {
//        const edge& e = edges[eI];
//
//        if (freezePoints.found(e[0]) && freezePoints.found(e[1]))
//        {
//            freezeEdges[eI] = true;
//            nFrozenEdges++;
//        }
//    }

    return freezePoints;
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

        if (isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch& pPatch =
                refCast<const processorPolyPatch>(patch);

            forAll(pPatch, fI)
            {
                const face& f = pPatch[fI];

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
#   include "addOverwriteOption.H"

    argList::validArgs.append("edge length [m]");
    argList::validArgs.append("merge angle (degrees)");
    argList::addOption("minLenFactor", "scalar", "edge length factor");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    scalar minLen  = args.argRead<scalar>(1);
    const scalar angle   = args.argRead<scalar>(2);
    const scalar minLenFactor
        = args.optionLookupOrDefault<scalar>("minLenFactor", 0.5);

    const bool overwrite = args.optionFound("overwrite");

    scalar maxCos = Foam::cos(degToRad(angle));

    Info<< "Merging:" << nl
        << "    edges with length less than " << minLen << " meters" << nl
        << "    edges split by a point with edges in line to within " << angle
        << " degrees" << nl
        << endl;

    Info<< "If an invalid mesh is generated then the edge length will be " << nl
        << "multiplied by a factor of " << minLenFactor << " and collapsing "
        << "will be reattempted" << nl << endl;

    bool meshChanged = false;

    checkMeshQuality(mesh);

    autoPtr<fvMesh> fvMeshPtr;

    scalarList freezeEdges(mesh.nEdges(), 1.0);

    do
    {
        label nIterations = 0;
        label nFrozenEdges = 0;

        fvMeshPtr.reset(new fvMesh(mesh));
        fvMesh& fvMeshRef = fvMeshPtr();

        // Contains new point label for original points
        labelList pointMap(identity(mesh.nPoints()));

        scalarList tmpFreezeEdges = freezeEdges;

        autoPtr<mapPolyMesh> morphMap;

        while (true)
        {
            Info<< "Iteration " << nIterations << incrIndent << endl;

            labelList boundaryPoint = findBoundaryPoints(fvMeshRef);

            List<pointEdgeCollapse> allPointInfo;

            // Collapse all edges that are too small.
            label nSmallCollapsed =
                collapseSmallEdges
                (
                    fvMeshRef,
                    tmpFreezeEdges,
                    boundaryPoint,
                    minLen,
                    allPointInfo
                );

            reduce(nSmallCollapsed, sumOp<label>());

            Info<< indent << "Collapsing " << nSmallCollapsed
                << " small edges" << endl;





            label nMerged = 0;

            // Remove midpoints on straight edges.
            if (nSmallCollapsed == 0)
            {
                //nMerged = mergeEdges(fvMeshRef, maxCos, allPointInfo);
            }

            reduce(nMerged, sumOp<label>());

            Info<< indent << "Collapsing " << nMerged << " in line edges"
                << endl;






            label nSliversCollapsed = 0;

            // Remove small sliver faces that can be collapsed to single edge
//            if (nSmallCollapsed == 0 && nMerged == 0)
//            {
//                nSliversCollapsed =
//                    collapseHighAspectFaces
//                    (
//                        mesh,
//                        boundaryPoint,
//                        processorPoints,
//                        1E-9,
//                        5,
//                        collapser
//                    );
//            }

            reduce(nSliversCollapsed, sumOp<label>());

            Info<< indent << "Collapsing " << nSliversCollapsed
                << " small high aspect ratio faces" << endl;


            // Simplify faces to quads wherever possible
            //if (nCollapsed == 0)
            //{
            //    nCollapsed =
            //        simplifyFaces
            //        (
            //            mesh,
            //            boundaryPoint,
            //            4,              // minimum size of face
            //            0.2,            // gap in edge lengths on face
            //            collapser
            //        );
            //    Info<< "Collapsing " << nCollapsed << " polygonal faces"
            //        << endl;
            //}


            label totalCollapsed =
                nSmallCollapsed
              + nMerged
              + nSliversCollapsed;

            polyTopoChange meshMod(fvMeshRef);

            // Insert mesh refinement into polyTopoChange.
            setRefinement(fvMeshRef, meshMod, allPointInfo);

            // Do all changes
            Info<< indent << "Applying changes to the mesh" << nl
                << decrIndent << endl;

            morphMap = meshMod.changeMesh(fvMeshRef, false);

//            // Contains new point label for old points
//            const labelList& reversePointMap = morphMap().reversePointMap();
//
//            forAll(pointMap, pI)
//            {
//                const label originalPoint = pI;
//                const label currentPoint = pointMap[pI];
//
//                if (currentPoint < reversePointMap.size())
//                {
//                    const label newPoint = reversePointMap[currentPoint];
//
//                    if (newPoint != -1)
//                    {
//                        pointMap[originalPoint] = newPoint;
//                    }
//                }
//            }

            if (totalCollapsed == 0)
            {
                labelHashSet freezePoints = checkMeshQuality(fvMeshRef);

                label nFreezePoints = freezePoints.size();
                reduce(nFreezePoints, sumOp<label>());

                nFrozenEdges = nFreezePoints;

                Info<< "Number of frozen points      : " << nFreezePoints
                    << endl;

                break;
            }

            if (morphMap().hasMotionPoints())
            {
                fvMeshRef.movePoints(morphMap().preMotionPoints());
            }

            meshChanged = true;

            nIterations++;
        }

        if (nFrozenEdges > 0)
        {
            minLen *= minLenFactor;
        }

        reduce(nFrozenEdges, sumOp<label>());

        Info<< "Number of frozen edges       : " << nFrozenEdges << nl
            << endl;

        if (nFrozenEdges == 0)
        {
            break;
        }

    } while (true);


    if (meshChanged)
    {
        // Write resulting mesh
        if (!overwrite)
        {
            runTime++;
        }
        else
        {
            fvMeshPtr().setInstance(oldInstance);
        }

        Info<< nl << "Writing collapsed mesh to time "
            << runTime.timeName() << nl << endl;

        fvMeshPtr().write();
    }

    Info<< "Final minimum length : " << minLen << " m" << nl << endl;

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
