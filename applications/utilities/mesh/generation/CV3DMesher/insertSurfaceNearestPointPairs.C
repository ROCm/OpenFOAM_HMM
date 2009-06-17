/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "CV3D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::CV3D::dualCellSurfaceIntersection
(
    const Triangulation::Finite_vertices_iterator& vit
) const
{
    std::list<Facet>  facets;
    incident_facets(vit, std::back_inserter(facets));

    for
    (
        std::list<Facet>::iterator fit=facets.begin();
        fit != facets.end();
        ++fit
    )
    {
        if
        (
            is_infinite(fit->first)
         || is_infinite(fit->first->neighbor(fit->second))
        )
        {
            return true;
        }

        point dE0 = topoint(dual(fit->first));

        // If edge end is outside bounding box then edge cuts boundary
        if (!qSurf_.bb().contains(dE0))
        {
            return true;
        }

        point dE1 = topoint(dual(fit->first->neighbor(fit->second)));

        // If other edge end is outside bounding box then edge cuts boundary
        if (!qSurf_.bb().contains(dE1))
        {
            return true;
        }

        if (magSqr(dE1 - dE0) > tols_.minEdgeLen2)
        {
            pointIndexHit pHit = qSurf_.tree().findLineAny(dE0, dE1);

            if (pHit.hit())
            {
                return true;
            }
        }
    }

    return false;
}


void Foam::CV3D::smoothEdgePositions
(
    DynamicList<point>& edgePoints,
    DynamicList<label>& edgeLabels
) const
{
    const pointField& localPts = qSurf_.localPoints();
    const edgeList& edges = qSurf_.edges();

    // Sort by edge label then by distance from the start of the edge

    SortableList<label> sortedEdgeLabels(edgeLabels);

    labelList sortingIndices = sortedEdgeLabels.indices();

    {
        labelList copyEdgeLabels(edgeLabels.size());

        forAll(sortingIndices, sI)
        {
            copyEdgeLabels[sI] = edgeLabels[sortingIndices[sI]];
        }

        edgeLabels.transfer(copyEdgeLabels);
    }

    {
        List<point> copyEdgePoints(edgePoints.size());

        forAll(sortingIndices, sI)
        {
            copyEdgePoints[sI] = edgePoints[sortingIndices[sI]];
        }

        edgePoints.transfer(copyEdgePoints);
    }

    List<scalar> edgeDistances(edgePoints.size());

    forAll(edgeDistances, eD)
    {
        const point& edgeStart = localPts[edges[edgeLabels[eD]].start()];

        edgeDistances[eD] = mag(edgeStart - edgePoints[eD]);
    }

    // Sort by edgeDistances in blocks of edgeLabel

    DynamicList<label> edgeLabelJumps;

    // Force first edgeLabel to be a jump
    edgeLabelJumps.append(0);

    for (label eL = 1; eL < edgeLabels.size(); eL++)
    {
        if (edgeLabels[eL] > edgeLabels[eL-1])
        {
            edgeLabelJumps.append(eL);
        }
    }

    edgeLabelJumps.shrink();

    forAll(edgeLabelJumps, eLJ)
    {
        label start = edgeLabelJumps[eLJ];

        label length;

        if (eLJ == edgeLabelJumps.size() - 1)
        {
            length = edgeLabels.size() - start;
        }
        else
        {
            length = edgeLabelJumps[eLJ + 1] - start;
        }

        SubList<scalar> edgeDistanceBlock(edgeDistances, length, start);

        SortableList<scalar> sortedEdgeDistanceBlock(edgeDistanceBlock);

        forAll(sortedEdgeDistanceBlock, sEDB)
        {
            sortingIndices[start + sEDB] =
                sortedEdgeDistanceBlock.indices()[sEDB] + start;
        }
    }

    {
        List<point> copyEdgePoints(edgePoints.size());

        forAll(sortingIndices, sI)
        {
            copyEdgePoints[sI] = edgePoints[sortingIndices[sI]];
        }

        edgePoints.transfer(copyEdgePoints);
    }

    {
        List<scalar> copyEdgeDistances(edgeDistances.size());

        forAll(sortingIndices, sI)
        {
            copyEdgeDistances[sI] = edgeDistances[sortingIndices[sI]];
        }

        edgeDistances.transfer(copyEdgeDistances);
    }

    // Create a List to hold each edge, process them individually, then
    // recombine the edges into a single list.

    List<List<point> > edgePointIndividualLists(edgeLabelJumps.size());

    List<List<scalar> > edgeDistanceIndividualLists(edgeLabelJumps.size());

    List<label> edgeLabelIndividualList(edgeLabelJumps.size());

    forAll(edgeLabelJumps, eLJ)
    {
        label start = edgeLabelJumps[eLJ];

        label edgeI = edgeLabels[start];

        label length;

        if (eLJ == edgeLabelJumps.size() - 1)
        {
            length = edgeLabels.size() - start;
        }
        else
        {
            length = edgeLabelJumps[eLJ + 1] - start;
        }

        edgePointIndividualLists[eLJ] = SubList<point>
        (
            edgePoints,
            length,
            start
        );

        edgeDistanceIndividualLists[eLJ] = SubList<scalar>
        (
            edgeDistances,
            length,
            start
        );

        edgeLabelIndividualList[eLJ] = edgeI;
    }

    edgePoints.clear();
    edgeDistances.clear();
    edgeLabels.clear();

    forAll(edgeLabelIndividualList, e)
    {
        label edgeI = edgeLabelIndividualList[e];

        smoothEdge
        (
            edgePointIndividualLists[e],
            edgeDistanceIndividualLists[e],
            edgeI
        );

        const List<point>& tempEdgePoints = edgePointIndividualLists[e];

        forAll(tempEdgePoints, tEP)
        {
            edgePoints.append(tempEdgePoints[tEP]);

            edgeLabels.append(edgeI);
        }
    }

    edgePoints.shrink();

    edgeLabels.shrink();
}



void Foam::CV3D::smoothEdge
(
    List<point>& edgePoints,
    List<scalar>& edgeDistances,
    const label edgeI
) const
{
    const pointField& localPts = qSurf_.localPoints();
    const edgeList& edges = qSurf_.edges();

    // Process the points along each edge (in blocks of edgeLabel) performing 3
    // functions:
    // 1: move points away from feature points
    // 2: aggregate tight groups of points into one point
    // 3: adjust the spacing of remaining points on a pair by pair basis to
    //    remove excess points and add points to long uncontrolled spans.


    const edge& e(edges[edgeI]);

    const point& eStart(localPts[e.start()]);

    const point& eEnd(localPts[e.end()]);

    scalar edgeLength = mag(eStart - eEnd);

    if (edgeLength < 2*tols_.featurePointGuard)
    {
        Info<< "edge " << edgeI
            << " is too short with respect to the featurePointGuard "
            << "distance to allow edge control points to be placed."
            << nl << "Edge length = " << edgeLength
            << nl <<endl;

        return;
    }

    // part 1
    {
        DynamicList<point> tempEdgePoints;

        bool startGuardPlaced = false;

        bool endGuardPlaced = false;

        forAll (edgePoints, eP)
        {
            const point& edgePoint = edgePoints[eP];

            const scalar& edgeDist = edgeDistances[eP];

            if
            (
                edgeDist < tols_.featurePointGuard
             && !startGuardPlaced
            )
            {
                tempEdgePoints.append
                (
                    eStart + (edgePoint - eStart)
                    * tols_.featurePointGuard/edgeDist
                );

                startGuardPlaced = true;
            }
            else if
            (
                edgeDist > (edgeLength - tols_.featurePointGuard)
             && !endGuardPlaced
            )
            {
                tempEdgePoints.append
                (
                    eEnd + (edgePoint - eEnd)
                    * tols_.featurePointGuard/(edgeLength - edgeDist)
                );

                endGuardPlaced = true;
            }
            else if
            (
                edgeDist > tols_.featurePointGuard
             && edgeDist < (edgeLength - tols_.featurePointGuard)
            )
            {
                tempEdgePoints.append(edgePoint);
            }
        }

        edgePoints.transfer(tempEdgePoints.shrink());
    }

    // Recalculate edge distances.

    edgeDistances.setSize(edgePoints.size());

    forAll(edgeDistances, eD)
    {
        edgeDistances[eD] = mag(eStart - edgePoints[eD]);
    }

    // part 2
    {
        DynamicList<point> tempEdgePoints;

        label groupSize = 0;

        point newEdgePoint(vector::zero);

        // if (edgePoints.size() == 1)
        // {
        //     tempEdgePoints.append(edgePoints[0]);
        // }
        // else if
        // (
        //     (edgeDistances[1] - edgeDistances[0]) > tols_.edgeGroupSpacing
        // )
        // {
        //     tempEdgePoints.append(edgePoints[0]);
        // }

        if (edgePoints.size() > 1)
        {
            if ((edgeDistances[1] - edgeDistances[0]) < tols_.edgeGroupSpacing)
            {
                // ...the first two points on the edge start a group

                newEdgePoint += edgePoints[0];

                groupSize++;
            }
            else
            {
                tempEdgePoints.append(edgePoints[0]);
            }
        }
        else
        {
            // ...add the first point by default

            tempEdgePoints.append(edgePoints[0]);
        }

        for (label eP = 1; eP < edgePoints.size(); eP++)
        {
            const scalar& edgeDist = edgeDistances[eP];
            const scalar& previousEdgeDist = edgeDistances[eP - 1];

            if ((edgeDist - previousEdgeDist) < tols_.edgeGroupSpacing)
            {
                newEdgePoint += edgePoints[eP];

                groupSize++;
            }
            else if (groupSize > 0)
            {
                // A point group has been formed and has finished

                newEdgePoint /= groupSize;

                tempEdgePoints.append(newEdgePoint);

                newEdgePoint = vector::zero;

                groupSize = 0;
            }
            else
            {
                tempEdgePoints.append(edgePoints[eP]);
            }
        }

        if (groupSize > 0)
        {
            // A point group has been formed at the end of the edge and needs to
            // be finished.

            newEdgePoint /= groupSize;

            tempEdgePoints.append(newEdgePoint);
        }

        edgePoints.transfer(tempEdgePoints.shrink());
    }

    // Recalculate edge distances.

    edgeDistances.setSize(edgePoints.size());

    forAll(edgeDistances, eD)
    {
        edgeDistances[eD] = mag(eStart - edgePoints[eD]);
    }

    // part 3
    {
        // Special treatment for gaps between closest point to start

        DynamicList<point> tempEdgePoints;

        if (edgeDistances[0] - tols_.featurePointGuard > tols_.maxEdgeSpacing)
        {
            scalar gap = edgeDistances[0] - tols_.featurePointGuard;

            label nInsertions = label(gap/tols_.maxEdgeSpacing);

            // Info<< "Gap at start of edge of "  << gap
            //     << ". Inserting " << nInsertions << " points" << endl;

            scalar spacing = gap / (nInsertions + 1);

            for (label nI = 1; nI <= nInsertions; nI++)
            {
                tempEdgePoints.append
                (
                    eStart + (eEnd - eStart)
                    * (nI * spacing + tols_.featurePointGuard) /edgeLength
                );
            }
        }

        // Identify gaps in middle of edges.
        // Insert first point by default.

        tempEdgePoints.append(edgePoints[0]);

        for (label eP = 1; eP < edgePoints.size(); eP++)
        {
            const scalar& edgeDist = edgeDistances[eP];
            const scalar& previousEdgeDist = edgeDistances[eP - 1];

            if ((edgeDist - previousEdgeDist) > tols_.maxEdgeSpacing)
            {
                scalar gap = edgeDist - previousEdgeDist;

                label nInsertions = label(gap/tols_.maxEdgeSpacing);

                // Info<< "Gap in edge of "  << gap
                //     << ". Inserting " << nInsertions << " points" << endl;

                scalar spacing = gap / (nInsertions + 1);

                for (label nI = 1; nI<= nInsertions; nI++)
                {
                    tempEdgePoints.append
                    (
                        eStart + (eEnd - eStart)
                        * (nI * spacing + previousEdgeDist) /edgeLength
                    );
                }


            }

            tempEdgePoints.append(edgePoints[eP]);
        }

        // Special treatment for gaps between closest point to end

        if
        (
            (edgeLength - edgeDistances[edgeDistances.size() - 1]
                - tols_.featurePointGuard)
                > tols_.maxEdgeSpacing
        )
        {
            scalar lastPointDist =  edgeDistances[edgeDistances.size() - 1];

            const point& lastPoint = edgePoints[edgePoints.size() - 1];

            scalar gap = edgeLength - lastPointDist - tols_.featurePointGuard;

            label nInsertions = label(gap/tols_.maxEdgeSpacing);

            // Info<< "Gap at end of edge of "  << gap
            //     << ". Inserting " << nInsertions << " points" << endl;

            scalar spacing = gap / (nInsertions + 1);

            for (label nI = 1; nI <= nInsertions; nI++)
            {
                tempEdgePoints.append
                (
                    lastPoint + (eEnd - lastPoint)
                        * nI * spacing / gap
                );
            }
        }

        edgePoints.transfer(tempEdgePoints.shrink());

        // Remove pairs of points that are too close together.

        label nPointsRemoved = 1;

        while (nPointsRemoved > 0)
        {
            nPointsRemoved = 0;

            // Recalculate edge distances.

            edgeDistances.setSize(edgePoints.size());

            forAll(edgeDistances, eD)
            {
                edgeDistances[eD] = mag(eStart - edgePoints[eD]);
            }

            // Insert first point
            tempEdgePoints.append(edgePoints[0]);

            bool previousPointMustBeKept = false;

            for (label eP = 1; eP < edgePoints.size(); eP++)
            {
                const scalar& edgeDist = edgeDistances[eP];
                const scalar& previousEdgeDist = edgeDistances[eP - 1];

                if ((edgeDist - previousEdgeDist) < tols_.minEdgeSpacing)
                {
                    if (!previousPointMustBeKept)
                    {
                        tempEdgePoints.remove();
                    }

                    const point& currentPoint = edgePoints[eP];

                    const point& previousPoint = edgePoints[eP - 1];

                    tempEdgePoints.append(0.5*(previousPoint + currentPoint));

                    nPointsRemoved++;

                    previousPointMustBeKept = true;
                }
                else
                {
                    tempEdgePoints.append(edgePoints[eP]);

                    previousPointMustBeKept = false;
                }
            }

            // Info<< edgeI << tab
            //<< nPointsRemoved << " points removed." << endl;

            edgePoints.transfer(tempEdgePoints.shrink());
        }
    }
}


void Foam::CV3D::insertPointPairs
(
    const DynamicList<point>& nearSurfacePoints,
    const DynamicList<point>& surfacePoints,
    const DynamicList<label>& surfaceTris,
    const fileName fName
)
{
    if (controls_.mirrorPoints)
    {
        forAll(surfacePoints, ppi)
        {
            insertMirrorPoint
            (
                nearSurfacePoints[ppi],
                surfacePoints[ppi]
            );
        }
    }
    else
    {
        forAll(surfacePoints, ppi)
        {
            insertPointPair
            (
                tols_.ppDist,
                surfacePoints[ppi],
                qSurf_.faceNormals()[surfaceTris[ppi]]
            );
        }
    }

    Info<< surfacePoints.size() << " point-pairs inserted" << endl;

    if (controls_.writeInsertedPointPairs)
    {
        OFstream str(fName);
        label vertI = 0;

        forAll(surfacePoints, ppi)
        {
            meshTools::writeOBJ(str, surfacePoints[ppi]);
            vertI++;
        }

        Info<< "insertPointPairs: Written " << surfacePoints.size()
            << " inserted point-pair locations to file "
            << str.name() << endl;
    }
}


void Foam::CV3D::insertEdgePointGroups
(
    const DynamicList<point>& edgePoints,
    const DynamicList<label>& edgeLabels,
    const fileName fName
)
{
    const pointField& localPts = qSurf_.localPoints();

    forAll(edgePoints, eP)
    {
        const point& edgePt = edgePoints[eP];

        const label edgeI = edgeLabels[eP];

        // Pick up the two faces adjacent to the feature edge
        const labelList& eFaces = qSurf_.edgeFaces()[edgeI];

        label faceA = eFaces[0];
        vector nA = qSurf_.faceNormals()[faceA];

        label faceB = eFaces[1];
        vector nB = qSurf_.faceNormals()[faceB];

        // Intersect planes parallel to faceA and faceB offset by ppDist
        // and the plane defined by edgePt and the edge vector.
        plane planeA(edgePt - tols_.ppDist*nA, nA);
        plane planeB(edgePt - tols_.ppDist*nB, nB);

        // Finding the nearest point on the intersecting line to the edge point.
        // Floating point errors often encountered using planePlaneIntersect
        // plane planeF(edgePt, (nA^nB));
        // point refPt = planeF.planePlaneIntersect(planeA,planeB);

        plane::ray planeIntersect(planeA.planeIntersect(planeB));

        pointHit refPtHit = linePointRef
        (
            planeIntersect.refPoint() + 2*tols_.span*planeIntersect.dir(),
            planeIntersect.refPoint() - 2*tols_.span*planeIntersect.dir()
        ).nearestDist(edgePt);

        point refPt = refPtHit.hitPoint();

        point faceAVert =
        localPts[triSurfaceTools::oppositeVertex(qSurf_, faceA, edgeI)];

        // Determine convex or concave angle
        if (((faceAVert - edgePt) & nB) < 0)
        {
            // Convex. So refPt will be inside domain and hence a master point

            // Insert the master point refering the the first slave
            label masterPtIndex = insertPoint(refPt, number_of_vertices() + 1);

            // Insert the slave points by reflecting refPt in both faces.
            // with each slave refering to the master

            point reflectedA = refPt + 2*((edgePt - refPt) & nA)*nA;
            insertPoint(reflectedA, masterPtIndex);

            point reflectedB = refPt + 2*((edgePt - refPt) & nB)*nB;
            insertPoint(reflectedB, masterPtIndex);
        }
        else
        {
            // Concave. master and reflected points inside the domain.
            // Generate reflected master to be outside.
            point reflMasterPt = refPt + 2*(edgePt - refPt);

            // Reflect refPt in both faces.
            point reflectedA =
            reflMasterPt + 2*((edgePt - reflMasterPt) & nA)*nA;

            point reflectedB =
            reflMasterPt + 2*((edgePt - reflMasterPt) & nB)*nB;

            scalar totalAngle =
            180*(mathematicalConstant::pi + acos(mag(nA & nB)))
            /mathematicalConstant::pi;

            // Number of quadrants the angle should be split into
            int nQuads = int(totalAngle/controls_.maxQuadAngle) + 1;

            // The number of additional master points needed to obtain the
            // required number of quadrants.
            int nAddPoints = min(max(nQuads - 2, 0), 2);

            // index of reflMaster
            label reflectedMaster = number_of_vertices() + 2 + nAddPoints;

            // Master A is inside.
            label reflectedAI = insertPoint(reflectedA, reflectedMaster);

            // Master B is inside.
            insertPoint(reflectedB, reflectedMaster);

            if (nAddPoints == 1)
            {
                // One additinal point is the reflection of the slave point,
                // i.e. the original reference point
                insertPoint(refPt, reflectedMaster);
            }
            else if (nAddPoints == 2)
            {
                point reflectedAa = refPt
                - ((edgePt - reflMasterPt) & nB)*nB;
                insertPoint(reflectedAa, reflectedMaster);

                point reflectedBb = refPt
                - ((edgePt - reflMasterPt) & nA)*nA;
                insertPoint(reflectedBb, reflectedMaster);
            }

            // Slave is outside.
            insertPoint(reflMasterPt, reflectedAI);
        }
    }

    Info<< edgePoints.size() << " edge-control locations inserted" << endl;

    if (controls_.writeInsertedPointPairs)
    {
        OFstream str(fName);
        label vertI = 0;

        forAll(edgePoints, eP)
        {
            meshTools::writeOBJ(str, edgePoints[eP]);
            vertI++;
        }

        Info<< "insertEdgePointGroups: Written " << edgePoints.size()
            << " inserted edge-control locations to file "
            << str.name() << endl;
    }
}


void Foam::CV3D::insertSurfaceNearestPointPairs()
{
    Info<< nl << "insertSurfaceNearestPointPairs: " << endl;

    label nSurfacePointsEst = number_of_vertices();

    scalar distanceFactor = 8.0;
    scalar distanceFactor2 = Foam::sqr(distanceFactor);

    if (qSurf_.features().featureEdges().size())
    {
        DynamicList<point> nearSurfacePoints(nSurfacePointsEst);

        for
        (
            Triangulation::Finite_vertices_iterator vit =
            finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->internalPoint())
            {
                point vert(topoint(vit->point()));

                pointIndexHit pHit = qSurf_.tree().findNearest
                (
                    vert,
                    distanceFactor2*controls_.minCellSize2
                );

                if (pHit.hit())
                {
                    vit->setNearBoundary();

                    if (dualCellSurfaceIntersection(vit))
                    {
                        nearSurfacePoints.append(vert);
                    }
                }
            }
        }

        pointField nearSurfacePointsForEdges(nearSurfacePoints.shrink());

        labelList allEdgeLabels;
        labelList allEdgeEndPoints;
        pointField allEdgePoints;

        qSurf_.features().nearestSurfEdge
        (
            qSurf_.features().featureEdges(),
            nearSurfacePointsForEdges,
            vector::one * distanceFactor * controls_.minCellSize,
            allEdgeLabels,
            allEdgeEndPoints,
            allEdgePoints
        );

        DynamicList<label> edgeLabels(allEdgeLabels.size());
        DynamicList<vector> edgePoints(allEdgePoints.size());

        forAll(allEdgePoints, eP)
        {
            if (allEdgeLabels[eP] >= 0 && allEdgeEndPoints[eP] < 0)
            {
                edgeLabels.append(allEdgeLabels[eP]);
                edgePoints.append(allEdgePoints[eP]);
            }
        }

        edgePoints.shrink();
        edgeLabels.shrink();

        if (edgePoints.size())
        {
            //Warning<< "Edge point insertion disabled." << endl;

            smoothEdgePositions(edgePoints, edgeLabels);

            insertEdgePointGroups
            (
                edgePoints,
                edgeLabels,
                "surfaceNearestEdgePoints.obj"
            );
        }
    }

    DynamicList<point> allNearSurfacePoints(nSurfacePointsEst);
    DynamicList<point> allSurfacePoints(nSurfacePointsEst);
    DynamicList<label> allSurfaceTris(nSurfacePointsEst);

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalPoint())
        {
            point vert(topoint(vit->point()));

            pointIndexHit pHit = qSurf_.tree().findNearest
            (
                vert,
                distanceFactor2*controls_.minCellSize2
            );

            if (pHit.hit())
            {
                vit->setNearBoundary();

                if (dualCellSurfaceIntersection(vit))
                {
                    allNearSurfacePoints.append(vert);
                    allSurfacePoints.append(pHit.hitPoint());
                    allSurfaceTris.append(pHit.index());
                }
            }
        }
    }

    pointField surfacePointsForEdges(allSurfacePoints.shrink());

    labelList allEdgeLabels;
    labelList allEdgeEndPoints;
    pointField allEdgePoints;

    qSurf_.features().nearestSurfEdge
    (
        qSurf_.features().featureEdges(),
        surfacePointsForEdges,
        vector::one * 2.0 * tols_.featureEdgeGuard,
        allEdgeLabels,
        allEdgeEndPoints,
        allEdgePoints
    );

    DynamicList<point> nearSurfacePoints(nSurfacePointsEst);
    DynamicList<point> surfacePoints(nSurfacePointsEst);
    DynamicList<label> surfaceTris(nSurfacePointsEst);

    forAll(allEdgePoints, eP)
    {
        if (allEdgeLabels[eP] == -1)
        {
            nearSurfacePoints.append(allNearSurfacePoints[eP]);
            surfacePoints.append(allSurfacePoints[eP]);
            surfaceTris.append(allSurfaceTris[eP]);
        }
    }

    Info<< nl << "Number of surface conformation points not placed because "
        << nl << "they were too close to a feature edge = "
        << allSurfacePoints.size() - surfacePoints.size()
        << endl;

    nearSurfacePoints.shrink();
    surfacePoints.shrink();
    surfaceTris.shrink();

    insertPointPairs
    (
        nearSurfacePoints,
        surfacePoints,
        surfaceTris,
        "surfaceNearestIntersections.obj"
    );
}


// ************************************************************************* //
