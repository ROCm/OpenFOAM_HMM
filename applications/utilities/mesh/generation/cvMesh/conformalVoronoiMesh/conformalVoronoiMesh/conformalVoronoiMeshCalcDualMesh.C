/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "conformalVoronoiMesh.H"
#include "motionSmoother.H"
#include "backgroundMeshDecomposition.H"
#include "polyMeshGeometry.H"
#include "indexedCellChecks.H"

#include "CGAL/Exact_predicates_exact_constructions_kernel.h"
#include "CGAL/Gmpq.h"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::checkCells()
{
    List<List<FixedList<Foam::point, 4> > > cellListList(Pstream::nProcs());

    List<FixedList<Foam::point, 4> > cells(number_of_finite_cells());

    globalIndex gIndex(number_of_vertices());

    label count = 0;
    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (tetrahedron(cit).volume() == 0)
        {
            Pout<< "ZERO VOLUME TET" << endl;
            Pout<< cit->info();
            Pout<< cit->dual();
        }

        if (cit->hasFarPoint())
        {
            continue;
        }

        List<labelPair> cellVerticesPair(4);
        List<Foam::point> cellVertices(4);

        for (label vI = 0; vI < 4; ++vI)
        {
            cellVerticesPair[vI] = labelPair
            (
                cit->vertex(vI)->procIndex(),
                cit->vertex(vI)->index()
            );
            cellVertices[vI] = topoint(cit->vertex(vI)->point());
        }

        List<Foam::point> cellVerticesOld(cellVertices);
        labelList oldToNew;
        sortedOrder(cellVerticesPair, oldToNew);
        oldToNew = invert(oldToNew.size(), oldToNew);
        inplaceReorder(oldToNew, cellVerticesPair);
        inplaceReorder(oldToNew, cellVertices);

//        FixedList<label, 4> globalTetCell
//        (
//            cit->globallyOrderedCellVertices(gIndex)
//        );
//
//        FixedList<Point, 4> cellVertices(Point(0,0,0));
//
//        forAll(globalTetCell, gvI)
//        {
//            label gI = globalTetCell[gvI];
//
//            cellVertices[gvI] = cit->vertex(gI)->point();
//        }

//        if (cit->hasFarPoint())
//        {
//            continue;
//        }

        for (label i = 0; i < 4; ++i)
        {
            //cells[count][i] = topoint(cit->vertex(i)->point());
            cells[count][i] = cellVertices[i];
        }

        count++;
    }

    cells.setSize(count);

    cellListList[Pstream::myProcNo()] = cells;

    Pstream::gatherList(cellListList);

    if (Pstream::master())
    {
        Info<< "Checking on master processor the cells of each " << nl
            << "processor point list against the master cell list." << nl
            << "There are " << cellListList.size() << " processors" << nl
            << "The size of each processor's cell list is:" << endl;

        forAll(cellListList, cfI)
        {
            Info<< "    Proc " << cfI << " has " << cellListList[cfI].size()
                << " cells" << endl;
        }

        label nMatches = 0, nMatchFoundDiffOrder = 0;

        forAll(cellListList[0], cmI)
        {
            const FixedList<Foam::point, 4>& masterCell = cellListList[0][cmI];

            bool matchFound = false;
            bool matchFoundDiffOrder = false;

            forAll(cellListList, cpI)
            {
                if (cpI == 0)
                {
                    continue;
                }

                forAll(cellListList[cpI], csI)
                {
                    const FixedList<Foam::point, 4>& slaveCell
                        = cellListList[cpI][csI];

                    if (masterCell == slaveCell)
                    {
                        matchFound = true;
                        break;
                    }
                    else
                    {
                        label samePt = 0;

                        forAll(masterCell, mI)
                        {
                            const Foam::point& mPt = masterCell[mI];

                            forAll(slaveCell, sI)
                            {
                                const Foam::point& sPt = slaveCell[sI];

                                if (mPt == sPt)
                                {
                                    samePt++;
                                }
                            }
                        }

                        if (samePt == 4)
                        {
                            matchFoundDiffOrder = true;

                            Pout<< masterCell << nl << slaveCell << endl;

                            break;
                        }
                    }
                }
            }

            if (matchFound)
            {
                nMatches++;
            }

            if (matchFoundDiffOrder)
            {
                nMatchFoundDiffOrder++;
            }
        }

        Info<< "Found " << nMatches << " matching cells and "
            << nMatchFoundDiffOrder << " matching cells with different "
            << "vertex ordering"<< endl;
    }
}


void Foam::conformalVoronoiMesh::checkDuals()
{
    List<List<Point> > pointFieldList(Pstream::nProcs());

    List<Point> duals(number_of_finite_cells());

    typedef CGAL::Exact_predicates_exact_constructions_kernel       EK2;
    typedef CGAL::Regular_triangulation_euclidean_traits_3<EK2>     EK;
    typedef CGAL::Cartesian_converter<baseK::Kernel, EK2>  To_exact;
    typedef CGAL::Cartesian_converter<EK2, baseK::Kernel>  Back_from_exact;

//    PackedBoolList bPoints(number_of_finite_cells());

//    indexDualVertices(duals, bPoints);

    label count = 0;//duals.size();

    duals.setSize(number_of_finite_cells());

    globalIndex gIndex(number_of_vertices());

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (cit->hasFarPoint())
        {
            continue;
        }

        duals[count++] = cit->circumcenter();

//        List<labelPair> cellVerticesPair(4);
//        List<Point> cellVertices(4);
//
//        for (label vI = 0; vI < 4; ++vI)
//        {
//            cellVerticesPair[vI] = labelPair
//            (
//                cit->vertex(vI)->procIndex(),
//                cit->vertex(vI)->index()
//            );
//            cellVertices[vI] = cit->vertex(vI)->point();
//        }
//
//        labelList oldToNew;
//        sortedOrder(cellVerticesPair, oldToNew);
//        oldToNew = invert(oldToNew.size(), oldToNew);
//        inplaceReorder(oldToNew, cellVerticesPair);
//        inplaceReorder(oldToNew, cellVertices);
//
//        duals[count++] = CGAL::circumcenter
//        (
//            cellVertices[0],
//            cellVertices[1],
//            cellVertices[2],
//            cellVertices[3]
//        );

//        To_exact to_exact;
//        Back_from_exact back_from_exact;
//        EK::Construct_circumcenter_3 exact_circumcenter =
//            EK().construct_circumcenter_3_object();
//
//        duals[count++] = topoint
//        (
//            back_from_exact
//            (
//                exact_circumcenter
//                (
//                    to_exact(cit->vertex(0)->point()),
//                    to_exact(cit->vertex(1)->point()),
//                    to_exact(cit->vertex(2)->point()),
//                    to_exact(cit->vertex(3)->point())
//                )
//            )
//        );
    }

    Pout<< "Duals Calculated " << count << endl;

    duals.setSize(count);

    pointFieldList[Pstream::myProcNo()] = duals;

    Pstream::gatherList(pointFieldList);

    if (Pstream::master())
    {
        Info<< "Checking on master processor the dual locations of each " << nl
            << "processor point list against the master dual list." << nl
            << "There are " << pointFieldList.size() << " processors" << nl
            << "The size of each processor's dual list is:" << endl;

        forAll(pointFieldList, pfI)
        {
            Info<< "    Proc " << pfI << " has " << pointFieldList[pfI].size()
                << " duals" << endl;
        }

        label nNonMatches = 0;
        label nNearMatches = 0;
        label nExactMatches = 0;

        forAll(pointFieldList[0], pI)
        {
            const Point& masterPoint = pointFieldList[0][pI];

            bool foundMatch = false;
            bool foundNearMatch = false;

            scalar minCloseness = GREAT;
            Point closestPoint(0, 0, 0);

            forAll(pointFieldList, pfI)
            {
                if (pfI == 0)
                {
                    continue;
                }

//                label pfI = 1;

                forAll(pointFieldList[pfI], pISlave)
                {
                    const Point& slavePoint
                        = pointFieldList[pfI][pISlave];

                    if (masterPoint == slavePoint)
                    {
                        foundMatch = true;
                        break;
                    }

                    const scalar closeness = mag
                    (
                        topoint(masterPoint) - topoint(slavePoint)
                    );

                    if (closeness < 1e-12)
                    {
                        foundNearMatch = true;
                    }
                    else
                    {
                        if (closeness < minCloseness)
                        {
                            minCloseness = closeness;
                            closestPoint = slavePoint;
                        }
                    }
                }

                if (!foundMatch)
                {
                    if (foundNearMatch)
                    {
                        CGAL::Gmpq x(CGAL::to_double(masterPoint.x()));
                        CGAL::Gmpq y(CGAL::to_double(masterPoint.y()));
                        CGAL::Gmpq z(CGAL::to_double(masterPoint.z()));

                        std::cout<< "master = " << x << " " << y << " " << z
                            << std::endl;

                        CGAL::Gmpq xs(CGAL::to_double(closestPoint.x()));
                        CGAL::Gmpq ys(CGAL::to_double(closestPoint.y()));
                        CGAL::Gmpq zs(CGAL::to_double(closestPoint.z()));
                        std::cout<< "slave  = " << xs << " " << ys << " " << zs
                            << std::endl;

                        nNearMatches++;
                    }
                    else
                    {
                        nNonMatches++;
                        Info<< "    Closest point to " << masterPoint << " is "
                            << closestPoint << nl
                            << "    Separation is " << minCloseness << endl;

                        CGAL::Gmpq x(CGAL::to_double(masterPoint.x()));
                        CGAL::Gmpq y(CGAL::to_double(masterPoint.y()));
                        CGAL::Gmpq z(CGAL::to_double(masterPoint.z()));

                        std::cout<< "master = " << x << " " << y << " " << z
                                 << std::endl;

                        CGAL::Gmpq xs(CGAL::to_double(closestPoint.x()));
                        CGAL::Gmpq ys(CGAL::to_double(closestPoint.y()));
                        CGAL::Gmpq zs(CGAL::to_double(closestPoint.z()));
                        std::cout<< "slave  = " << xs << " " << ys << " " << zs
                                 << std::endl;
                    }
                }
                else
                {
                    nExactMatches++;
                }
            }
        }

        Info<< "Found " << nNonMatches << " non-matching duals" << nl
            << " and " << nNearMatches << " near matches"
            << " and " << nExactMatches << " exact matches" << endl;
    }
}


void Foam::conformalVoronoiMesh::checkVertices()
{
    List<pointField> pointFieldList(Pstream::nProcs());

    pointField points(number_of_vertices());

    labelPairHashSet duplicateVertices;

    label count = 0;
    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (duplicateVertices.found(labelPair(vit->procIndex(), vit->index())))
        {
            Pout<< "DUPLICATE " << vit->procIndex() << vit->index() << endl;
        }
        else
        {
            duplicateVertices.insert(labelPair(vit->procIndex(), vit->index()));
        }

        points[count++] = topoint(vit->point());
    }

    pointFieldList[Pstream::myProcNo()] = points;

    Pstream::gatherList(pointFieldList);

    OFstream str("missingPoints.obj");

    if (Pstream::master())
    {
        Info<< "Checking on master processor the point locations of each " << nl
            << "processor point list against the master point list." << nl
            << "There are " << pointFieldList.size() << " processors" << nl
            << "The size of each processor's point list is:" << endl;

        forAll(pointFieldList, pfI)
        {
            Info<< "    Proc " << pfI << " has " << pointFieldList[pfI].size()
                << " points" << endl;
        }

        label nNonMatches = 0;

        forAll(pointFieldList[0], pI)
        {
            const Foam::point& masterPoint = pointFieldList[0][pI];

            forAll(pointFieldList, pfI)
            {
                if (pI == 0)
                {
                    continue;
                }

                bool foundMatch = false;

                forAll(pointFieldList[pfI], pISlave)
                {
                    const Foam::point& slavePoint
                        = pointFieldList[pfI][pISlave];

                    if (masterPoint == slavePoint)
                    {
                        foundMatch = true;
                        break;
                    }
                }

                if (!foundMatch)
                {
                    Info<< "    Proc " << pfI << " Master != Slave -> "
                        << masterPoint << endl;

                    meshTools::writeOBJ(str, masterPoint);

                    nNonMatches++;
                }
            }
        }

        Info<< "Found a total of " << nNonMatches << " non-matching points"
            << endl;
    }
}


void Foam::conformalVoronoiMesh::calcDualMesh
(
    pointField& points,
    labelList& boundaryPts,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchTypes,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    labelList& procNeighbours,
    pointField& cellCentres,
    labelList& cellToDelaunayVertex,
    labelListList& patchToDelaunayVertex,
    PackedBoolList& boundaryFacesToRemove
)
{
    timeCheck("Start calcDualMesh");

//    if (debug)
//    {
//        Pout<< nl << "Perfoming some checks . . ." << nl << nl
//            << "Total number of vertices = " << number_of_vertices() << nl
//            << "Total number of cells    = " << number_of_finite_cells()
//            << endl;
//
//        checkVertices();
//        checkCells();
//        checkDuals();
//
//        Info<< nl << "Finished checks" << nl << endl;
//    }

    setVertexSizeAndAlignment();

    timeCheck("After setVertexSizeAndAlignment");

    indexDualVertices(points, boundaryPts);

    {
        Info<< nl << "Merging identical points" << endl;

        // There is no guarantee that a merge of close points is no-risk
        mergeIdenticalDualVertices(points, boundaryPts);
    }

    // Final dual face and owner neighbour construction

    timeCheck("Before createFacesOwnerNeighbourAndPatches");

    createFacesOwnerNeighbourAndPatches
    (
        faces,
        owner,
        neighbour,
        patchTypes,
        patchNames,
        patchSizes,
        patchStarts,
        procNeighbours,
        patchToDelaunayVertex,  // from patch face to Delaunay vertex (slavePp)
        boundaryFacesToRemove,
        false
    );

    // deferredCollapseFaceSet(owner, neighbour, deferredCollapseFaces);

    cellCentres = allPoints();

    cellToDelaunayVertex = removeUnusedCells(owner, neighbour);

    cellCentres = pointField(cellCentres, cellToDelaunayVertex);

    removeUnusedPoints(faces, points, boundaryPts);

    timeCheck("End of calcDualMesh");
}


void Foam::conformalVoronoiMesh::calcTetMesh
(
    pointField& points,
    labelList& pointToDelaunayVertex,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchTypes,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    labelList vertexMap(number_of_vertices());

    label vertI = 0;

    points.setSize(number_of_vertices());
    pointToDelaunayVertex.setSize(number_of_vertices());

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() || vit->boundaryPoint())
        {
            vertexMap[vit->index()] = vertI;
            points[vertI] = topoint(vit->point());
            pointToDelaunayVertex[vertI] = vit->index();
            vertI++;
        }
    }

    points.setSize(vertI);
    pointToDelaunayVertex.setSize(vertI);

    label cellI = 0;

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (cit->internalOrBoundaryDualVertex())
        {
             cit->cellIndex() = cellI++;
        }
        else
        {
            cit->cellIndex() = Cb::ctFar;
        }
    }

    patchNames = geometryToConformTo_.patchNames();

    patchNames.setSize(patchNames.size() + 1);

    patchNames[patchNames.size() - 1] = "cvMesh_defaultPatch";
    patchTypes.setSize(patchNames.size(), wallPolyPatch::typeName);

    label nPatches = patchNames.size();

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_finite_facets());

    owner.setSize(number_of_finite_facets());

    neighbour.setSize(number_of_finite_facets());

    label faceI = 0;

    labelList verticesOnTriFace(3, -1);

    face newFace(verticesOnTriFace);

    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const int oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        if (c1->hasFarPoint() && c2->hasFarPoint())
        {
            // Both tets are outside, skip
            continue;
        }

        label c1I = c1->cellIndex();
        label c2I = c2->cellIndex();

        label ownerCell = -1;
        label neighbourCell = -1;

        for (label i = 0; i < 3; i++)
        {
            verticesOnTriFace[i] = vertexMap
            [
                c1->vertex(vertex_triple_index(oppositeVertex, i))->index()
            ];
        }

        newFace = face(verticesOnTriFace);

        if (c1->hasFarPoint() || c2->hasFarPoint())
        {
            // Boundary face...
            if (c1->hasFarPoint())
            {
                //... with c1 outside
                ownerCell = c2I;
            }
            else
            {
                // ... with c2 outside
                ownerCell = c1I;

                reverse(newFace);
            }

            label patchIndex = geometryToConformTo_.findPatch
            (
                newFace.centre(points)
            );

            if (patchIndex == -1)
            {
                patchIndex = patchNames.size() - 1;

                WarningIn("Foam::conformalVoronoiMesh::calcTetMesh")
                    << "Tet face centre at  " << nl
                    << newFace.centre(points) << nl
                    << "did not find a surface patch. Adding to "
                    << patchNames[patchIndex]
                    << endl;
            }

            patchFaces[patchIndex].append(newFace);
            patchOwners[patchIndex].append(ownerCell);
        }
        else
        {
            // Internal face...
            if (c1I < c2I)
            {
                // ...with c1 as the ownerCell
                ownerCell = c1I;
                neighbourCell = c2I;

                reverse(newFace);
            }
            else
            {
                // ...with c2 as the ownerCell
                ownerCell = c2I;
                neighbourCell = c1I;
            }

            faces[faceI] = newFace;
            owner[faceI] = ownerCell;
            neighbour[faceI] = neighbourCell;
            faceI++;
        }
    }

    label nInternalFaces = faceI;

    faces.setSize(nInternalFaces);
    owner.setSize(nInternalFaces);
    neighbour.setSize(nInternalFaces);

    sortFaces(faces, owner, neighbour);

//    addPatches
//    (
//        nInternalFaces,
//        faces,
//        owner,
//        patchSizes,
//        patchStarts,
//        patchFaces,
//        patchOwners
//    );
}


void Foam::conformalVoronoiMesh::mergeIdenticalDualVertices
(
    const pointField& pts,
    const labelList& boundaryPts
)
{
    // Assess close points to be merged

    label nPtsMerged = 0;
    label nPtsMergedSum = 0;

    do
    {
        Map<label> dualPtIndexMap;

        nPtsMerged = mergeIdenticalDualVertices
        (
            pts,
            boundaryPts,
            dualPtIndexMap
        );

        reindexDualVertices(dualPtIndexMap);

        reduce(nPtsMerged, sumOp<label>());

        nPtsMergedSum += nPtsMerged;

    } while (nPtsMerged > 0);

    if (nPtsMergedSum > 0)
    {
        Info<< "    Merged " << nPtsMergedSum << " points " << endl;
    }
}


Foam::label Foam::conformalVoronoiMesh::mergeIdenticalDualVertices
(
    const pointField& pts,
    const labelList& boundaryPts,
    Map<label>& dualPtIndexMap
) const
{
    label nPtsMerged = 0;

    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const int oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        if (is_infinite(c1) || is_infinite(c2))
        {
            continue;
        }

        label& c1I = c1->cellIndex();
        label& c2I = c2->cellIndex();

        if ((c1I != c2I) && !c1->hasFarPoint() && !c2->hasFarPoint())
        {
            const Foam::point& p1 = pts[c1I];
            const Foam::point& p2 = pts[c2I];

            if (p1 == p2)
            {
                if (c1I < c2I)
                {
                    dualPtIndexMap.insert(c1I, c1I);
                    dualPtIndexMap.insert(c2I, c1I);
                }
                else
                {
                    dualPtIndexMap.insert(c1I, c2I);
                    dualPtIndexMap.insert(c2I, c2I);
                }

                nPtsMerged++;
            }
        }
    }

    if (debug)
    {
        Info<< "mergeIdenticalDualVertices:" << endl
            << "    zero-length edges     : "
            << returnReduce(nPtsMerged, sumOp<label>()) << endl
            << endl;
    }

    return nPtsMerged;
}


//void Foam::conformalVoronoiMesh::smoothSurface
//(
//    pointField& pts,
//    const labelList& boundaryPts
//)
//{
//    label nCollapsedFaces = 0;
//
//    label iterI = 0;
//
//    do
//    {
//        Map<label> dualPtIndexMap;
//
//        nCollapsedFaces = smoothSurfaceDualFaces
//        (
//            pts,
//            boundaryPts,
//            dualPtIndexMap
//        );
//
//        reduce(nCollapsedFaces, sumOp<label>());
//
//        reindexDualVertices(dualPtIndexMap);
//
//        mergeIdenticalDualVertices(pts, boundaryPts);
//
//        if (nCollapsedFaces > 0)
//        {
//            Info<< "    Collapsed " << nCollapsedFaces << " boundary faces"
//                << endl;
//        }
//
//        if (++iterI > cvMeshControls().maxCollapseIterations())
//        {
//            Info<< "    maxCollapseIterations reached, stopping collapse"
//                << endl;
//
//            break;
//        }
//
//    } while (nCollapsedFaces > 0);
//
//    // Force all points of boundary faces to be on the surface
////    for
////    (
////        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
////        cit != finite_cells_end();
////        ++cit
////    )
////    {
////        label ptI = cit->cellIndex();
////
////        label fC = cit->filterCount();
////
////        if (fC > cvMeshControls().filterCountSkipThreshold())
////        {
////            // This vertex has been limited too many times, skip
////            continue;
////        }
////
////        // Only cells with indices > -1 are valid
////        if (ptI > -1)
////        {
////            if (boundaryPts[ptI] != -1)
////            {
////                Foam::point& pt = pts[ptI];
////
////                pointIndexHit surfHit;
////                label hitSurface;
////
////                geometryToConformTo_.findSurfaceNearest
////                (
////                    pt,
////                    sqr(GREAT),
////                    surfHit,
////                    hitSurface
////                );
////
////                if (surfHit.hit())
////                {
////                    pt += (surfHit.hitPoint() - pt)
////                         *pow
////                          (
////                              cvMeshControls().filterErrorReductionCoeff(),
////                              fC
////                          );
////                }
////            }
////        }
////    }
////
////    mergeCloseDualVertices(pts, boundaryPts);
//}
//
//
//Foam::label Foam::conformalVoronoiMesh::smoothSurfaceDualFaces
//(
//    pointField& pts,
//    const labelList& boundaryPts,
//    Map<label>& dualPtIndexMap
//) const
//{
//    label nCollapsedFaces = 0;
//
//    const scalar cosPerpendicularToleranceAngle = cos
//    (
//        degToRad(cvMeshControls().surfaceStepFaceAngle())
//    );
//
//    for
//    (
//        Delaunay::Finite_edges_iterator eit = finite_edges_begin();
//        eit != finite_edges_end();
//        ++eit
//    )
//    {
//        Cell_circulator ccStart = incident_cells(*eit);
//        Cell_circulator cc = ccStart;
//
//        bool skipFace = false;
//
//        do
//        {
//            if (dualPtIndexMap.found(cc->cellIndex()))
//            {
//                // One of the points of this face has already been
//                // collapsed this sweep, leave for next sweep
//
//                skipFace = true;
//
//                break;
//            }
//
//        } while (++cc != ccStart);
//
//        if (skipFace)
//        {
//            continue;
//        }
//
//        if (isBoundaryDualFace(eit))
//        {
//            face dualFace = buildDualFace(eit);
//
//            if (dualFace.size() < 3)
//            {
//                // This face has been collapsed already
//                continue;
//            }
//
//            label maxFC = maxFilterCount(eit);
//
//            if (maxFC > cvMeshControls().filterCountSkipThreshold())
//            {
//                // A vertex on this face has been limited too many
//                // times, skip
//                continue;
//            }
//
//            Cell_handle c = eit->first;
//            Vertex_handle vA = c->vertex(eit->second);
//            Vertex_handle vB = c->vertex(eit->third);
//
//            if
//            (
//                vA->internalBoundaryPoint() && vA->surfacePoint()
//             && vB->externalBoundaryPoint() && vB->surfacePoint()
//            )
//            {
//                if (vA->index() == vB->index() - 1)
//                {
//                    continue;
//                }
//            }
//            else if
//            (
//                vA->externalBoundaryPoint() && vA->surfacePoint()
//             && vB->internalBoundaryPoint() && vB->surfacePoint()
//            )
//            {
//                if (vA->index() == vB->index() + 1)
//                {
//                    continue;
//                }
//            }
////            else if
////            (
////                vA->internalBoundaryPoint() && vA->featureEdgePoint()
////             && vB->externalBoundaryPoint() && vB->featureEdgePoint()
////            )
////            {
////                if (vA->index() == vB->index() - 1)
////                {
////                    continue;
////                }
////            }
////            else if
////            (
////                vA->externalBoundaryPoint() && vA->featureEdgePoint()
////             && vB->internalBoundaryPoint() && vB->featureEdgePoint()
////            )
////            {
////                if (vA->index() == vB->index() + 1)
////                {
////                    continue;
////                }
////            }
////            else if
////            (
////                vA->internalBoundaryPoint() && vA->featurePoint()
////             && vB->externalBoundaryPoint() && vB->featurePoint()
////            )
////            {
////                if (vA->index() == vB->index() - 1)
////                {
////                    continue;
////                }
////            }
////            else if
////            (
////                vA->externalBoundaryPoint() && vA->featurePoint()
////             && vB->internalBoundaryPoint() && vB->featurePoint()
////            )
////            {
////                if (vA->index() == vB->index() + 1)
////                {
////                    continue;
////                }
////            }
//
//
//            if ((faceNormal & surfaceNormal) < cosPerpendicularToleranceAngle)
//            {
//                scalar targetFaceSize = averageAnyCellSize(vA, vB);
//
//                // Selecting faces to collapse based on angle to
//                // surface, so set collapseSizeLimitCoeff to GREAT to
//                // allow collapse of all faces
//
//                faceCollapseMode mode = collapseFace
//                (
//                    dualFace,
//                    pts,
//                    boundaryPts,
//                    dualPtIndexMap,
//                    targetFaceSize,
//                    GREAT,
//                    maxFC
//                );
//
//                if (mode == fcmPoint || mode == fcmEdge)
//                {
//                    nCollapsedFaces++;
//                }
//            }
//        }
//    }
//
//    return nCollapsedFaces;
//}


void Foam::conformalVoronoiMesh::deferredCollapseFaceSet
(
    labelList& owner,
    labelList& neighbour,
    const HashSet<labelPair, labelPair::Hash<> >& deferredCollapseFaces
) const
{
    DynamicList<label> faceLabels;

    forAll(neighbour, nI)
    {
        if (deferredCollapseFaces.found(Pair<label>(owner[nI], neighbour[nI])))
        {
            faceLabels.append(nI);
        }
    }

    Pout<< "facesToCollapse" << nl << faceLabels << endl;
}


Foam::autoPtr<Foam::polyMesh>
Foam::conformalVoronoiMesh::createPolyMeshFromPoints
(
    const pointField& pts
) const
{
    faceList faces;
    labelList owner;
    labelList neighbour;
    wordList patchTypes;
    wordList patchNames;
    labelList patchSizes;
    labelList patchStarts;
    labelList procNeighbours;
    pointField cellCentres;
    labelListList patchToDelaunayVertex;
    PackedBoolList boundaryFacesToRemove;

    timeCheck("Start of checkPolyMeshQuality");

    Info<< nl << "Creating polyMesh to assess quality" << endl;

    createFacesOwnerNeighbourAndPatches
    (
        faces,
        owner,
        neighbour,
        patchTypes,
        patchNames,
        patchSizes,
        patchStarts,
        procNeighbours,
        patchToDelaunayVertex,
        boundaryFacesToRemove,
        false
    );

    //createCellCentres(cellCentres);
    cellCentres = allPoints();

    labelList cellToDelaunayVertex(removeUnusedCells(owner, neighbour));
    cellCentres = pointField(cellCentres, cellToDelaunayVertex);

    autoPtr<polyMesh> meshPtr
    (
        new polyMesh
        (
            IOobject
            (
                "cvMesh_temporary",
                runTime_.timeName(),
                runTime_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferCopy(pts),
            xferMove(faces),
            xferMove(owner),
            xferMove(neighbour)
        )
    );

    polyMesh& pMesh = meshPtr();

    List<polyPatch*> patches(patchStarts.size());

    label nValidPatches = 0;

    forAll(patches, p)
    {
       if (patchTypes[p] == processorPolyPatch::typeName)
       {
           // Do not create empty processor patches
           if (patchSizes[p] > 0)
           {
               patches[nValidPatches] = new processorPolyPatch
               (
                   patchNames[p],
                   patchSizes[p],
                   patchStarts[p],
                   nValidPatches,
                   pMesh.boundaryMesh(),
                   Pstream::myProcNo(),
                   procNeighbours[p],
                   coupledPolyPatch::COINCIDENTFULLMATCH
               );

               nValidPatches++;
           }
       }
       else
       {
           patches[nValidPatches] = polyPatch::New
           (
               patchTypes[p],
               patchNames[p],
               patchSizes[p],
               patchStarts[p],
               nValidPatches,
               pMesh.boundaryMesh()
           ).ptr();

           nValidPatches++;
       }
    }

    patches.setSize(nValidPatches);

    pMesh.addPatches(patches);

    // Info<< "ADDPATCHES NOT IN PARALLEL" << endl;

    // forAll(patches, p)
    // {
    //     patches[p] = new polyPatch
    //     (
    //         patchNames[p],
    //         patchSizes[p],
    //         patchStarts[p],
    //         p,
    //         pMesh.boundaryMesh()
    //     );
    // }

    // pMesh.addPatches(patches, false);

    // pMesh.overrideCellCentres(cellCentres);

    return meshPtr;
}


void Foam::conformalVoronoiMesh::checkCellSizing()
{
    Info<< "Checking cell sizes..."<< endl;

    timeCheck("Start of Cell Sizing");

    labelList boundaryPts(number_of_finite_cells(), -1);
    pointField ptsField;

    indexDualVertices(ptsField, boundaryPts);

    // Merge close dual vertices.
    mergeIdenticalDualVertices(ptsField, boundaryPts);

    autoPtr<polyMesh> meshPtr = createPolyMeshFromPoints(ptsField);
    const polyMesh& pMesh = meshPtr();

    //pMesh.write();

    // Find cells with poor quality
    DynamicList<label> checkFaces(identity(pMesh.nFaces()));
    labelHashSet wrongFaces(pMesh.nFaces()/100);

    Info<< "Running checkMesh on mesh with " << pMesh.nCells()
        << " cells "<< endl;

    const dictionary& dict
        = cvMeshControls().cvMeshDict().subDict("meshQualityControls");

    const scalar maxNonOrtho = readScalar(dict.lookup("maxNonOrtho", true));

    label nWrongFaces = 0;

    if (maxNonOrtho < 180.0 - SMALL)
    {
        polyMeshGeometry::checkFaceDotProduct
        (
            false,
            maxNonOrtho,
            pMesh,
            pMesh.cellCentres(),
            pMesh.faceAreas(),
            checkFaces,
            List<labelPair>(),
            &wrongFaces
        );

        label nNonOrthogonal = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    non-orthogonality > " << maxNonOrtho
            << " degrees : " << nNonOrthogonal << endl;

        nWrongFaces += nNonOrthogonal;
    }

    labelHashSet protrudingCells = findOffsetPatchFaces(pMesh, 0.25);

    label nProtrudingCells = protrudingCells.size();

    Info<< "    protruding/intruding cells : " << nProtrudingCells << endl;

    nWrongFaces += nProtrudingCells;

//    motionSmoother::checkMesh
//    (
//        false,
//        pMesh,
//        cvMeshControls().cvMeshDict().subDict("meshQualityControls"),
//        checkFaces,
//        wrongFaces
//    );

    Info<< "    Found total of " << nWrongFaces << " bad faces" << endl;

    {
        labelHashSet cellsToResizeMap(pMesh.nFaces()/100);

        // Find cells that are attached to the faces in wrongFaces.
        forAllConstIter(labelHashSet, wrongFaces, iter)
        {
            const label faceOwner = pMesh.faceOwner()[iter.key()];
            const label faceNeighbour = pMesh.faceNeighbour()[iter.key()];

            if (!cellsToResizeMap.found(faceOwner))
            {
                cellsToResizeMap.insert(faceOwner);
            }

            if (!cellsToResizeMap.found(faceNeighbour))
            {
                cellsToResizeMap.insert(faceNeighbour);
            }
        }

        cellsToResizeMap += protrudingCells;

        pointField cellsToResize(cellsToResizeMap.size());

        label count = 0;
        for (label cellI = 0; cellI < pMesh.nCells(); ++cellI)
        {
            if (cellsToResizeMap.found(cellI))
            {
                cellsToResize[count++] = pMesh.cellCentres()[cellI];
            }
        }

        Info<< "    DISABLED: Automatically re-sizing " << cellsToResize.size()
            << " cells that are attached to the bad faces: " << endl;

        //cellSizeControl_.setCellSizes(cellsToResize);
    }

    timeCheck("End of Cell Sizing");

    Info<< "Finished checking cell sizes"<< endl;
}


Foam::labelHashSet Foam::conformalVoronoiMesh::findOffsetPatchFaces
(
    const polyMesh& mesh,
    const scalar allowedOffset
) const
{
    timeCheck("Start findRemainingProtrusionSet");

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    cellSet offsetBoundaryCells
    (
        mesh,
        "cvMesh_protrudingCells",
        mesh.nCells()/1000
    );

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        const faceList& localFaces = patch.localFaces();
        const pointField& localPoints = patch.localPoints();

        const labelList& fCell = patch.faceCells();

        forAll(localFaces, pLFI)
        {
            const face& f = localFaces[pLFI];

            const Foam::point& faceCentre = f.centre(localPoints);

            const scalar targetSize = targetCellSize(faceCentre);

            pointIndexHit pHit;
            label surfHit = -1;

            geometryToConformTo_.findSurfaceNearest
            (
                faceCentre,
                sqr(targetSize),
                pHit,
                surfHit
            );

            if
            (
                pHit.hit()
             && (mag(pHit.hitPoint() - faceCentre) > allowedOffset*targetSize)
            )
            {
                offsetBoundaryCells.insert(fCell[pLFI]);
            }
        }
    }

    if (cvMeshControls().objOutput())
    {
        offsetBoundaryCells.write();
    }

    return offsetBoundaryCells;
}


Foam::labelHashSet Foam::conformalVoronoiMesh::checkPolyMeshQuality
(
    const pointField& pts
) const
{
    autoPtr<polyMesh> meshPtr = createPolyMeshFromPoints(pts);
    polyMesh& pMesh = meshPtr();

    timeCheck("polyMesh created, checking quality");

    labelHashSet wrongFaces(pMesh.nFaces()/100);

    DynamicList<label> checkFaces(pMesh.nFaces());

    const vectorField& fAreas = pMesh.faceAreas();

    scalar faceAreaLimit = SMALL;

    forAll(fAreas, fI)
    {
        if (mag(fAreas[fI]) > faceAreaLimit)
        {
            checkFaces.append(fI);
        }
    }

    Info<< nl << "Excluding "
        << returnReduce(fAreas.size() - checkFaces.size(), sumOp<label>())
        << " faces from check, < " << faceAreaLimit << " area" << endl;

    motionSmoother::checkMesh
    (
        false,
        pMesh,
        cvMeshControls().cvMeshDict().subDict("meshQualityControls"),
        checkFaces,
        wrongFaces
    );

    {
        // Check for cells with more than 1 but fewer than 4 faces
        label nInvalidPolyhedra = 0;

        const cellList& cells = pMesh.cells();

        forAll(cells, cI)
        {
            if (cells[cI].size() < 4 && cells[cI].size() > 0)
            {
                // Pout<< "cell " << cI << " " << cells[cI]
                //     << " has " << cells[cI].size() << " faces."
                //     << endl;

                nInvalidPolyhedra++;

                forAll(cells[cI], cFI)
                {
                    wrongFaces.insert(cells[cI][cFI]);
                }
            }
        }

        Info<< "    cells with more than 1 but fewer than 4 faces          : "
            << returnReduce(nInvalidPolyhedra, sumOp<label>())
            << endl;

        // Check for cells with one internal face only

        labelList nInternalFaces(pMesh.nCells(), 0);

        for (label fI = 0; fI < pMesh.nInternalFaces(); fI++)
        {
            nInternalFaces[pMesh.faceOwner()[fI]]++;
            nInternalFaces[pMesh.faceNeighbour()[fI]]++;
        }

        const polyBoundaryMesh& patches = pMesh.boundaryMesh();

        forAll(patches, patchI)
        {
            if (patches[patchI].coupled())
            {
                const labelUList& owners = patches[patchI].faceCells();

                forAll(owners, i)
                {
                    nInternalFaces[owners[i]]++;
                }
            }
        }

        label oneInternalFaceCells = 0;

        forAll(nInternalFaces, cI)
        {
            if (nInternalFaces[cI] <= 1)
            {
                oneInternalFaceCells++;

                forAll(cells[cI], cFI)
                {
                    wrongFaces.insert(cells[cI][cFI]);
                }
            }
        }

        Info<< "    cells with with zero or one non-boundary face          : "
            << returnReduce(oneInternalFaceCells, sumOp<label>())
            << endl;
    }


    PackedBoolList ptToBeLimited(pts.size(), false);

    forAllConstIter(labelHashSet, wrongFaces, iter)
    {
        const face f = pMesh.faces()[iter.key()];

        forAll(f, fPtI)
        {
            ptToBeLimited[f[fPtI]] = true;
        }
    }

    // // Limit connected cells

    // labelHashSet limitCells(pMesh.nCells()/100);

    // const labelListList& ptCells = pMesh.pointCells();

    // forAllConstIter(labelHashSet, wrongFaces, iter)
    // {
    //     const face f = pMesh.faces()[iter.key()];

    //     forAll(f, fPtI)
    //     {
    //         label ptI = f[fPtI];

    //         const labelList& pC = ptCells[ptI];

    //         forAll(pC, pCI)
    //         {
    //             limitCells.insert(pC[pCI]);
    //         }
    //     }
    // }

    // const labelListList& cellPts = pMesh.cellPoints();

    // forAllConstIter(labelHashSet, limitCells, iter)
    // {
    //     label cellI = iter.key();

    //     const labelList& cP = cellPts[cellI];

    //     forAll(cP, cPI)
    //     {
    //         ptToBeLimited[cP[cPI]] = true;
    //     }
    // }


    // Apply Delaunay cell filterCounts and determine the maximum
    // overall filterCount

    label maxFilterCount = 0;

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        label cI = cit->cellIndex();

        if (cI >= 0)
        {
            if (ptToBeLimited[cI] == true)
            {
                cit->filterCount()++;
            }

            if (cit->filterCount() > maxFilterCount)
            {
                maxFilterCount = cit->filterCount();
            }
        }
    }

    Info<< nl << "Maximum number of filter limits applied: "
        << returnReduce(maxFilterCount, maxOp<label>()) << endl;

    return wrongFaces;
}


void Foam::conformalVoronoiMesh::indexDualVertices
(
    pointField& pts,
    labelList& boundaryPts
)
{
    // Indexing Delaunay cells, which are the dual vertices

    this->resetCellCount();

    pts.setSize(number_of_finite_cells());

    boundaryPts.setSize(number_of_finite_cells(), -1);

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
//        if (tetrahedron(cit).volume() == 0)
//        {
//            Pout<< "ZERO VOLUME TET" << endl;
//            Pout<< cit->info();
//            Pout<< "Dual = " << cit->dual();
//        }

        if (!cit->hasFarPoint())
        {
            cit->cellIndex() = getNewCellIndex();

            // For nearly coplanar Delaunay cells that are present on different
            // processors the result of the circumcentre calculation depends on
            // the ordering of the vertices, so synchronise it across processors

            if (Pstream::parRun() && cit->parallelDualVertex())
            {
                typedef CGAL::Exact_predicates_exact_constructions_kernel Exact;
                typedef CGAL::Point_3<Exact> ExactPoint;

                List<labelPair> cellVerticesPair(4);
                List<ExactPoint> cellVertices(4);

                for (label vI = 0; vI < 4; ++vI)
                {
                    cellVerticesPair[vI] = labelPair
                    (
                        cit->vertex(vI)->procIndex(),
                        cit->vertex(vI)->index()
                    );

                    cellVertices[vI] = ExactPoint
                    (
                        cit->vertex(vI)->point().x(),
                        cit->vertex(vI)->point().y(),
                        cit->vertex(vI)->point().z()
                    );
                }

                // Sort the vertices so that they will be in the same order on
                // each processor
                labelList oldToNew;
                sortedOrder(cellVerticesPair, oldToNew);
                oldToNew = invert(oldToNew.size(), oldToNew);
                inplaceReorder(oldToNew, cellVertices);

                ExactPoint synchronisedDual = CGAL::circumcenter
                (
                    cellVertices[0],
                    cellVertices[1],
                    cellVertices[2],
                    cellVertices[3]
                );

                pts[cit->cellIndex()] = Foam::point
                (
                    CGAL::to_double(synchronisedDual.x()),
                    CGAL::to_double(synchronisedDual.y()),
                    CGAL::to_double(synchronisedDual.z())
                );
            }
            else
            {
                pts[cit->cellIndex()] = cit->dual();
            }

            if (cit->boundaryDualVertex())
            {
                if (cit->featureEdgeDualVertex())
                {
                    boundaryPts[cit->cellIndex()] = 1;
                }
                else
                {
                    boundaryPts[cit->cellIndex()] = 0;
                }
            }
        }
        else
        {
            cit->cellIndex() = Cb::ctFar;
        }
    }

    pts.setSize(this->cellCount());

    boundaryPts.setSize(this->cellCount());
}


void Foam::conformalVoronoiMesh::reindexDualVertices
(
    const Map<label>& dualPtIndexMap
)
{
    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (dualPtIndexMap.found(cit->cellIndex()))
        {
            cit->cellIndex() = dualPtIndexMap[cit->cellIndex()];
        }
    }
}


Foam::label Foam::conformalVoronoiMesh::createPatchInfo
(
    wordList& patchNames,
    wordList& patchTypes,
    labelList& procNeighbours
) const
{
    patchNames = geometryToConformTo_.patchNames();
    patchTypes.setSize(patchNames.size() + 1, wallPolyPatch::typeName);
    procNeighbours.setSize(patchNames.size() + 1, -1);

    const PtrList<dictionary>& patchInfo = geometryToConformTo_.patchInfo();

    forAll(patchNames, patchI)
    {
        if (patchInfo.set(patchI))
        {
            patchTypes[patchI] =
                patchInfo[patchI].lookupOrDefault<word>
                (
                    "type",
                    wallPolyPatch::typeName
                );
        }
    }

    patchNames.setSize(patchNames.size() + 1);
    label defaultPatchIndex = patchNames.size() - 1;
    patchNames[defaultPatchIndex] = "cvMesh_defaultPatch";

    label nProcPatches = 0;

    if (Pstream::parRun())
    {
        List<boolList> procUsedList
        (
            Pstream::nProcs(),
            boolList(Pstream::nProcs(), false)
        );

        boolList& procUsed = procUsedList[Pstream::myProcNo()];

        // Determine which processor patches are required
        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            // This test is not sufficient if one of the processors does
            // not receive a referred vertex from another processor, but does
            // send one to the other processor.
            if (vit->referred())
            {
                procUsed[vit->procIndex()] = true;
            }
        }

        // Because the previous test was insufficient, combine the lists.
        Pstream::gatherList(procUsedList);
        Pstream::scatterList(procUsedList);

        forAll(procUsedList, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                if (procUsedList[procI][Pstream::myProcNo()])
                {
                    procUsed[procI] = true;
                }
            }
        }

        forAll(procUsed, pUI)
        {
            if (procUsed[pUI])
            {
                nProcPatches++;
            }
        }

        label nNonProcPatches = patchNames.size();

        patchNames.setSize(nNonProcPatches + nProcPatches);
        patchTypes.setSize(nNonProcPatches + nProcPatches);
        procNeighbours.setSize(nNonProcPatches + nProcPatches, -1);

        label procAddI = 0;

        forAll(procUsed, pUI)
        {
            if (procUsed[pUI])
            {
                patchTypes[nNonProcPatches + procAddI] =
                    processorPolyPatch::typeName;

                patchNames[nNonProcPatches + procAddI] =
                    "procBoundary"
                   + name(Pstream::myProcNo())
                   + "to"
                   + name(pUI);

                procNeighbours[nNonProcPatches + procAddI] = pUI;

                procAddI++;
            }
        }
    }

    return defaultPatchIndex;
}


void Foam::conformalVoronoiMesh::createFacesOwnerNeighbourAndPatches
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchTypes,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    labelList& procNeighbours,
    labelListList& patchPointPairSlaves,
    PackedBoolList& boundaryFacesToRemove,
    bool includeEmptyPatches
) const
{
    const label defaultPatchIndex = createPatchInfo
    (
        patchNames,
        patchTypes,
        procNeighbours
    );

    const label nPatches = patchNames.size();

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));
    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));
    // Per patch face the index of the slave node of the point pair
    List<DynamicList<label> > patchPPSlaves(nPatches, DynamicList<label>(0));

    List<DynamicList<bool> > indirectPatchFace(nPatches, DynamicList<bool>(0));

    faces.setSize(number_of_finite_edges());
    owner.setSize(number_of_finite_edges());
    neighbour.setSize(number_of_finite_edges());

    labelPairPairDynListList procPatchSortingIndex(nPatches);

    label dualFaceI = 0;

    for
    (
        Delaunay::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_handle c = eit->first;
        Vertex_handle vA = c->vertex(eit->second);
        Vertex_handle vB = c->vertex(eit->third);

        if
        (
            (vA->internalOrBoundaryPoint() && !vA->referred())
         || (vB->internalOrBoundaryPoint() && !vB->referred())
        )
        {
            face newDualFace = buildDualFace(eit);

            if (newDualFace.size() >= 3)
            {
                label own = -1;
                label nei = -1;

                if (ownerAndNeighbour(vA, vB, own, nei))
                {
                    reverse(newDualFace);
                }

                if (nei == -1)
                {
                    // boundary face

                    pointFromPoint ptA = topoint(vA->point());
                    pointFromPoint ptB = topoint(vB->point());

                    label patchIndex = -1;

                    if (isProcBoundaryEdge(eit))
                    {
                        // One (and only one) of the points is an internal
                        // point from another processor

                        label procIndex = max(vA->procIndex(), vB->procIndex());

                        patchIndex = max
                        (
                            findIndex(procNeighbours, vA->procIndex()),
                            findIndex(procNeighbours, vB->procIndex())
                        );

                        // The lower processor index is the owner of the
                        // two for the purpose of sorting the patch faces.

                        if (Pstream::myProcNo() < procIndex)
                        {
                            // Use this processor's vertex index as the master
                            // for sorting

                            DynamicList<Pair<labelPair> >& sortingIndex =
                                procPatchSortingIndex[patchIndex];

                            if (vB->internalOrBoundaryPoint() && vB->referred())
                            {
                                sortingIndex.append
                                (
                                    Pair<labelPair>
                                    (
                                        labelPair(vA->index(), vA->procIndex()),
                                        labelPair(vB->index(), vB->procIndex())
                                    )
                                );
                            }
                            else
                            {
                                sortingIndex.append
                                (
                                    Pair<labelPair>
                                    (
                                        labelPair(vB->index(), vB->procIndex()),
                                        labelPair(vA->index(), vA->procIndex())
                                    )
                                );
                            }
                        }
                        else
                        {
                            // Use the other processor's vertex index as the
                            // master for sorting

                            DynamicList<Pair<labelPair> >& sortingIndex =
                                procPatchSortingIndex[patchIndex];

                            if (vA->internalOrBoundaryPoint() && vA->referred())
                            {
                                sortingIndex.append
                                (
                                    Pair<labelPair>
                                    (
                                        labelPair(vA->index(), vA->procIndex()),
                                        labelPair(vB->index(), vB->procIndex())
                                    )
                                );
                            }
                            else
                            {
                                sortingIndex.append
                                (
                                    Pair<labelPair>
                                    (
                                        labelPair(vB->index(), vB->procIndex()),
                                        labelPair(vA->index(), vA->procIndex())
                                    )
                                );
                            }
                        }

//                        Pout<< ptA << " " << ptB
//                            << " proc indices "
//                            << vA->procIndex() << " " << vB->procIndex()
//                            << " indices " << vA->index()
//                            << " " << vB->index()
//                            << " my proc " << Pstream::myProcNo()
//                            << " addedIndex "
//                            << procPatchSortingIndex[patchIndex].last()
//                            << endl;
                    }
                    else
                    {
                        patchIndex = geometryToConformTo_.findPatch(ptA, ptB);
                    }

                    if (patchIndex == -1)
                    {
                        // Did not find a surface patch between
                        // between Dv pair, finding nearest patch

//                         Pout<< "Did not find a surface patch between "
//                             << "for face, finding nearest patch to"
//                             << 0.5*(ptA + ptB) << endl;

                        patchIndex = geometryToConformTo_.findPatch
                        (
                            0.5*(ptA + ptB)
                        );
                    }

                    patchFaces[patchIndex].append(newDualFace);
                    patchOwners[patchIndex].append(own);

                    // If the two vertices are a pair, then the patch face is
                    // a desired one.
                    if (!isPointPair(vA, vB))
                    {
                        indirectPatchFace[patchIndex].append(true);
                    }
                    else
                    {
                        indirectPatchFace[patchIndex].append(false);
                    }

                    // Store the non-internal or boundary point
                    if (vA->internalOrBoundaryPoint())
                    {
                        patchPPSlaves[patchIndex].append(vB->index());
                    }
                    else
                    {
                        patchPPSlaves[patchIndex].append(vA->index());
                    }
                }
                else
                {
                    // internal face
                    faces[dualFaceI] = newDualFace;
                    owner[dualFaceI] = own;
                    neighbour[dualFaceI] = nei;

                    dualFaceI++;
                }
            }
        }
    }

    if (!patchFaces[defaultPatchIndex].empty())
    {
        Pout<< nl << patchFaces[defaultPatchIndex].size()
            << " faces were not able to have their patch determined from "
            << "the surface. "
            << nl <<  "Adding to patch " << patchNames[defaultPatchIndex]
            << endl;
    }

    label nInternalFaces = dualFaceI;

    faces.setSize(nInternalFaces);
    owner.setSize(nInternalFaces);
    neighbour.setSize(nInternalFaces);

    timeCheck("polyMesh quality checked");

    sortFaces(faces, owner, neighbour);

    sortProcPatches
    (
        patchFaces,
        patchOwners,
        patchPPSlaves,
        procPatchSortingIndex
    );

    timeCheck("faces, owner, neighbour sorted");

    addPatches
    (
        nInternalFaces,
        faces,
        owner,
        patchSizes,
        patchStarts,
        boundaryFacesToRemove,
        patchFaces,
        patchOwners,
        indirectPatchFace
    );

    // Return     patchPointPairSlaves.setSize(nPatches);
    patchPointPairSlaves.setSize(nPatches);
    forAll(patchPPSlaves, patchI)
    {
        patchPointPairSlaves[patchI].transfer(patchPPSlaves[patchI]);
    }

//    if (cvMeshControls().objOutput())
    {
        Info<< "Writing processor interfaces" << endl;

        forAll(procNeighbours, nbI)
        {
            if (patchFaces[nbI].size() > 0)
            {
                const label neighbour = procNeighbours[nbI];

                faceList procPatchFaces = patchFaces[nbI];

                // Reverse faces as it makes it easier to analyse the output
                // using a diff
                if (neighbour < Pstream::myProcNo())
                {
                    forAll(procPatchFaces, fI)
                    {
                        procPatchFaces[fI] = procPatchFaces[fI].reverseFace();
                    }
                }

                if (neighbour != -1)
                {
                    word fName =
                        "processor_"
                      + name(Pstream::myProcNo())
                      + "_to_"
                      + name(neighbour)
                      + "_interface.obj";

                    writeProcessorInterface(fName, procPatchFaces);
                }
            }
        }
    }
}


void Foam::conformalVoronoiMesh::createCellCentres
(
    pointField& cellCentres
) const
{
    cellCentres.setSize(number_of_vertices(), point::max);

    label vertI = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            cellCentres[vertI++] = topoint(vit->point());
        }
    }

    cellCentres.setSize(vertI);
}


Foam::tmp<Foam::pointField> Foam::conformalVoronoiMesh::allPoints() const
{
    tmp<pointField> tpts(new pointField(number_of_vertices(), point::max));
    pointField& pts = tpts();

    label nVert = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            pts[nVert++] = topoint(vit->point());
        }
    }

    return tpts;
}


void Foam::conformalVoronoiMesh::sortFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour
) const
{
    // Upper triangular order:
    // + owner is sorted in ascending cell order
    // + within each block of equal value for owner, neighbour is sorted in
    //   ascending cell order.
    // + faces sorted to correspond
    // e.g.
    // owner | neighbour
    // 0     | 2
    // 0     | 23
    // 0     | 71
    // 1     | 23
    // 1     | 24
    // 1     | 91

    List<labelPair> ownerNeighbourPair(owner.size());

    forAll(ownerNeighbourPair, oNI)
    {
        ownerNeighbourPair[oNI] = labelPair(owner[oNI], neighbour[oNI]);
    }

    Info<< nl
        << "Sorting faces, owner and neighbour into upper triangular order"
        << endl;

    labelList oldToNew;

    sortedOrder(ownerNeighbourPair, oldToNew);

    oldToNew = invert(oldToNew.size(), oldToNew);

    inplaceReorder(oldToNew, faces);
    inplaceReorder(oldToNew, owner);
    inplaceReorder(oldToNew, neighbour);
}


void Foam::conformalVoronoiMesh::sortProcPatches
(
    List<DynamicList<face> >& patchFaces,
    List<DynamicList<label> >& patchOwners,
    List<DynamicList<label> >& patchPointPairSlaves,
    labelPairPairDynListList& patchSortingIndices
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    forAll(patchSortingIndices, patchI)
    {
        faceList& faces = patchFaces[patchI];
        labelList& owner = patchOwners[patchI];
        DynamicList<label>& slaves = patchPointPairSlaves[patchI];

        DynamicList<Pair<labelPair> >& sortingIndices
            = patchSortingIndices[patchI];

        if (!sortingIndices.empty())
        {
            if
            (
                faces.size() != sortingIndices.size()
             || owner.size() != sortingIndices.size()
             || slaves.size() != sortingIndices.size()
            )
            {
                FatalErrorIn
                (
                    "void Foam::conformalVoronoiMesh::sortProcPatches"
                    "("
                        "List<DynamicList<face> >& patchFaces, "
                        "List<DynamicList<label> >& patchOwners, "
                        "const List<DynamicList<label> >& patchSortingIndices"
                    ") const"
                )
                    << "patch size and size of sorting indices is inconsistent "
                    << " for patch " << patchI << nl
                    << " faces.size() " << faces.size() << nl
                    << " owner.size() " << owner.size() << nl
                    << " slaves.size() " << slaves.size() << nl
                    << " sortingIndices.size() "
                    << sortingIndices.size()
                    << exit(FatalError) << endl;
            }

            labelList oldToNew;

            sortedOrder(sortingIndices, oldToNew);

            oldToNew = invert(oldToNew.size(), oldToNew);

            inplaceReorder(oldToNew, sortingIndices);
            inplaceReorder(oldToNew, faces);
            inplaceReorder(oldToNew, owner);
            inplaceReorder(oldToNew, slaves);
        }
    }
}


void Foam::conformalVoronoiMesh::addPatches
(
    const label nInternalFaces,
    faceList& faces,
    labelList& owner,
    labelList& patchSizes,
    labelList& patchStarts,
    PackedBoolList& boundaryFacesToRemove,
    const List<DynamicList<face> >& patchFaces,
    const List<DynamicList<label> >& patchOwners,
    const List<DynamicList<bool> >& indirectPatchFace
) const
{
    label nPatches = patchFaces.size();

    patchSizes.setSize(nPatches, -1);
    patchStarts.setSize(nPatches, -1);

    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        patchSizes[p] = patchFaces[p].size();
        patchStarts[p] = nInternalFaces + nBoundaryFaces;

        nBoundaryFaces += patchSizes[p];
    }

    faces.setSize(nInternalFaces + nBoundaryFaces);
    owner.setSize(nInternalFaces + nBoundaryFaces);

    label faceI = nInternalFaces;

    forAll(patchFaces, p)
    {
        forAll(patchFaces[p], f)
        {
            faces[faceI] = patchFaces[p][f];
            owner[faceI] = patchOwners[p][f];
            boundaryFacesToRemove[faceI] = indirectPatchFace[p][f];

            faceI++;
        }
    }
}


void Foam::conformalVoronoiMesh::removeUnusedPoints
(
    faceList& faces,
    pointField& pts,
    labelList& boundaryPts
) const
{
    Info<< nl << "Removing unused points" << endl;

    PackedBoolList ptUsed(pts.size(), false);

    // Scan all faces to find all of the points that are used

    forAll(faces, fI)
    {
        const face& f = faces[fI];

        forAll(f, fPtI)
        {
            ptUsed[f[fPtI]] = true;
        }
    }

    label pointI = 0;

    labelList oldToNew(pts.size(), -1);

    // Move all of the used points to the start of the pointField and
    // truncate it

    forAll(ptUsed, ptUI)
    {
        if (ptUsed[ptUI] == true)
        {
            oldToNew[ptUI] = pointI++;
        }
    }

    inplaceReorder(oldToNew, pts);
    inplaceReorder(oldToNew, boundaryPts);

    Info<< "    Removing "
        << returnReduce(pts.size() - pointI, sumOp<label>())
        << " unused points"
        << endl;

    pts.setSize(pointI);
    boundaryPts.setSize(pointI);

    // Renumber the faces to use the new point numbers

    forAll(faces, fI)
    {
        inplaceRenumber(oldToNew, faces[fI]);
    }
}


Foam::labelList Foam::conformalVoronoiMesh::removeUnusedCells
(
    labelList& owner,
    labelList& neighbour
) const
{
    Info<< nl << "Removing unused cells" << endl;

    PackedBoolList cellUsed(number_of_vertices(), false);

    // Scan all faces to find all of the cells that are used

    forAll(owner, oI)
    {
        cellUsed[owner[oI]] = true;
    }

    forAll(neighbour, nI)
    {
        cellUsed[neighbour[nI]] = true;
    }

    label cellI = 0;

    labelList oldToNew(cellUsed.size(), -1);

    // Move all of the used cellCentres to the start of the pointField and
    // truncate it

    forAll(cellUsed, cellUI)
    {
        if (cellUsed[cellUI] == true)
        {
            oldToNew[cellUI] = cellI++;
        }
    }

    labelList newToOld(invert(cellI, oldToNew));

    // Find all of the unused cells, create a list of them, then
    // subtract one from each owner and neighbour entry for each of
    // the unused cell indices that it is above.

    DynamicList<label> unusedCells;

    forAll(cellUsed, cUI)
    {
        if (cellUsed[cUI] == false)
        {
            unusedCells.append(cUI);
        }
    }

    if (unusedCells.size() > 0)
    {
        Info<< "    Removing "
            << returnReduce(unusedCells.size(), sumOp<label>())
            <<  " unused cell labels" << endl;

        forAll(owner, oI)
        {
            label& o = owner[oI];

            o -= findLower(unusedCells, o) + 1;
        }

        forAll(neighbour, nI)
        {
            label& n = neighbour[nI];

            n -= findLower(unusedCells, n) + 1;
        }
    }

    return newToOld;
}


// ************************************************************************* //
