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

\*---------------------------------------------------------------------------*/

#include "conformalVoronoiMesh.H"
#include "initialPointsMethod.H"
#include "relaxationModel.H"
#include "faceAreaWeightModel.H"
#include "meshSearch.H"
#include "vectorTools.H"
#include "IOmanip.H"
#include "indexedCellChecks.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(conformalVoronoiMesh, 0);

}

const Foam::scalar Foam::conformalVoronoiMesh::tolParallel = 1e-3;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::cellSizeMeshOverlapsBackground() const
{
    const cellShapeControlMesh& cellSizeMesh =
        cellShapeControl_.shapeControlMesh();

    DynamicList<Foam::point> pts(number_of_vertices());

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            pts.append(topoint(vit->point()));
        }
    }

    boundBox bb(pts);

    boundBox cellSizeMeshBb = cellSizeMesh.bounds();

    bool fullyContained = true;

    if (!cellSizeMeshBb.contains(bb))
    {
        Pout<< "Triangulation not fully contained in cell size mesh."
            << endl;

        Pout<< "Cell Size Mesh Bounds = " << cellSizeMesh.bounds() << endl;
        Pout<< "cvMesh Bounds         = " << bb << endl;

        fullyContained = false;
    }

    reduce(fullyContained, andOp<unsigned int>());

    Info<< "Triangulation is "
        << (fullyContained ? "fully" : "not fully")
        << " contained in the cell size mesh"
        << endl;
}


Foam::scalar Foam::conformalVoronoiMesh::requiredSize
(
    const Foam::point& pt
) const
{
    pointIndexHit surfHit;
    label hitSurface;

    DynamicList<scalar> cellSizeHits;

    geometryToConformTo_.findSurfaceNearest
    (
        pt,
        sqr(GREAT),
        surfHit,
        hitSurface
    );

    if (!surfHit.hit())
    {
        FatalErrorIn
        (
            "Foam::tensor Foam::conformalVoronoiMesh::requiredAlignment"
        )
            << "findSurfaceNearest did not find a hit across the surfaces."
            << exit(FatalError) << endl;
    }

    cellSizeHits.append(cellShapeControls().cellSize(pt));

    // Primary alignment

    vectorField norm(1);

    allGeometry_[hitSurface].getNormal
    (
        List<pointIndexHit>(1, surfHit),
        norm
    );

    const vector np = norm[0];

    // Generate equally spaced 'spokes' in a circle normal to the
    // direction from the vertex to the closest point on the surface
    // and look for a secondary intersection.

    const vector d = surfHit.hitPoint() - pt;

    const tensor Rp = rotationTensor(vector(0,0,1), np);

    const label s = cvMeshControls().alignmentSearchSpokes();

    const scalar spanMag = geometryToConformTo_.globalBounds().mag();

    scalar totalDist = 0;

    for (label i = 0; i < s; i++)
    {
        vector spoke
        (
            Foam::cos(i*constant::mathematical::twoPi/s),
            Foam::sin(i*constant::mathematical::twoPi/s),
            0
        );

        spoke *= spanMag;

        spoke = Rp & spoke;

        pointIndexHit spokeHit;

        label spokeSurface = -1;

        // internal spoke

        geometryToConformTo_.findSurfaceNearestIntersection
        (
            pt,
            pt + spoke,
            spokeHit,
            spokeSurface
        );

        if (spokeHit.hit())
        {
            const Foam::point& hitPt = spokeHit.hitPoint();

            scalar spokeHitDistance = mag(hitPt - pt);

            cellSizeHits.append
            (
                cellShapeControls().cellSize(hitPt)
            );

            totalDist += spokeHitDistance;
        }

        //external spoke

        Foam::point mirrorPt = pt + 2*d;

        geometryToConformTo_.findSurfaceNearestIntersection
        (
            mirrorPt,
            mirrorPt + spoke,
            spokeHit,
            spokeSurface
        );

        if (spokeHit.hit())
        {
            const Foam::point& hitPt = spokeHit.hitPoint();

            scalar spokeHitDistance = mag(hitPt - mirrorPt);

            cellSizeHits.append
            (
                cellShapeControls().cellSize(hitPt)
            );

            totalDist += spokeHitDistance;
        }
    }

    scalar cellSize = 0;

    forAll(cellSizeHits, hitI)
    {
        cellSize += cellSizeHits[hitI];
    }

    return cellSize/cellSizeHits.size();
    //return cellShapeControls().cellSize(pt);
}


Foam::tensor Foam::conformalVoronoiMesh::requiredAlignment
(
    const Foam::point& pt
) const
{
    pointIndexHit surfHit;
    label hitSurface;

    geometryToConformTo_.findSurfaceNearest
    (
        pt,
        sqr(GREAT),
        surfHit,
        hitSurface
    );

    if (!surfHit.hit())
    {
        FatalErrorIn
        (
            "Foam::tensor Foam::conformalVoronoiMesh::requiredAlignment"
        )
            << "findSurfaceNearest did not find a hit across the surfaces."
            << exit(FatalError) << endl;
    }

    // Primary alignment

    vectorField norm(1);

    allGeometry_[hitSurface].getNormal
    (
        List<pointIndexHit>(1, surfHit),
        norm
    );

    const vector np = norm[0];

    // Generate equally spaced 'spokes' in a circle normal to the
    // direction from the vertex to the closest point on the surface
    // and look for a secondary intersection.

    const vector d = surfHit.hitPoint() - pt;

    const tensor Rp = rotationTensor(vector(0,0,1), np);

    const label s = cvMeshControls().alignmentSearchSpokes();

    scalar closestSpokeHitDistance = GREAT;

    pointIndexHit closestSpokeHit;

    label closestSpokeSurface = -1;

    const scalar spanMag = geometryToConformTo_.globalBounds().mag();

    for (label i = 0; i < s; i++)
    {
        vector spoke
        (
            Foam::cos(i*constant::mathematical::twoPi/s),
            Foam::sin(i*constant::mathematical::twoPi/s),
            0
        );

        spoke *= spanMag;

        spoke = Rp & spoke;

        pointIndexHit spokeHit;

        label spokeSurface = -1;

        // internal spoke

        geometryToConformTo_.findSurfaceNearestIntersection
        (
            pt,
            pt + spoke,
            spokeHit,
            spokeSurface
        );

        if (spokeHit.hit())
        {
            scalar spokeHitDistance = mag
            (
                spokeHit.hitPoint() - pt
            );

            if (spokeHitDistance < closestSpokeHitDistance)
            {
                closestSpokeHit = spokeHit;
                closestSpokeSurface = spokeSurface;
                closestSpokeHitDistance = spokeHitDistance;
            }
        }

        //external spoke

        Foam::point mirrorPt = pt + 2*d;

        geometryToConformTo_.findSurfaceNearestIntersection
        (
            mirrorPt,
            mirrorPt + spoke,
            spokeHit,
            spokeSurface
        );

        if (spokeHit.hit())
        {
            scalar spokeHitDistance = mag
            (
                spokeHit.hitPoint() - mirrorPt
            );

            if (spokeHitDistance < closestSpokeHitDistance)
            {
                closestSpokeHit = spokeHit;
                closestSpokeSurface = spokeSurface;
                closestSpokeHitDistance = spokeHitDistance;
            }
        }
    }

    if (closestSpokeSurface == -1)
    {
        WarningIn
        (
            "conformalVoronoiMesh::requiredAlignment"
            "("
                "const Foam::point& pt"
            ") const"
        )   << "No secondary surface hit found in spoke search "
            << "using " << s
            << " spokes, try increasing alignmentSearchSpokes."
            << endl;

        return I;
    }

    // Auxiliary alignment generated by spoke intersection normal.

    allGeometry_[closestSpokeSurface].getNormal
    (
        List<pointIndexHit>(1, closestSpokeHit),
        norm
    );

    const vector& na = norm[0];

    // Secondary alignment
    vector ns = np ^ na;

    if (mag(ns) < SMALL)
    {
        FatalErrorIn("conformalVoronoiMesh::requiredAlignment")
            << "Parallel normals detected in spoke search." << nl
            << "point: " << pt << nl
            << "closest surface point: " << surfHit.hitPoint() << nl
            << "closest spoke hit: " << closestSpokeHit.hitPoint() << nl
            << "np: " << surfHit.hitPoint() + np << nl
            << "ns: " << closestSpokeHit.hitPoint() + na << nl
            << exit(FatalError);
    }

    ns /= mag(ns);

    tensor Rs = rotationTensor((Rp & vector(0,1,0)), ns);

    return (Rs & Rp);
}


void Foam::conformalVoronoiMesh::insertInternalPoints
(
    List<Point>& points,
    bool distribute
)
{
    label nPoints = points.size();

    if (Pstream::parRun())
    {
        reduce(nPoints, sumOp<label>());
    }

    Info<< "    " << nPoints << " points to insert..." << endl;

    if (Pstream::parRun() && distribute)
    {
        List<Foam::point> transferPoints(points.size());

        forAll(points, pI)
        {
            transferPoints[pI] = topoint(points[pI]);
        }

        // Send the points that are not on this processor to the appropriate
        // place
        Foam::autoPtr<Foam::mapDistribute> map
        (
            decomposition_().distributePoints(transferPoints)
        );

        map().distribute(points);
    }

    label nVert = number_of_vertices();

    // using the range insert (faster than inserting points one by one)
    insert(points.begin(), points.end());

    label nInserted(number_of_vertices() - nVert);

    if (Pstream::parRun())
    {
        reduce(nInserted, sumOp<label>());
    }

    Info<< "    " << nInserted << " points inserted"
        << ", failed to insert " << nPoints - nInserted
        << " ("
        << 100.0*(nPoints - nInserted)/nInserted
        << " %)"<< endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->uninitialised())
        {
            vit->index() = getNewVertexIndex();
            vit->type() = Vb::vtInternal;
        }
    }
}


void Foam::conformalVoronoiMesh::insertPoints
(
    List<Vb>& vertices,
    bool distribute
)
{
    if (Pstream::parRun() && distribute)
    {
        const label preDistributionSize = vertices.size();

        List<Foam::point> pts(preDistributionSize);

        forAll(vertices, vI)
        {
            const Foam::point& pt = topoint(vertices[vI].point());

            pts[vI] = pt;
        }

        // Distribute points to their appropriate processor
        autoPtr<mapDistribute> map
        (
            decomposition_().distributePoints(pts)
        );

        map().distribute(vertices);

        forAll(vertices, vI)
        {
            vertices[vI].procIndex() = Pstream::myProcNo();
        }
    }

    rangeInsertWithInfo
    (
        vertices.begin(),
        vertices.end(),
        true
    );
}


void Foam::conformalVoronoiMesh::insertSurfacePointPairs
(
    const pointIndexHitAndFeatureList& surfaceHits,
    const fileName fName
)
{
    DynamicList<Vb> pts(2.0*surfaceHits.size());

    forAll(surfaceHits, i)
    {
        vectorField norm(1);

        const pointIndexHit surfaceHit = surfaceHits[i].first();
        const label featureIndex = surfaceHits[i].second();

        allGeometry_[featureIndex].getNormal
        (
            List<pointIndexHit>(1, surfaceHit),
            norm
        );

        const vector& normal = norm[0];

        const Foam::point& surfacePt(surfaceHit.hitPoint());

        if (geometryToConformTo_.isBaffle(featureIndex))
        {
            createBafflePointPair
            (
                pointPairDistance(surfacePt),
                surfacePt,
                normal,
                pts
            );
        }
        else
        {
            createPointPair
            (
                pointPairDistance(surfacePt),
                surfacePt,
                normal,
                pts
            );
        }
    }

    insertPoints(pts, true);

    if (cvMeshControls().objOutput() && fName != fileName::null)
    {
        writePoints(fName, pts);
    }
}


void Foam::conformalVoronoiMesh::insertEdgePointGroups
(
    const pointIndexHitAndFeatureList& edgeHits,
    const fileName fName
)
{
    DynamicList<Vb> pts(3.0*edgeHits.size());

    forAll(edgeHits, i)
    {
        const extendedFeatureEdgeMesh& feMesh
        (
            geometryToConformTo_.features()[edgeHits[i].second()]
        );

        createEdgePointGroup(feMesh, edgeHits[i].first(), pts);
    }

    pts.shrink();

    insertPoints(pts, true);

    if (cvMeshControls().objOutput() && fName != fileName::null)
    {
        writePoints(fName, pts);
    }
}


bool Foam::conformalVoronoiMesh::nearFeaturePt(const Foam::point& pt) const
{
    scalar exclusionRangeSqr = featurePointExclusionDistanceSqr(pt);

    pointIndexHit info;
    label featureHit;

    geometryToConformTo_.findFeaturePointNearest
    (
        pt,
        exclusionRangeSqr,
        info,
        featureHit
    );

    return info.hit();
}


void Foam::conformalVoronoiMesh::insertInitialPoints()
{
    Info<< nl << "Inserting initial points" << endl;

    timeCheck("Before initial points call");

    List<Point> initPts = initialPointsMethod_->initialPoints();

    timeCheck("After initial points call");

    // Assume that the initial points method made the correct decision for
    // which processor each point should be on, so give distribute = false
    insertInternalPoints(initPts);
}


Foam::scalar Foam::conformalVoronoiMesh::calculateLoadUnbalance() const
{
    label nRealVertices = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        // Only store real vertices that are not feature vertices
        if (vit->real() && !vit->featurePoint())
        {
            nRealVertices++;
        }
    }

    scalar globalNRealVertices = returnReduce
    (
        nRealVertices,
        sumOp<label>()
    );

    scalar unbalance = returnReduce
    (
        mag(1.0 - nRealVertices/(globalNRealVertices/Pstream::nProcs())),
        maxOp<scalar>()
    );

    Info<< "    Processor unbalance " << unbalance << endl;

    return unbalance;
}


bool Foam::conformalVoronoiMesh::distributeBackground()
{
    if (!Pstream::parRun())
    {
        return false;
    }

    Info<< nl << "Redistributing points" << endl;

    timeCheck("Before distribute");

    label iteration = 0;

    scalar previousLoadUnbalance = 0;

    while (true)
    {
        scalar maxLoadUnbalance = calculateLoadUnbalance();

        if
        (
            maxLoadUnbalance <= cvMeshControls().maxLoadUnbalance()
         || maxLoadUnbalance <= previousLoadUnbalance
        )
        {
            // If this is the first iteration, return false, if it was a
            // subsequent one, return true;
            return iteration != 0;
        }

        previousLoadUnbalance = maxLoadUnbalance;

        Info<< "    Total number of vertices before redistribution "
            << returnReduce(label(number_of_vertices()), sumOp<label>())
            << endl;

        const fvMesh& bMesh = decomposition_().mesh();

        volScalarField cellWeights
        (
            IOobject
            (
                "cellWeights",
                bMesh.time().timeName(),
                bMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bMesh,
            dimensionedScalar("weight", dimless, 1e-2),
            zeroGradientFvPatchScalarField::typeName
        );

        meshSearch cellSearch(bMesh, polyMesh::FACEPLANES);

        labelList cellVertices(bMesh.nCells(), 0);

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            // Only store real vertices that are not feature vertices
            if (vit->real() && !vit->featurePoint())
            {
                pointFromPoint v = topoint(vit->point());

                label cellI = cellSearch.findCell(v);

                if (cellI == -1)
                {
//                     Pout<< "findCell conformalVoronoiMesh::distribute "
//                         << "findCell "
//                         << vit->type() << " "
//                         << vit->index() << " "
//                         << v << " "
//                         << cellI
//                         << " find nearest cellI ";

                    cellI = cellSearch.findNearestCell(v);
                }

                cellVertices[cellI]++;
            }
        }

        forAll(cellVertices, cI)
        {
            // Give a small but finite weight for empty cells.  Some
            // decomposition methods have difficulty with integer overflows in
            // the sum of the normalised weight field.
            cellWeights.internalField()[cI] = max
            (
                cellVertices[cI],
                1e-2
            );
        }

        autoPtr<mapDistributePolyMesh> mapDist = decomposition_().distribute
        (
            cellWeights
        );

        cellShapeControl_.shapeControlMesh().distribute(decomposition_);

        distribute();

        timeCheck("After distribute");

        iteration++;
    }

    return true;
}


void Foam::conformalVoronoiMesh::distribute()
{
    if (!Pstream::parRun())
    {
        return ;
    }

    autoPtr<mapDistribute> mapDist =
        DistributedDelaunayMesh<Delaunay>::distribute(decomposition_());

    DynamicList<Foam::point> points(number_of_vertices());
    DynamicList<Foam::indexedVertexEnum::vertexType> types
    (
        number_of_vertices()
    );
    DynamicList<scalar> sizes(number_of_vertices());
    DynamicList<tensor> alignments(number_of_vertices());

    for
    (
        Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            points.append(topoint(vit->point()));
            types.append(vit->type());
            sizes.append(vit->targetCellSize());
            alignments.append(vit->alignment());
        }
    }

    mapDist().distribute(points);
    mapDist().distribute(types);
    mapDist().distribute(sizes);
    mapDist().distribute(alignments);

    // Reset the entire tessellation
    DelaunayMesh<Delaunay>::reset();

    Info<< nl << "    Inserting distributed tessellation" << endl;

    // Internal points have to be inserted first

    DynamicList<Vb> verticesToInsert(points.size());

    forAll(points, pI)
    {
        verticesToInsert.append
        (
            Vb
            (
                toPoint<Point>(points[pI]),
                -1,
                types[pI],
                Pstream::myProcNo()
            )
        );

        verticesToInsert.last().targetCellSize() = sizes[pI];
        verticesToInsert.last().alignment() = alignments[pI];
    }

    this->rangeInsertWithInfo
    (
        verticesToInsert.begin(),
        verticesToInsert.end(),
        true
    );

    Info<< "    Total number of vertices after redistribution "
        << returnReduce
           (
               label(number_of_vertices()), sumOp<label>()
           )
        << endl;
}


void Foam::conformalVoronoiMesh::buildCellSizeAndAlignmentMesh()
{
    cellShapeControl_.initialMeshPopulation(decomposition_);

    cellShapeControlMesh& cellSizeMesh = cellShapeControl_.shapeControlMesh();

    if (Pstream::parRun())
    {
        cellSizeMesh.distribute(decomposition_);
    }

    label nMaxIter = 2;

    for (label i = 0; i < nMaxIter; ++i)
    {
        label nAdded = cellShapeControl_.refineMesh(decomposition_);
        reduce(nAdded, sumOp<label>());

        if (Pstream::parRun())
        {
            cellSizeMesh.distribute(decomposition_);
        }

        if (nAdded == 0)
        {
            break;
        }

        Info<< "    Iteration " << i << ": Added = " << nAdded << " points"
            << endl;
    }

    cellShapeControl_.smoothMesh();

    Info<< "Background cell size and alignment mesh:" << endl;
    cellSizeMesh.printInfo(Info);

//    cellSizeMesh.write();
}


void Foam::conformalVoronoiMesh::storeSizesAndAlignments()
{
    DynamicList<Point> storePts(number_of_vertices());

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalPoint())
        {
            storePts.append(vit->point());
        }
    }

    storePts.shrink();

    storeSizesAndAlignments(storePts);
}


void Foam::conformalVoronoiMesh::storeSizesAndAlignments
(
    const List<Point>& storePts
)
{
//    timeCheck("Start of storeSizesAndAlignments");
//
//    Info << nl << "Store size and alignment" << endl;
//
//    sizeAndAlignmentLocations_.setSize(storePts.size());
//
//    storedSizes_.setSize(sizeAndAlignmentLocations_.size());
//
//    storedAlignments_.setSize(sizeAndAlignmentLocations_.size());
//
//    label i = 0;
//
//    //checkCellSizing();
//
//    for
//    (
//        List<Point>::const_iterator pit = storePts.begin();
//        pit != storePts.end();
//        ++pit
//    )
//    {
//        pointFromPoint pt = topoint(*pit);
//
////        storedAlignments_[i] = requiredAlignment(pt);
////
////        storedSizes_[i] = cellShapeControls().cellSize(pt);
//
//        cellShapeControls().cellSizeAndAlignment
//        (
//            pt,
//            storedSizes_[i],
//            storedAlignments_[i]
//        );
//
//        i++;
//    }
//
//    timeCheck("Sizes and alignments calculated, build tree");
//
//    buildSizeAndAlignmentTree();
//
//    timeCheck("Size and alignment tree built");
}


void Foam::conformalVoronoiMesh::updateSizesAndAlignments
(
    const List<Point>& storePts
)
{
    // This function is only used in serial, the background redistribution
    // triggers this when unbalance is detected in parallel.

    if
    (
        !Pstream::parRun()
     && runTime_.run()
     && runTime_.timeIndex()
      % cvMeshControls().sizeAndAlignmentRebuildFrequency() == 0
    )
    {
        storeSizesAndAlignments(storePts);

        timeCheck("Updated sizes and alignments");
    }
}


const Foam::indexedOctree<Foam::treeDataPoint>&
Foam::conformalVoronoiMesh::sizeAndAlignmentTree() const
{
    if (sizeAndAlignmentTreePtr_.empty())
    {
        buildSizeAndAlignmentTree();
    }

    return sizeAndAlignmentTreePtr_();
}


void Foam::conformalVoronoiMesh::setVertexSizeAndAlignment()
{
//    Info<< nl << "Looking up target cell alignment and size" << endl;
//
//    const indexedOctree<treeDataPoint>& tree = sizeAndAlignmentTree();
//
//    for
//    (
//        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
//        vit != finite_vertices_end();
//        vit++
//    )
//    {
//        if
//        (
//            vit->internalOrBoundaryPoint()
//         || vit->referredInternalOrBoundaryPoint()
//        )
//        {
//            pointFromPoint pt = topoint(vit->point());
//
//            pointIndexHit info = tree.findNearest(pt, sqr(GREAT));
//
//            if (info.hit())
//            {
//                vit->alignment() = storedAlignments_[info.index()];
//
//                vit->targetCellSize() = storedSizes_[info.index()];
//            }
//            else
//            {
//                WarningIn
//                (
//                    "void "
//                    "Foam::conformalVoronoiMesh::setVertexSizeAndAlignment()"
//                )
//                    << "Point " << pt << " did not find a nearest point "
//                    << " for alignment and size lookup." << endl;
//
//                vit->alignment() = cellShapeControls().cellAlignment(pt);
//
//                vit->targetCellSize() = cellShapeControls().cellSize(pt);
//            }
//        }
//    }

    Info<< nl << "Calculating target cell alignment and size" << endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            pointFromPoint pt = topoint(vit->point());

            cellShapeControls().cellSizeAndAlignment
            (
                pt,
                vit->targetCellSize(),
                vit->alignment()
            );

            //vit->alignment() = tensor(1,0,0,0,1,0,0,0,1);
            //vit->alignment() = requiredAlignment(pt);

            //vit->targetCellSize() = cellShapeControls().cellSize(pt);
        }
    }
}


Foam::face Foam::conformalVoronoiMesh::buildDualFace
(
    const Delaunay::Finite_edges_iterator& eit
) const
{
    Cell_circulator ccStart = incident_cells(*eit);
    Cell_circulator cc1 = ccStart;
    Cell_circulator cc2 = cc1;

    // Advance the second circulator so that it always stays on the next
    // cell around the edge;
    cc2++;

    DynamicList<label> verticesOnFace;

    label nUniqueVertices = 0;

    do
    {
        if
        (
            cc1->hasFarPoint() || cc2->hasFarPoint()
         || is_infinite(cc1) || is_infinite(cc2)
        )
        {
            Cell_handle c = eit->first;
            Vertex_handle vA = c->vertex(eit->second);
            Vertex_handle vB = c->vertex(eit->third);

            drawDelaunayCell(Pout, cc1);
            drawDelaunayCell(Pout, cc2);

            FatalErrorIn("Foam::conformalVoronoiMesh::buildDualFace")
                << "Dual face uses circumcenter defined by a "
                << "Delaunay tetrahedron with no internal "
                << "or boundary points.  Defining Delaunay edge ends: "
                << topoint(vA->point()) << " "
                << topoint(vB->point()) << nl
                << exit(FatalError);
        }

        label cc1I = cc1->cellIndex();

        label cc2I = cc2->cellIndex();

        if (cc1I != cc2I)
        {
            if (findIndex(verticesOnFace, cc1I) == -1)
            {
                nUniqueVertices++;
            }

            verticesOnFace.append(cc1I);
        }

        cc1++;

        cc2++;

    } while (cc1 != ccStart);

    verticesOnFace.shrink();

    if (verticesOnFace.size() >= 3 && nUniqueVertices < 3)
    {
        // There are not enough unique vertices on this face to
        // justify its size, it may have a form like:

        // Vertices:
        // A                                  B
        // A                                  B

        // Face:
        // ABAB

        // Setting the size to be below 3, so that it will not be
        // created

        verticesOnFace.setSize(nUniqueVertices);
    }

    return face(verticesOnFace);
}


Foam::label Foam::conformalVoronoiMesh::maxFilterCount
(
    const Delaunay::Finite_edges_iterator& eit
) const
{
    Cell_circulator ccStart = incident_cells(*eit);
    Cell_circulator cc = ccStart;

    label maxFC = 0;

    do
    {
        if (cc->hasFarPoint())
        {
            Cell_handle c = eit->first;
            Vertex_handle vA = c->vertex(eit->second);
            Vertex_handle vB = c->vertex(eit->third);

            FatalErrorIn("Foam::conformalVoronoiMesh::buildDualFace")
                << "Dual face uses circumcenter defined by a "
                << "Delaunay tetrahedron with no internal "
                << "or boundary points.  Defining Delaunay edge ends: "
                << topoint(vA->point()) << " "
                << topoint(vB->point()) << nl
                << exit(FatalError);
        }

        if (cc->filterCount() > maxFC)
        {
            maxFC = cc->filterCount();
        }

        cc++;

    } while (cc != ccStart);

    return maxFC;
}


bool Foam::conformalVoronoiMesh::ownerAndNeighbour
(
    Vertex_handle vA,
    Vertex_handle vB,
    label& owner,
    label& neighbour
) const
{
    bool reverse = false;

    owner = -1;
    neighbour = -1;

    label dualCellIndexA = vA->index();

    if (!vA->internalOrBoundaryPoint() || vA->referred())
    {
        dualCellIndexA = -1;
    }

    label dualCellIndexB = vB->index();

    if (!vB->internalOrBoundaryPoint() || vB->referred())
    {
        dualCellIndexB = -1;
    }

    if (dualCellIndexA == -1 && dualCellIndexB == -1)
    {
        FatalErrorIn
        (
            "bool Foam::conformalVoronoiMesh::ownerAndNeighbour"
            "("
                "Vertex_handle vA,"
                "Vertex_handle vB,"
                "label& owner,"
                "label& neighbour"
            ") const"
        )
            << "Attempting to create a face joining "
            << "two unindexed dual cells "
            << exit(FatalError);
    }
    else if (dualCellIndexA == -1 || dualCellIndexB == -1)
    {
        // boundary face, find which is the owner

        if (dualCellIndexA == -1)
        {
            owner = dualCellIndexB;

            reverse = true;
        }
        else
        {
            owner = dualCellIndexA;
        }
    }
    else
    {
        // internal face, find the lower cell to be the owner

        if (dualCellIndexB > dualCellIndexA)
        {
            owner = dualCellIndexA;
            neighbour = dualCellIndexB;
        }
        else
        {
            owner = dualCellIndexB;
            neighbour = dualCellIndexA;

            // reverse face order to correctly orientate normal
            reverse = true;
        }
    }

    return reverse;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::conformalVoronoiMesh
(
    const Time& runTime,
    const dictionary& cvMeshDict
)
:
    DistributedDelaunayMesh<Delaunay>(),
    runTime_(runTime),
    rndGen_(64293*Pstream::myProcNo()),
    cvMeshControls_(cvMeshDict),
    allGeometry_
    (
        IOobject
        (
            "cvSearchableSurfaces",
            runTime_.constant(),
            "triSurface",
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        cvMeshDict.subDict("geometry")
    ),
    geometryToConformTo_
    (
        runTime_,
        rndGen_,
        allGeometry_,
        cvMeshDict.subDict("surfaceConformation")
    ),
    cellShapeControl_
    (
        runTime_,
        cvMeshDict.subDict("motionControl"),
        allGeometry_,
        geometryToConformTo_
    ),
    limitBounds_(),
    featureVertices_(),
    featurePointLocations_(),
    edgeLocationTreePtr_(),
    surfacePtLocationTreePtr_(),
    sizeAndAlignmentLocations_(),
    storedSizes_(),
    storedAlignments_(),
    sizeAndAlignmentTreePtr_(),
    surfaceConformationVertices_(),
    initialPointsMethod_
    (
        initialPointsMethod::New
        (
            cvMeshDict.subDict("initialPoints"),
            *this
        )
    ),
    relaxationModel_
    (
        relaxationModel::New
        (
            cvMeshDict.subDict("motionControl"),
            runTime_
        )
    ),
    faceAreaWeightModel_
    (
        faceAreaWeightModel::New
        (
            cvMeshDict.subDict("motionControl")
        )
    ),
    decomposition_()
{
    if (cvMeshControls().objOutput())
    {
        geometryToConformTo_.writeFeatureObj("cvMesh");
    }

    if (Pstream::parRun())
    {
        decomposition_.reset
        (
            new backgroundMeshDecomposition
            (
                runTime_,
                rndGen_,
                geometryToConformTo_,
                cvMeshDict.subDict("backgroundMeshDecomposition")
            )
        );
    }

    buildCellSizeAndAlignmentMesh();

    insertInitialPoints();

    insertFeaturePoints();

    setVertexSizeAndAlignment();

    cellSizeMeshOverlapsBackground();

    // Improve the guess that the backgroundMeshDecomposition makes with the
    // initial positions.  Use before building the surface conformation to
    // better balance the surface conformation load.
    distributeBackground();

    buildSurfaceConformation();

    // The introduction of the surface conformation may have distorted the
    // balance of vertices, distribute if necessary.
    distributeBackground();

    if (Pstream::parRun())
    {
        sync(decomposition_().procBounds());
    }

    // Do not store the surface conformation until after it has been
    // (potentially) redistributed.
    storeSurfaceConformation();

    // Use storeSizesAndAlignments with no feed points because all background
    // points may have been distributed.
    storeSizesAndAlignments();

    // Report any Delaunay vertices that do not think that they are in the
    // domain the processor they are on.
    // reportProcessorOccupancy();

    cellSizeMeshOverlapsBackground();

    printVertexInfo();

    if (cvMeshControls().objOutput())
    {
        writePoints
        (
            "internalPoints_" + runTime_.timeName() + ".obj",
            Foam::indexedVertexEnum::vtUnassigned,
            Foam::indexedVertexEnum::vtExternalFeaturePoint
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::~conformalVoronoiMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::move()
{
    timeCheck("Start of move");

    scalar relaxation = relaxationModel_->relaxation();

    Info<< nl << "Relaxation = " << relaxation << endl;

    pointField dualVertices(number_of_finite_cells());

    this->resetCellCount();

    // Find the dual point of each tetrahedron and assign it an index.
    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        cit->cellIndex() = Cb::ctUnassigned;

        if (cit->anyInternalOrBoundaryDualVertex())
        {
            cit->cellIndex() = getNewCellIndex();

            dualVertices[cit->cellIndex()] = cit->dual();
        }

        if (cit->hasFarPoint())
        {
            cit->cellIndex() = Cb::ctFar;
        }
    }

    dualVertices.setSize(cellCount());

    setVertexSizeAndAlignment();

    timeCheck("Determined sizes and alignments");

    Info<< nl << "Determining vertex displacements" << endl;

    vectorField cartesianDirections(3);

    cartesianDirections[0] = vector(1, 0, 0);
    cartesianDirections[1] = vector(0, 1, 0);
    cartesianDirections[2] = vector(0, 0, 1);

    vectorField displacementAccumulator
    (
        number_of_vertices(),
        vector::zero
    );

    PackedBoolList pointToBeRetained
    (
        number_of_vertices(),
        true
    );

    DynamicList<Point> pointsToInsert(number_of_vertices());

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
            (
                vA->internalPoint() && !vA->referred()
             && vB->internalOrBoundaryPoint()
            )
         || (
                vB->internalPoint() && !vB->referred()
             && vA->internalOrBoundaryPoint()
            )
        )
        {
            pointFromPoint dVA = topoint(vA->point());
            pointFromPoint dVB = topoint(vB->point());

            Field<vector> alignmentDirsA
            (
                vA->alignment().T() & cartesianDirections
            );
            Field<vector> alignmentDirsB
            (
                vB->alignment().T() & cartesianDirections
            );

            Field<vector> alignmentDirs(3);

            forAll(alignmentDirsA, aA)
            {
                const vector& a = alignmentDirsA[aA];

                scalar maxDotProduct = 0.0;

                forAll(alignmentDirsB, aB)
                {
                    const vector& b = alignmentDirsB[aB];

                    const scalar dotProduct = a & b;

                    if (mag(dotProduct) > maxDotProduct)
                    {
                        maxDotProduct = mag(dotProduct);

                        alignmentDirs[aA] = a + sign(dotProduct)*b;

                        alignmentDirs[aA] /= mag(alignmentDirs[aA]);
                    }
                }
            }

            vector rAB = dVA - dVB;

            scalar rABMag = mag(rAB);

            if (rABMag < SMALL)
            {
                // Removal of close points

                if
                (
                    vA->internalPoint() && !vA->referred()
                 && vB->internalPoint() && !vB->referred()
                )
                {
                    // Only insert a point at the midpoint of
                    // the short edge if neither attached
                    // point has already been identified to be
                    // removed.

                    if
                    (
                        pointToBeRetained[vA->index()] == true
                     && pointToBeRetained[vB->index()] == true
                    )
                    {
                        pointsToInsert.append
                        (
                            toPoint<Point>(0.5*(dVA + dVB))
                        );
                    }
                }

                if (vA->internalPoint() && !vA->referred())
                {
                    pointToBeRetained[vA->index()] = false;
                }

                if (vB->internalPoint() && !vB->referred())
                {
                    pointToBeRetained[vB->index()] = false;
                }

                // Do not consider this Delaunay edge any further

                continue;
            }

            forAll(alignmentDirs, aD)
            {
                vector& alignmentDir = alignmentDirs[aD];

                scalar dotProd = rAB & alignmentDir;

                if (dotProd < 0)
                {
                    // swap the direction of the alignment so that has the
                    // same sense as rAB
                    alignmentDir *= -1;
                    dotProd *= -1;
                }

                const scalar alignmentDotProd = dotProd/rABMag;

                if
                (
                    alignmentDotProd
                  > cvMeshControls().cosAlignmentAcceptanceAngle()
                )
                {
                    scalar targetCellSize = averageCellSize(vA, vB);

                    scalar targetFaceArea = sqr(targetCellSize);

                    const vector originalAlignmentDir = alignmentDir;

                    // Update cell size and face area
                    cellShapeControls().aspectRatio().updateCellSizeAndFaceArea
                    (
                        alignmentDir,
                        targetFaceArea,
                        targetCellSize
                    );

                    // Vector to move end points around middle of vector
                    // to align edge (i.e. dual face normal) with alignment
                    // directions.
                    vector delta = alignmentDir - 0.5*rAB;

                    face dualFace = buildDualFace(eit);

//                    Pout<< dualFace << endl;
//                    Pout<< "    " << vA->info() << endl;
//                    Pout<< "    " << vB->info() << endl;

                    const scalar faceArea = dualFace.mag(dualVertices);

                    // Update delta vector
                    cellShapeControls().aspectRatio().updateDeltaVector
                    (
                        originalAlignmentDir,
                        targetCellSize,
                        rABMag,
                        delta
                    );

                    if (targetFaceArea == 0)
                    {
                        Pout<< vA->info() << vB->info();

                        Cell_handle ch = locate(vA->point());
                        if (is_infinite(ch))
                        {
                            Pout<< "vA " << vA->targetCellSize() << endl;
                        }

                        ch = locate(vB->point());
                        if (is_infinite(ch))
                        {
                            Pout<< "vB " << vB->targetCellSize() << endl;
                        }
                    }

                    delta *= faceAreaWeightModel_->faceAreaWeight
                    (
                        faceArea/targetFaceArea
                    );

                    if
                    (
                        (
                            (vA->internalPoint() && vB->internalPoint())
                         && (!vA->referred() || !vB->referred())
//                         ||
//                            (
//                                vA->referredInternalPoint()
//                             && vB->referredInternalPoint()
//                            )
                        )
                     && rABMag
                      > cvMeshControls().insertionDistCoeff()*targetCellSize
                     && faceArea
                      > cvMeshControls().faceAreaRatioCoeff()*targetFaceArea
                     && alignmentDotProd
                      > cvMeshControls().cosInsertionAcceptanceAngle()
                    )
                    {
                        // Point insertion
                        if
                        (
                            !geometryToConformTo_.findSurfaceAnyIntersection
                            (
                                dVA,
                                dVB
                            )
                        )
                        {
                            const Foam::point& newPt = 0.5*(dVA + dVB);

                            if (positionOnThisProc(newPt))
                            {
                                // Prevent insertions spanning surfaces
                                pointsToInsert.append(toPoint<Point>(newPt));
                            }
                        }
                    }
                    else if
                    (
                        (
                            (vA->internalPoint() && !vA->referred())
                         || (vB->internalPoint() && !vB->referred())
                        )
                     && rABMag
                      < cvMeshControls().removalDistCoeff()
                       *targetCellSize
                    )
                    {
                        // Point removal
                        if
                        (
                            vA->internalPoint() && !vA->referred()
                         && vB->internalPoint() && !vB->referred()
                        )
                        {
                            // Only insert a point at the midpoint of
                            // the short edge if neither attached
                            // point has already been identified to be
                            // removed.
                            if
                            (
                                pointToBeRetained[vA->index()] == true
                             && pointToBeRetained[vB->index()] == true
                            )
                            {
                                pointsToInsert.append
                                (
                                    toPoint<Point>(0.5*(dVA + dVB))
                                );
                            }
                        }

                        if (vA->internalPoint() && !vA->referred())
                        {
                            pointToBeRetained[vA->index()] = false;
                        }

                        if (vB->internalPoint() && !vB->referred())
                        {
                            pointToBeRetained[vB->index()] = false;
                        }
                    }
                    else
                    {
                        if (vA->internalPoint() && !vA->referred())
                        {
                            displacementAccumulator[vA->index()] += delta;
                        }

                        if (vB->internalPoint() && !vB->referred())
                        {
                            displacementAccumulator[vB->index()] -= delta;
                        }
                    }
                }
            }
        }
    }

    Info<< "Limit displacements" << endl;

    // Limit displacements that pierce, or get too close to the surface
    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() && !vit->referred())
        {
            if (pointToBeRetained[vit->index()] == true)
            {
                limitDisplacement
                (
                    vit,
                    displacementAccumulator[vit->index()]
                );
            }
        }
    }

    vector totalDisp = gSum(displacementAccumulator);
    scalar totalDist = gSum(mag(displacementAccumulator));

    displacementAccumulator *= relaxation;

    Info<< "Sum displacements" << endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() && !vit->referred())
        {
            if (pointToBeRetained[vit->index()] == true)
            {
                // Convert vit->point() to FOAM vector (double) to do addition,
                // avoids memory increase because a record of the constructions
                // would be kept otherwise.
                // See cgal-discuss@lists-sop.inria.fr:
                // "Memory issue with openSUSE 11.3, exact kernel, adding
                //  points/vectors"
                // 14/1/2011.
                // Only necessary if using an exact constructions kernel
                // (extended precision)

                pointsToInsert.append
                (
                    toPoint<Point>
                    (
                        topoint(vit->point())
                      + displacementAccumulator[vit->index()]
                    )
                );
            }
        }
    }

    pointsToInsert.shrink();

    // Save displacements to file.
    if (cvMeshControls().objOutput() && runTime_.outputTime())
    {
        Pout<< "Writing point displacement vectors to file." << endl;
        OFstream str("displacements_" + runTime_.timeName() + ".obj");

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            if (vit->internalPoint() && !vit->referred())
            {
                if (pointToBeRetained[vit->index()] == true)
                {
                    meshTools::writeOBJ(str, topoint(vit->point()));

                    str << "vn "
                        << displacementAccumulator[vit->index()][0] << " "
                        << displacementAccumulator[vit->index()][1] << " "
                        << displacementAccumulator[vit->index()][2] << " "
                        << endl;
                }
            }
        }
    }

    // Remove the entire tessellation
    DelaunayMesh<Delaunay>::reset();

    timeCheck("Displacement calculated");

    Info<< nl << "Inserting displaced tessellation" << endl;

    insertInternalPoints(pointsToInsert, true);

    reinsertFeaturePoints(true);

    // Remove internal points that have been inserted outside the surface.
//    label internalPtIsOutside = 0;
//
//    for
//    (
//        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
//        vit != finite_vertices_end();
//        ++vit
//    )
//    {
//        if (vit->internalPoint() && !vit->referred())
//        {
//            bool inside = geometryToConformTo_.inside
//            (
//                topoint(vit->point())
//            );
//
//            if (!inside)
//            {
//                remove(vit);
//                internalPtIsOutside++;
//            }
//        }
//    }
//
//    Info<< "    " << internalPtIsOutside
//        << " internal points were inserted outside the domain. "
//        << "They have been removed." << endl;

    // Fix points that have not been significantly displaced
//    for
//    (
//        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
//        vit != finite_vertices_end();
//        ++vit
//    )
//    {
//        if (vit->internalPoint())
//        {
//            if
//            (
//                mag(displacementAccumulator[vit->index()])
//              < 0.1*targetCellSize(topoint(vit->point()))
//            )
//            {
//                vit->setVertexFixed();
//            }
//        }
//    }

    timeCheck("Internal points inserted");

    {
        // Check that no index is shared between any of the local points
        labelHashSet usedIndices;
        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            if (!vit->referred() && !usedIndices.insert(vit->index()))
            {
                FatalErrorIn("Foam::conformalVoronoiMesh::move()")
                    << "Index already used! Could not insert: " << nl
                    << vit->info()
                    << abort(FatalError);
            }
        }
    }

    conformToSurface();

    if (cvMeshControls().objOutput())
    {
        writePoints
        (
            "internalPoints_" + runTime_.timeName() + ".obj",
            Foam::indexedVertexEnum::vtInternal
        );
    }

    if (cvMeshControls().objOutput() && runTime_.outputTime())
    {
        writeBoundaryPoints("boundaryPoints_" + runTime_.timeName() + ".obj");
    }

    timeCheck("After conformToSurface");

    printVertexInfo();

    // Write the intermediate mesh, do not filter the dual faces.
    if (runTime_.outputTime())
    {
        writeMesh(runTime_.timeName());
    }

    updateSizesAndAlignments(pointsToInsert);

    Info<< nl
        << "Total displacement = " << totalDisp << nl
        << "Total distance = " << totalDist << nl
        << endl;
}


bool Foam::conformalVoronoiMesh::positionOnThisProc
(
    const Foam::point& pt
) const
{
    if (Pstream::parRun())
    {
        return decomposition_().positionOnThisProcessor(pt);
    }

    return true;
}


Foam::boolList Foam::conformalVoronoiMesh::positionOnThisProc
(
    const Foam::List<Foam::point>& pts
) const
{
    if (Pstream::parRun())
    {
        return decomposition_().positionOnThisProcessor(pts);
    }

    return boolList(pts.size(), true);
}


Foam::labelList Foam::conformalVoronoiMesh::positionProc
(
    const Foam::List<Foam::point>& pts
) const
{
    if (!Pstream::parRun())
    {
        return labelList(pts.size(), -1);
    }

    return decomposition_().processorPosition(pts);
}


Foam::List<Foam::List<Foam::pointIndexHit> >
Foam::conformalVoronoiMesh::intersectsProc
(
    const List<Foam::point>& starts,
    const List<Foam::point>& ends
) const
{
    if (!Pstream::parRun())
    {
        return List<List<pointIndexHit> >(starts.size());
    }

    return decomposition_().intersectsProcessors(starts, ends, false);
}


//Foam::labelListList Foam::conformalVoronoiMesh::overlapsProc
//(
//    const List<Foam::point>& centres,
//    const List<scalar>& radiusSqrs
//) const
//{
//    if (!Pstream::parRun())
//    {
//        return labelListList(centres.size(), labelList(0));
//    }
//
////    DynamicList<Foam::point> pts(number_of_vertices());
//
////    for
////    (
////        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
////        vit != finite_vertices_end();
////        vit++
////    )
////    {
////        pts.append(topoint(vit->point()));
////    }
////
////    dynamicIndexedOctree<dynamicTreeDataPoint> vertexOctree
////    (
////        dynamicTreeDataPoint(pts),
////        treeBoundBox(min(pts), max(pts)),
////        10, // maxLevel
////        10, // leafSize
////        3.0 // duplicity
////    );
//
//    return decomposition_().overlapsProcessors
//    (
//        centres,
//        radiusSqrs,
//        *this,
//        false//,
////        vertexOctree
//    );
//}


void Foam::conformalVoronoiMesh::checkCoPlanarCells() const
{
    typedef CGAL::Exact_predicates_exact_constructions_kernel   Kexact;
    typedef CGAL::Point_3<Kexact>                               PointExact;

    if (!is_valid())
    {
        Pout<< "Triangulation is invalid!" << endl;
    }

    OFstream str("badCells.obj");

    label badCells = 0;

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        const scalar quality = cvMeshChecks::coplanarTet(cit, 1e-16);

        if (quality == 0)
        {
            Pout<< "COPLANAR: " << cit->info() << nl
                << "    quality = " << quality << nl
                << "    dual    = " << topoint(cit->dual()) << endl;

            drawDelaunayCell(str, cit, badCells++);

            FixedList<PointExact, 4> cellVerticesExact(PointExact(0,0,0));
            forAll(cellVerticesExact, vI)
            {
                cellVerticesExact[vI] = PointExact
                (
                    cit->vertex(vI)->point().x(),
                    cit->vertex(vI)->point().y(),
                    cit->vertex(vI)->point().z()
                );
            }

            PointExact synchronisedDual = CGAL::circumcenter<Kexact>
            (
                cellVerticesExact[0],
                cellVerticesExact[1],
                cellVerticesExact[2],
                cellVerticesExact[3]
            );

            Foam::point exactPt
            (
                CGAL::to_double(synchronisedDual.x()),
                CGAL::to_double(synchronisedDual.y()),
                CGAL::to_double(synchronisedDual.z())
            );

            Info<< "inexact = " << cit->dual() << nl
                << "exact   = " << exactPt << endl;
        }
    }

    Pout<< "There are " << badCells << " bad cells out of "
        << number_of_finite_cells() << endl;


    label nNonGabriel = 0;
    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        if (!is_Gabriel(*fit))
        {
            nNonGabriel++;//Pout<< "Non-gabriel face" << endl;
        }
    }

    Pout<< "There are " << nNonGabriel << " non-Gabriel faces out of "
        << number_of_finite_facets() << endl;
}


// ************************************************************************* //
