/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "conformalVoronoiMesh.H"
#include "initialPointsMethod.H"
#include "relaxationModel.H"
#include "faceAreaWeightModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::conformalVoronoiMesh
(
    const Time& runTime,
    const IOdictionary& cvMeshDict
)
:
    HTriangulation(),
    runTime_(runTime),
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
        *this,
        allGeometry_,
        cvMeshDict.subDict("surfaceConformation")
    ),
    cellSizeControl_
    (
        *this,
        allGeometry_,
        cvMeshDict.subDict("motionControl")
    ),
    cvMeshControls_(*this, cvMeshDict),
    startOfInternalPoints_(0),
    startOfSurfacePointPairs_(0),
    featureVertices_(),
    featurePointLocations_(),
    featurePointTree_(),
    sizeAndAlignmentLocations_(),
    storedSizes_(),
    storedAlignments_(),
    sizeAndAlignmentTree_(),
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
            *this
        )
    ),
    faceAreaWeightModel_
    (
        faceAreaWeightModel::New
        (
            cvMeshDict.subDict("motionControl"),
            *this
        )
    )
{
    createFeaturePoints();

    if (cvMeshControls().objOutput())
    {
        geometryToConformTo_.writeFeatureObj("cvMesh");
    }

    insertInitialPoints();

    buildSurfaceConformation(rmCoarse);

    if(cvMeshControls().objOutput())
    {
        writePoints("allInitialPoints.obj", false);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::~conformalVoronoiMesh()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tensor Foam::conformalVoronoiMesh::requiredAlignment
(
    const point& pt
) const
{
    pointIndexHit surfHit;
    label hitSurface;

    geometryToConformTo_.findSurfaceNearest
    (
        pt,
        cvMeshControls().spanSqr(),
        surfHit,
        hitSurface
    );

    if (!surfHit.hit())
    {
        FatalErrorIn("conformalVoronoiMesh::requiredAlignment")
            << "findSurfaceNearest did not find a hit across the span of the "
            << "surfaces."
            << nl << exit(FatalError) << endl;
    }

    // Primary alignment

    vectorField norm(1);

    allGeometry_[hitSurface].getNormal
    (
        List<pointIndexHit>(1, surfHit),
        norm
    );

    vector np = norm[0];

    // Generate equally spaced 'spokes' in a circle normal to the
    // direction from the vertex to the closest point on the surface
    // and look for a secondary intersection.

    vector d = surfHit.hitPoint() - pt;

    tensor Rp = rotationTensor(vector(0,0,1), np);

    label s = cvMeshControls().alignmentSearchSpokes();

    scalar closestSpokeHitDistance = GREAT;

    pointIndexHit closestSpokeHit;

    label closestSpokeSurface = -1;

    for(label i = 0; i < s; i++)
    {
        vector spoke
        (
            Foam::cos(i*constant::mathematical::twoPi/s),
            Foam::sin(i*constant::mathematical::twoPi/s),
            0
        );

        spoke *= cvMeshControls().span();

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

        point mirrorPt = pt + 2*d;

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
        FatalErrorIn("conformalVoronoiMesh::requiredAlignment")
            << "No secondary surface hit found in spoke search."
            << nl << exit(FatalError);
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

    if (mag(ns) <  SMALL)
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


void Foam::conformalVoronoiMesh::insertSurfacePointPairs
(
    const List<pointIndexHit>& surfaceHits,
    const List<label>& hitSurfaces,
    const fileName fName
)
{
    if (surfaceHits.size() != hitSurfaces.size())
    {
        FatalErrorIn("Foam::conformalVoronoiMesh::insertPointPairs")
            << "surfaceHits and hitSurfaces are not the same size. Sizes "
            << surfaceHits.size() << ' '
            << hitSurfaces.size()
            << exit(FatalError);
    }

    forAll(surfaceHits, i)
    {
        vectorField norm(1);

        allGeometry_[hitSurfaces[i]].getNormal
        (
            List<pointIndexHit>(1, surfaceHits[i]),
            norm
        );

        const vector& normal = norm[0];

        const point& surfacePt(surfaceHits[i].hitPoint());

        insertPointPair
        (
            pointPairDistance(surfacePt),
            surfacePt,
            normal
        );
    }

    if(cvMeshControls().objOutput() && fName != fileName::null)
    {
        List<point> surfacePts(surfaceHits.size());

        forAll(surfacePts, i)
        {
            surfacePts[i] = surfaceHits[i].hitPoint();
        }

        writePoints(fName, surfacePts);
    }
}


void Foam::conformalVoronoiMesh::insertEdgePointGroups
(
    const List<pointIndexHit>& edgeHits,
    const labelList& featuresHit,
    const fileName fName
)
{
    if (edgeHits.size() != featuresHit.size())
    {
        FatalErrorIn("Foam::conformalVoronoiMesh::insertEdgePointGroups")
            << "edgeHits and featuresHit are not the same size. Sizes "
            << edgeHits.size() << ' '
            << featuresHit.size()
            << exit(FatalError);
    }

    forAll(edgeHits, i)
    {
        const featureEdgeMesh& feMesh
        (
            geometryToConformTo_.features()[featuresHit[i]]
        );

        insertEdgePointGroup(feMesh, edgeHits[i]);
    }

    if(cvMeshControls().objOutput() && fName != fileName::null)
    {
        List<point> edgePts(edgeHits.size());

        forAll(edgePts, i)
        {
            edgePts[i] = edgeHits[i].hitPoint();
        }

        writePoints(fName, edgePts);
    }
}


void Foam::conformalVoronoiMesh::insertEdgePointGroup
(
    const featureEdgeMesh& feMesh,
    const pointIndexHit& edHit
)
{
    label edgeI = edHit.index();

    featureEdgeMesh::edgeStatus edStatus = feMesh.getEdgeStatus(edgeI);

    if (edStatus == featureEdgeMesh::EXTERNAL)
    {
        insertExternalEdgePointGroup(feMesh, edHit);
    }
    else if (edStatus == featureEdgeMesh::INTERNAL)
    {
        insertInternalEdgePointGroup(feMesh, edHit);
    }
    else if (edStatus == featureEdgeMesh::FLAT)
    {
        insertFlatEdgePointGroup(feMesh, edHit);
    }
    else if (edStatus == featureEdgeMesh::OPEN)
    {
        insertOpenEdgePointGroup(feMesh, edHit);
    }
    else if (edStatus == featureEdgeMesh::MULTIPLE)
    {
        insertMultipleEdgePointGroup(feMesh, edHit);
    }
}


void Foam::conformalVoronoiMesh::insertExternalEdgePointGroup
(
    const featureEdgeMesh& feMesh,
    const pointIndexHit& edHit
)
{
    const point& edgePt = edHit.hitPoint();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    // As this is an external edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    // Convex. So refPt will be inside domain and hence a master point
    point refPt = edgePt - ppDist*(nA + nB)/(1 + (nA & nB) + VSMALL);

    // Insert the master point referring the the first slave
    label masterPtIndex = insertPoint(refPt, number_of_vertices() + 1);

    // Insert the slave points by reflecting refPt in both faces.
    // with each slave refering to the master

    point reflectedA = refPt + 2*ppDist*nA;
    insertPoint(reflectedA, masterPtIndex);

    point reflectedB = refPt + 2*ppDist*nB;
    insertPoint(reflectedB, masterPtIndex);
}


void Foam::conformalVoronoiMesh::insertInternalEdgePointGroup
(
    const featureEdgeMesh& feMesh,
    const pointIndexHit& edHit
)
{
    const point& edgePt = edHit.hitPoint();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    // As this is an external edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    // Concave. master and reflected points inside the domain.
    point refPt = edgePt - ppDist*(nA + nB)/(1 + (nA & nB) + VSMALL);

    // Generate reflected master to be outside.
    point reflMasterPt = refPt + 2*(edgePt - refPt);

    // Reflect reflMasterPt in both faces.
    point reflectedA = reflMasterPt - 2*ppDist*nA;

    point reflectedB = reflMasterPt - 2*ppDist*nB;

    scalar totalAngle =
        radToDeg(constant::mathematical::pi + acos(mag(nA & nB)));

    // Number of quadrants the angle should be split into
    int nQuads = int(totalAngle/cvMeshControls().maxQuadAngle()) + 1;

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
        point reflectedAa = refPt + ppDist*nB;
        insertPoint(reflectedAa, reflectedMaster);

        point reflectedBb = refPt + ppDist*nA;
        insertPoint(reflectedBb, reflectedMaster);
    }

    // Slave is outside.
    insertPoint(reflMasterPt, reflectedAI);
}


void Foam::conformalVoronoiMesh::insertFlatEdgePointGroup
(
    const featureEdgeMesh& feMesh,
    const pointIndexHit& edHit
)
{
    Info<< "    NOT INSERTING FLAT EDGE POINT GROUP, NOT IMPLEMENTED" << endl;
}


void Foam::conformalVoronoiMesh::insertOpenEdgePointGroup
(
    const featureEdgeMesh& feMesh,
    const pointIndexHit& edHit
)
{
    Info<< "    NOT INSERTING OPEN EDGE POINT GROUP, NOT IMPLEMENTED" << endl;
}


void Foam::conformalVoronoiMesh::insertMultipleEdgePointGroup
(
    const featureEdgeMesh& feMesh,
    const pointIndexHit& edHit
)
{
    Info<< "    NOT INSERTING MULTIPLE EDGE POINT GROUP, NOT IMPLEMENTED"
        << endl;
}


void Foam::conformalVoronoiMesh::createFeaturePoints()
{
    Info<< nl << "Creating bounding points" << endl;

    scalar bigSpan = 10*cvMeshControls().span();

    insertPoint(point(-bigSpan, -bigSpan, -bigSpan), Vb::ptFarPoint);
    insertPoint(point(-bigSpan, -bigSpan,  bigSpan), Vb::ptFarPoint);
    insertPoint(point(-bigSpan,  bigSpan, -bigSpan), Vb::ptFarPoint);
    insertPoint(point(-bigSpan,  bigSpan,  bigSpan), Vb::ptFarPoint);
    insertPoint(point( bigSpan, -bigSpan, -bigSpan), Vb::ptFarPoint);
    insertPoint(point( bigSpan, -bigSpan,  bigSpan), Vb::ptFarPoint);
    insertPoint(point( bigSpan,  bigSpan, -bigSpan), Vb::ptFarPoint);
    insertPoint(point( bigSpan,  bigSpan , bigSpan), Vb::ptFarPoint);

    Info<< nl << "Conforming to feature points" << endl;

    insertConvexFeaturesPoints();

    insertConcaveFeaturePoints();

    insertMixedFeaturePoints();

    Info<< "    Inserted " << number_of_vertices() << " vertices" << endl;

    featureVertices_.setSize(number_of_vertices());

    label featPtI = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        featureVertices_[featPtI] = Vb(vit->point());

        featureVertices_[featPtI].index() = vit->index();

        featureVertices_[featPtI].type() = vit->type();

        featPtI++;
    }

    constructFeaturePointLocations();

    if(cvMeshControls().objOutput())
    {
        writePoints("featureVertices.obj", false);
    }
}


void Foam::conformalVoronoiMesh::insertConvexFeaturesPoints()
{
    const PtrList<featureEdgeMesh>& feMeshes(geometryToConformTo_.features());

    forAll(feMeshes, i)
    {
        const featureEdgeMesh& feMesh(feMeshes[i]);

        for
        (
            label ptI = feMesh.convexStart();
            ptI < feMesh.concaveStart();
            ptI++
        )
        {
            vectorField featPtNormals = feMesh.featurePointNormals(ptI);

            const point& featPt = feMesh.points()[ptI];

            vector cornerNormal = sum(featPtNormals);
            cornerNormal /= mag(cornerNormal);

            point internalPt =  featPt - pointPairDistance(featPt)*cornerNormal;

            label internalPtIndex =
                insertPoint(internalPt, number_of_vertices() + 1);

            forAll (featPtNormals, nI)
            {
                const vector& n = featPtNormals[nI];

                plane planeN = plane(featPt, n);

                point externalPt =
                    internalPt + 2.0 * planeN.distance(internalPt) * n;

                insertPoint(externalPt, internalPtIndex);
            }
        }
    }
}


void Foam::conformalVoronoiMesh::insertConcaveFeaturePoints()
{
    const PtrList<featureEdgeMesh>& feMeshes(geometryToConformTo_.features());

    forAll(feMeshes, i)
    {
        const featureEdgeMesh& feMesh(feMeshes[i]);

        for
        (
            label ptI = feMesh.concaveStart();
            ptI < feMesh.mixedStart();
            ptI++
        )
        {
            vectorField featPtNormals = feMesh.featurePointNormals(ptI);

            const point& featPt = feMesh.points()[ptI];

            vector cornerNormal = sum(featPtNormals);
            cornerNormal /= mag(cornerNormal);

            point externalPt = featPt + pointPairDistance(featPt)*cornerNormal;

            label externalPtIndex = number_of_vertices() + featPtNormals.size();

            label internalPtIndex = -1;

            forAll (featPtNormals, nI)
            {
                const vector& n = featPtNormals[nI];

                plane planeN = plane(featPt, n);

                point internalPt =
                    externalPt - 2.0 * planeN.distance(externalPt) * n;

                internalPtIndex = insertPoint(internalPt, externalPtIndex);
            }

            insertPoint(externalPt, internalPtIndex);
        }
    }
}


void Foam::conformalVoronoiMesh::insertMixedFeaturePoints()
{
    const PtrList<featureEdgeMesh>& feMeshes(geometryToConformTo_.features());

    forAll(feMeshes, i)
    {
        const featureEdgeMesh& feMesh(feMeshes[i]);

        for
        (
            label ptI = feMesh.mixedStart();
            ptI < feMesh.nonFeatureStart();
            ptI++
        )
        {
            labelList pEds(feMesh.pointEdges()[ptI]);

            // Skipping unsupported mixed feature point types

            bool skipEdge = false;

            forAll(pEds, e)
            {
                label edgeI = pEds[e];

                featureEdgeMesh::edgeStatus edStatus =
                    feMesh.getEdgeStatus(edgeI);

                if
                (
                    edStatus == featureEdgeMesh::FLAT
                 || edStatus == featureEdgeMesh::OPEN
                 || edStatus == featureEdgeMesh::MULTIPLE
                )
                {
                    Info<< "    Edge type " << edStatus
                        << " found for mixed feature point " << ptI
                        << ". Not supported."
                        << endl;

                    skipEdge = true;
                }

            }

            if(skipEdge)
            {
                Info<< "    Skipping point " << ptI << nl << endl;

                continue;
            }

            // Inserting mixed internal and external feature points

            const point& pt(feMesh.points()[ptI]);

            scalar edgeGroupDistance = mixedFeaturePointDistance(pt);

            forAll(pEds, e)
            {
                label edgeI = pEds[e];

                point edgePt =
                    pt + edgeGroupDistance*feMesh.edgeDirection(edgeI, ptI);

                pointIndexHit edgeHit(true, edgePt, edgeI);

                insertEdgePointGroup(feMesh, edgeHit);
            }
        }
    }
}


void Foam::conformalVoronoiMesh::constructFeaturePointLocations()
{
    DynamicList<point> ftPtLocs;

    const PtrList<featureEdgeMesh>& feMeshes(geometryToConformTo_.features());

    forAll(feMeshes, i)
    {
        const featureEdgeMesh& feMesh(feMeshes[i]);

        for
        (
            label ptI = feMesh.convexStart();
            ptI < feMesh.nonFeatureStart();
            ptI++
        )
        {
            ftPtLocs.append(feMesh.points()[ptI]);
        }
    }

    featurePointLocations_.transfer(ftPtLocs);
}


void Foam::conformalVoronoiMesh::reinsertFeaturePoints()
{
    Info<< nl << "    Reinserting stored feature points" << endl;

    forAll(featureVertices_, f)
    {
        insertVb(featureVertices_[f]);
    }
}


const Foam::indexedOctree<Foam::treeDataPoint>&
Foam::conformalVoronoiMesh::featurePointTree() const
{
    if (featurePointTree_.empty())
    {
        treeBoundBox overallBb(geometryToConformTo_.bounds());

        Random rndGen(92561);

        overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        featurePointTree_.reset
        (
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint(featurePointLocations_),
                overallBb,  // overall search domain
                10,         // max levels
                10.0,       // maximum ratio of cubes v.s. cells
                100.0       // max. duplicity; n/a since no bounding boxes.
            )
        );
    }

    return featurePointTree_();
}


bool Foam::conformalVoronoiMesh::nearFeaturePt(const point& pt) const
{
    const indexedOctree<treeDataPoint>& tree = featurePointTree();

    scalar exclusionRangeSqr = featurePointExclusionDistanceSqr(pt);

    pointIndexHit info = tree.findNearest(pt, exclusionRangeSqr);

    return info.hit();
}


void Foam::conformalVoronoiMesh::insertInitialPoints()
{
    startOfInternalPoints_ = number_of_vertices();

    Info<< nl << "Inserting initial points" << endl;

    std::vector<Point> initPts = initialPointsMethod_->initialPoints();

    insertPoints(initPts);

    if(cvMeshControls().objOutput())
    {
        writePoints("initialPoints.obj", true);
    }

    storeSizesAndAlignments(initPts);
}


void Foam::conformalVoronoiMesh::storeSizesAndAlignments
(
    const std::vector<Point>& initPts
)
{
    timeCheck();

    Info << nl << "    Initialise stored size and alignment data" << endl;

    sizeAndAlignmentLocations_.setSize(initPts.size());

    storedSizes_.setSize(sizeAndAlignmentLocations_.size());

    storedAlignments_.setSize(sizeAndAlignmentLocations_.size());

    forAll(sizeAndAlignmentLocations_, i)
    {
        sizeAndAlignmentLocations_[i] = topoint(initPts[i]);

        storedSizes_[i] = cellSizeControl().cellSize
        (
            sizeAndAlignmentLocations_[i],
            false
        );

        storedAlignments_[i] = requiredAlignment(sizeAndAlignmentLocations_[i]);
    }

    timeCheck();

    buildSizeAndAlignmentTree();

    timeCheck();
}


const Foam::indexedOctree<Foam::treeDataPoint>&
Foam::conformalVoronoiMesh::sizeAndAlignmentTree() const
{
    if (sizeAndAlignmentTree_.empty())
    {
        buildSizeAndAlignmentTree();
    }

    return sizeAndAlignmentTree_();
}

void Foam::conformalVoronoiMesh::setVertexSizeAndAlignment()
{
    Info<< nl << "    Looking up target cell alignment and size" << endl;

    scalar spanSqr = cvMeshControls().spanSqr();

    const indexedOctree<treeDataPoint>& tree = sizeAndAlignmentTree();

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            point pt(topoint(vit->point()));

            pointIndexHit info = tree.findNearest(pt, spanSqr);

            vit->alignment() = storedAlignments_[info.index()];

            vit->targetCellSize() = storedSizes_[info.index()];
        }
    }

    // Info<< nl << "    Calculating target cell alignment and size" << endl;

    // for
    // (
    //     Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
    //     vit != finite_vertices_end();
    //     vit++
    // )
    // {
    //     if (vit->internalOrBoundaryPoint())
    //     {
    //         point pt(topoint(vit->point()));

    //         vit->alignment() = requiredAlignment(pt);

    //         vit->targetCellSize() = targetCellSize(pt);
    //     }
    // }
}


Foam::face Foam::conformalVoronoiMesh::buildDualFace
(
    const Triangulation::Finite_edges_iterator& eit
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
        label cc1I = cc1->cellIndex();

        label cc2I = cc2->cellIndex();

        if (cc1I < 0 || cc2I < 0)
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


Foam::scalar Foam::conformalVoronoiMesh::minFilterLimit
(
    const Triangulation::Finite_edges_iterator& eit
) const
{
    Cell_circulator ccStart = incident_cells(*eit);
    Cell_circulator cc = ccStart;

    scalar minFilterLimit = GREAT;

    do
    {
        if (cc->cellIndex() < 0)
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

        if (cc->filterLimit() < minFilterLimit)
        {
            minFilterLimit = cc->filterLimit();
        }

        cc++;

    } while (cc != ccStart);

    return minFilterLimit;
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

    if (!vA->internalOrBoundaryPoint())
    {
        dualCellIndexA = -1;
    }

    label dualCellIndexB = vB->index();

    if (!vB->internalOrBoundaryPoint())
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


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::move()
{
    timeCheck();

    scalar relaxation = relaxationModel_->relaxation();

    Info<< nl << "    Relaxation = " << relaxation << endl;

    pointField dualVertices(number_of_cells());

    label dualVertI = 0;

    // Find the dual point of each tetrahedron and assign it an index.
    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        cit->cellIndex() = -1;

        if
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
        {
            cit->cellIndex() = dualVertI;

            dualVertices[dualVertI] = topoint(dual(cit));

            dualVertI++;
        }
    }

    dualVertices.setSize(dualVertI);

    timeCheck();

    setVertexSizeAndAlignment();

    timeCheck();

    Info<< nl << "    Determining vertex displacements" << endl;

    vectorField cartesianDirections(3);

    cartesianDirections[0] = vector(0,0,1);
    cartesianDirections[1] = vector(0,1,0);
    cartesianDirections[2] = vector(1,0,0);

    vectorField displacementAccumulator
    (
        startOfSurfacePointPairs_,
        vector::zero
    );

    PackedBoolList pointToBeRetained(startOfSurfacePointPairs_, true);

    std::vector<Point> pointsToInsert;

    label pointsAdded = 0;

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        if
        (
            eit->first->vertex(eit->second)->internalOrBoundaryPoint()
         && eit->first->vertex(eit->third)->internalOrBoundaryPoint()
        )
        {
            face dualFace = buildDualFace(eit);

            Cell_handle c = eit->first;
            Vertex_handle vA = c->vertex(eit->second);
            Vertex_handle vB = c->vertex(eit->third);

            point dVA = topoint(vA->point());
            point dVB = topoint(vB->point());

            Field<vector> alignmentDirsA =
                vA->alignment() & cartesianDirections;
            Field<vector> alignmentDirsB =
                vB->alignment() & cartesianDirections;

            Field<vector> alignmentDirs(3);

            forAll(alignmentDirsA, aA)
            {
                const vector& a(alignmentDirsA[aA]);

                scalar maxDotProduct = 0.0;

                forAll(alignmentDirsB, aB)
                {
                    const vector& b(alignmentDirsB[aB]);

                    scalar dotProduct = a & b;

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

            forAll(alignmentDirs, aD)
            {
                vector& alignmentDir = alignmentDirs[aD];

                if ((rAB & alignmentDir) < 0)
                {
                    // swap the direction of the alignment so that has the
                    // same sense as rAB
                    alignmentDir *= -1;
                }

                scalar alignmentDotProd = ((rAB/rABMag) & alignmentDir);

                if
                (
                    alignmentDotProd
                  > cvMeshControls().cosAlignmentAcceptanceAngle()
                )
                {
                    scalar targetCellSize = averageCellSize(vA, vB);

                    scalar targetFaceArea = sqr(targetCellSize);

                    alignmentDir *= 0.5*targetCellSize;

                    vector delta = alignmentDir - 0.5*rAB;

                    scalar faceArea = dualFace.mag(dualVertices);

                    delta *= faceAreaWeightModel_->faceAreaWeight
                    (
                        faceArea/targetFaceArea
                    );

                    if
                    (
                        vA->internalPoint()
                     && vB->internalPoint()
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
                            // Prevent insertions spanning surfaces

                            pointsToInsert.push_back
                            (
                                toPoint(0.5*(dVA + dVB))
                            );

                            pointsAdded++;
                        }
                    }
                    else if
                    (
                        (vA->internalPoint() || vB->internalPoint())
                     && rABMag
                          < cvMeshControls().removalDistCoeff()*targetCellSize
                    )
                    {
                        // Point removal

                        // Only insert a point at the midpoint of the short edge
                        // if neither attached point has already been identified
                        // to be removed.
                        if
                        (
                            pointToBeRetained[vA->index()] == true
                         && pointToBeRetained[vB->index()] == true
                        )
                        {
                            pointsToInsert.push_back
                            (
                                toPoint(0.5*(dVA + dVB))
                            );
                        }

                        if (vA->internalPoint())
                        {
                            pointToBeRetained[vA->index()] = false;
                        }

                        if (vB->internalPoint())
                        {
                            pointToBeRetained[vB->index()] = false;
                        }
                    }
                    else
                    {
                        if (vA->internalPoint())
                        {
                            displacementAccumulator[vA->index()] += delta;
                        }

                        if (vB->internalPoint())
                        {
                            displacementAccumulator[vB->index()] += -delta;
                        }
                    }
                }
            }
        }
    }

    // Limit displacements that pierce, or get too close to the surface
    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            limitDisplacement
            (
                vit,
                displacementAccumulator[vit->index()]
            );
        }
    }

    vector totalDisp = sum(displacementAccumulator);
    scalar totalDist = sum(mag(displacementAccumulator));

    // Relax the calculated displacement
    displacementAccumulator *= relaxation;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            if (pointToBeRetained[vit->index()] == true)
            {
                pointsToInsert.push_back
                (
                    vit->point()
                  + toCGALVector(displacementAccumulator[vit->index()])
                );
            }
        }
    }

    // Write the mesh before clearing it.  Beware that writeMesh destroys the
    // indexing of the tessellation.
    if (runTime_.outputTime())
    {
        writeMesh(false);

        writeTargetCellSize();
    }

    // Remove the entire tessellation
    this->clear();

    reinsertFeaturePoints();

    startOfInternalPoints_ = number_of_vertices();

    timeCheck();

    Info<< nl << "    Inserting displaced tessellation" << endl;

    insertPoints(pointsToInsert);

    startOfSurfacePointPairs_ = number_of_vertices();

    label pointsRemoved =
        displacementAccumulator.size()
      - number_of_vertices()
      + pointsAdded;

    timeCheck();

    conformToSurface();

    timeCheck();

    Info<< nl
        << "    Total displacement = " << totalDisp << nl
        << "    Total distance = " << totalDist << nl
        << "    Points added = " << pointsAdded << nl
        << "    Points removed = " << pointsRemoved
        << endl;
}


// ************************************************************************* //
