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


void Foam::conformalVoronoiMesh::conformToSurface()
{
    reconformationMode reconfMode = reconformationControl();

    if (reconfMode == rmNone)
    {
        // Reinsert stored surface conformation
        reinsertSurfaceConformation();
    }
    else
    {
        // Rebuild, insert and store new surface conformation
        buildSurfaceConformation(reconfMode);
    }
}

Foam::conformalVoronoiMesh::reconformationMode
Foam::conformalVoronoiMesh::reconformationControl() const
{
    if (!runTime_.run())
    {
        Info<< nl << "    Rebuilding surface conformation for final output"
            << endl;

        return rmFine;
    }
    else if
    (
        runTime_.timeIndex()
      % cvMeshControls().surfaceConformationRebuildFrequency()
     == 0
    )
    {
        Info<< nl << "    Rebuilding surface conformation for more iterations"
            << endl;

        return rmCoarse;
    }

    return rmNone;
}

void Foam::conformalVoronoiMesh::buildSurfaceConformation
(
    reconformationMode reconfMode
)
{
    if (reconfMode == rmCoarse)
    {
        Info<< nl << "    Build coarse surface conformation" << endl;
    }
    else if (reconfMode == rmFine)
    {
        Info<< nl << "    Build fine surface conformation" << endl;
    }
    else if (reconfMode == rmNone)
    {
        WarningIn("buildSurfaceConformation(reconformationMode reconfMode)")
            << "reconformationMode rmNone specified, not building conformation"
            << endl;

        return;
    }
    else
    {
        WarningIn("buildSurfaceConformation(reconformationMode reconfMode)")
            << "Unknown reconformationMode " << reconfMode << endl;

        return;
    }

    timeCheck();

    startOfSurfacePointPairs_ = number_of_vertices();

    // Initialise containers to store the edge conformation locations
    DynamicList<point> newEdgeLocations;

    pointField existingEdgeLocations(0);

    autoPtr<indexedOctree<treeDataPoint> > edgeLocationTree;

    // Initialise the edgeLocationTree
    buildEdgeLocationTree(edgeLocationTree, existingEdgeLocations);

    label initialTotalHits = 0;

    // Initial surface protrusion conformation - nearest surface point
    {
        Info<< "    EDGE DISTANCE COEFFS HARD-CODED." << endl;
        scalar edgeSearchDistCoeffSqr = sqr(1.1);
        scalar surfacePtReplaceDistCoeffSqr = sqr(0.5);

        DynamicList<pointIndexHit> surfaceHits;
        DynamicList<label> hitSurfaces;

        DynamicList<pointIndexHit> featureEdgeHits;
        DynamicList<label> featureEdgeFeaturesHit;

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
                scalar searchDistanceSqr = surfaceSearchDistanceSqr(vert);
                pointIndexHit surfHit;
                label hitSurface;

                geometryToConformTo_.findSurfaceNearest
                (
                    vert,
                    searchDistanceSqr,
                    surfHit,
                    hitSurface
                );

                if (surfHit.hit())
                {
                    vit->setNearBoundary();

                    if (dualCellSurfaceAnyIntersection(vit))
                    {
                        addSurfaceAndEdgeHits
                        (
                            vit,
                            vert,
                            surfHit,
                            hitSurface,
                            surfacePtReplaceDistCoeffSqr,
                            edgeSearchDistCoeffSqr,
                            surfaceHits,
                            hitSurfaces,
                            featureEdgeHits,
                            featureEdgeFeaturesHit,
                            newEdgeLocations,
                            existingEdgeLocations,
                            edgeLocationTree
                        );
                    }
                }
            }
        }

        Info<< nl <<"    Initial conformation" << nl
            << "        Number of vertices " << number_of_vertices() << nl
            << "        Number of surface hits " << surfaceHits.size() << nl
            << "        Number of edge hits " << featureEdgeHits.size()
            << endl;

        insertSurfacePointPairs
        (
            surfaceHits,
            hitSurfaces,
            "surfaceConformationLocations_initial.obj"
        );

        insertEdgePointGroups
        (
            featureEdgeHits,
            featureEdgeFeaturesHit,
            "edgeConformationLocations_initial.obj"
        );

        initialTotalHits = surfaceHits.size() + featureEdgeHits.size();
    }

    label iterationNo = 0;

    label maxIterations = 10;

    Info<< nl << "    MAX ITERATIONS HARD CODED TO "<< maxIterations << endl;

    scalar iterationToIntialHitRatioLimit = 0.01;

    label hitLimit = label(iterationToIntialHitRatioLimit*initialTotalHits);

    Info<< "    STOPPING ITERATIONS WHEN TOTAL NUMBER OF HITS DROPS BELOW "
        << iterationToIntialHitRatioLimit << " (HARD CODED) OF INITIAL HITS ("
        << hitLimit << ")"
        << endl;

    // Set totalHits to a large enough positive value to enter the while loop on
    // the first iteration
    label totalHits = initialTotalHits;

    while
    (
        totalHits > 0
     && totalHits > hitLimit
     && iterationNo < maxIterations
    )
    {
        Info<< "    EDGE DISTANCE COEFFS HARD-CODED." << endl;
        scalar edgeSearchDistCoeffSqr = sqr(1.25);
        scalar surfacePtReplaceDistCoeffSqr = sqr(0.7);

        DynamicList<pointIndexHit> surfaceHits;
        DynamicList<label> hitSurfaces;

        DynamicList<pointIndexHit> featureEdgeHits;
        DynamicList<label> featureEdgeFeaturesHit;

        for
        (
            Triangulation::Finite_vertices_iterator vit =
            finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            // The initial surface conformation has already identified the
            // nearBoundary set of vertices.  Previously inserted boundary
            // points can also generate protrusions and must be assessed too.

            if (vit->nearBoundary() || vit->ppMaster())
            {
                point vert(topoint(vit->point()));
                pointIndexHit surfHit;
                label hitSurface;

                dualCellLargestSurfaceProtrusion(vit, surfHit, hitSurface);

                if (surfHit.hit())
                {
                    addSurfaceAndEdgeHits
                    (
                        vit,
                        vert,
                        surfHit,
                        hitSurface,
                        surfacePtReplaceDistCoeffSqr,
                        edgeSearchDistCoeffSqr,
                        surfaceHits,
                        hitSurfaces,
                        featureEdgeHits,
                        featureEdgeFeaturesHit,
                        newEdgeLocations,
                        existingEdgeLocations,
                        edgeLocationTree
                    );
                }
            }
        }

        Info<< nl <<"    Conformation iteration " << iterationNo << nl
            << "        Number of vertices " << number_of_vertices() << nl
            << "        Number of surface hits " << surfaceHits.size() << nl
            << "        Number of edge hits " << featureEdgeHits.size()
            << endl;

        totalHits = surfaceHits.size() + featureEdgeHits.size();

        if (totalHits > 0)
        {
            insertSurfacePointPairs
            (
                surfaceHits,
                hitSurfaces,
                fileName
                (
                    "surfaceConformationLocations_" + name(iterationNo) + ".obj"
                )
            );

            insertEdgePointGroups
            (
                featureEdgeHits,
                featureEdgeFeaturesHit,
                "edgeConformationLocations_" + name(iterationNo) + ".obj"
            );
        }

        iterationNo++;

        if (iterationNo == maxIterations)
        {
            WarningIn("conformalVoronoiMesh::conformToSurface()")
                << "Maximum surface conformation iterations ("
                << maxIterations <<  ") reached." << endl;
        }

        if (totalHits < hitLimit)
        {
            Info<< nl << "    Total hits (" << totalHits
                << ") less than limit (" << hitLimit
                << "), stopping iterations" << endl;
        }
    }

    // Info<< nl << "    After iterations, check penetrations" << endl;

    // for
    // (
    //     Triangulation::Finite_vertices_iterator vit =
    //     finite_vertices_begin();
    //     vit != finite_vertices_end();
    //     vit++
    // )
    // {
    //     if (vit->internalOrBoundaryPoint())
    //     {
    //         point vert(topoint(vit->point()));
    //         pointIndexHit surfHit;
    //         label hitSurface;

    //         dualCellLargestSurfaceProtrusion(vit, surfHit, hitSurface);

    //         if (surfHit.hit())
    //         {
    //             Info<< nl << "Residual penetration: " << nl
    //                 << vit->index() << nl
    //                 << vit->type() << nl
    //                 << vit->ppMaster() << nl
    //                 << "nearFeaturePt "
    //                 << nearFeaturePt(surfHit.hitPoint()) << nl
    //                 << vert << nl
    //                 << surfHit.hitPoint()
    //                 << endl;
    //         }
    //     }
    // }

    storeSurfaceConformation();
}


bool Foam::conformalVoronoiMesh::dualCellSurfaceAnyIntersection
(
    const Triangulation::Finite_vertices_iterator& vit
) const
{
    std::list<Facet> facets;
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
        if (!geometryToConformTo_.bounds().contains(dE0))
        {
            return true;
        }

        point dE1 = topoint(dual(fit->first->neighbor(fit->second)));

        // If other edge end is outside bounding box then edge cuts boundary
        if (!geometryToConformTo_.bounds().contains(dE1))
        {
            return true;
        }

        // Check for the edge passing through a surface
        if (geometryToConformTo_.findSurfaceAnyIntersection(dE0, dE1))
        {
            return true;
        }
    }

    return false;
}


void Foam::conformalVoronoiMesh::dualCellLargestSurfaceProtrusion
(
    const Triangulation::Finite_vertices_iterator& vit,
    pointIndexHit& surfHitLargest,
    label& hitSurfaceLargest
) const
{
    std::list<Facet> facets;
    incident_facets(vit, std::back_inserter(facets));

    point vert(topoint(vit->point()));

    scalar maxProtrusionDistance = maxSurfaceProtrusion(vert);

    for
    (
        std::list<Facet>::iterator fit=facets.begin();
        fit != facets.end();
        ++fit
    )
    {
        if
        (
            !is_infinite(fit->first)
         && !is_infinite(fit->first->neighbor(fit->second))
        )
        {
            point edgeMid =
                0.5
               *(
                    topoint(dual(fit->first))
                  + topoint(dual(fit->first->neighbor(fit->second)))
                );

            pointIndexHit surfHit;
            label hitSurface;

            geometryToConformTo_.findSurfaceAnyIntersection
            (
                vert,
                edgeMid,
                surfHit,
                hitSurface
            );

            if (surfHit.hit())
            {
                vectorField norm(1);

                allGeometry_[hitSurface].getNormal
                (
                    List<pointIndexHit>(1, surfHit),
                    norm
                );

                const vector& n = norm[0];

                scalar normalProtrusionDistance =
                    (edgeMid - surfHit.hitPoint()) & n;

                if (normalProtrusionDistance > maxProtrusionDistance)
                {
                    surfHitLargest = surfHit;
                    hitSurfaceLargest = hitSurface;

                    maxProtrusionDistance = normalProtrusionDistance;
                }
            }
        }
    }
}


void Foam::conformalVoronoiMesh::limitDisplacement
(
    const Triangulation::Finite_vertices_iterator& vit,
    vector& displacement
) const
{
    point pt = topoint(vit->point());
    point dispPt = pt + displacement;

    bool limit = false;

    pointIndexHit surfHit;
    label hitSurface;

    if (!geometryToConformTo_.bounds().contains(dispPt))
    {
        // If dispPt is outside bounding box then displacement cuts boundary
        limit = true;

        // Info<< "    bb limit" << endl;
    }
    else if (geometryToConformTo_.findSurfaceAnyIntersection(pt, dispPt))
    {
        // Full surface penetration test
        limit = true;

        // Info<< "    intersection limit" << endl;
    }
    else
    {
        // Testing if the displaced position is too close to the surface.
        // Within twice the local surface point pair insertion distance is
        // considered "too close"

        scalar searchDistanceSqr = sqr
        (
            2*vit->targetCellSize()
           *cvMeshControls().pointPairDistanceCoeff()
        );

        geometryToConformTo_.findSurfaceNearest
        (
            dispPt,
            searchDistanceSqr,
            surfHit,
            hitSurface
        );

        if (surfHit.hit())
        {
            // Info<< "    proximity limit" << endl;

            limit = true;

            if (magSqr(pt - surfHit.hitPoint()) <= searchDistanceSqr)
            {
                // Info<< "    Cannot limit displacement, point " << pt
                //     << " closer than tolerance" << endl;

                return;
            }
        }
    }

    if (limit)
    {
        // Halve the displacement and call this function again.  Will continue
        // recursively until the displacement is small enough.

        displacement *= 0.5;

        // Info<< "    Limiting displacement of point " << pt << endl;

        limitDisplacement(vit, displacement);
    }
}


bool Foam::conformalVoronoiMesh::nearFeatureEdgeLocation
(
    const point& pt,
    DynamicList<point>& newEdgeLocations,
    pointField& existingEdgeLocations,
    autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree
) const
{
    scalar exclusionRangeSqr = featureEdgeExclusionDistanceSqr(pt);

    // 0.01 and 1000 determined from speed tests, varying the indexedOctree
    // rebuild frequency and balance of additions to queries.

    if
    (
        newEdgeLocations.size()
     >= max(0.01*existingEdgeLocations.size(), 1000)
    )
    {
        existingEdgeLocations.append(newEdgeLocations);

        buildEdgeLocationTree(edgeLocationTree, existingEdgeLocations);

        newEdgeLocations.clear();
    }
    else
    {
        // Search for the nearest point in newEdgeLocations.
        // Searching here first, because the intention is that the value of
        // newEdgeLocationsSizeLimit should make this faster by design.

        if (min(magSqr(newEdgeLocations - pt)) <= exclusionRangeSqr)
        {
            return true;
        }
    }

    // Searching for the nearest point in existingEdgeLocations using the
    // indexedOctree

    pointIndexHit info = edgeLocationTree().findNearest(pt, exclusionRangeSqr);

    return info.hit();
}


void Foam::conformalVoronoiMesh::buildEdgeLocationTree
(
    autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree,
    const pointField& existingEdgeLocations
) const
{
    treeBoundBox overallBb(geometryToConformTo_.bounds());

    Random rndGen(72953);

    overallBb.extend(rndGen, 1E-4);
    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    edgeLocationTree.reset
    (
        new indexedOctree<treeDataPoint>
        (
            treeDataPoint(existingEdgeLocations),
            overallBb,  // overall search domain
            10,         // max levels
            10.0,       // maximum ratio of cubes v.s. cells
            100.0       // max. duplicity; n/a since no bounding boxes.
        )
    );
}


void Foam::conformalVoronoiMesh::buildSizeAndAlignmentTree() const
{
    treeBoundBox overallBb(geometryToConformTo_.bounds());

    Random rndGen(627391);

    overallBb.extend(rndGen, 1E-4);
    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    sizeAndAlignmentTree_.reset
    (
        new indexedOctree<treeDataPoint>
        (
            treeDataPoint(sizeAndAlignmentLocations_),
            overallBb,  // overall search domain
            10,         // max levels
            10.0,       // maximum ratio of cubes v.s. cells
            100.0       // max. duplicity; n/a since no bounding boxes.
        )
    );
}


void Foam::conformalVoronoiMesh::addSurfaceAndEdgeHits
(
    const Triangulation::Finite_vertices_iterator& vit,
    const point& vert,
    const pointIndexHit& surfHit,
    label hitSurface,
    scalar surfacePtReplaceDistCoeffSqr,
    scalar edgeSearchDistCoeffSqr,
    DynamicList<pointIndexHit>& surfaceHits,
    DynamicList<label>& hitSurfaces,
    DynamicList<pointIndexHit>& featureEdgeHits,
    DynamicList<label>& featureEdgeFeaturesHit,
    DynamicList<point>& newEdgeLocations,
    pointField& existingEdgeLocations,
    autoPtr<indexedOctree<treeDataPoint> >& edgeLocationTree
) const
{
    bool keepSurfacePoint = true;

    if (nearFeaturePt(surfHit.hitPoint()))
    {
        keepSurfacePoint = false;

        if (vit->index() < startOfInternalPoints_)
        {
            surfaceHits.append(surfHit);

            hitSurfaces.append(hitSurface);
        }
    }

    List<pointIndexHit> edHits;

    labelList featuresHit;

    scalar targetCellSizeSqr = sqr(targetCellSize(vert));

    geometryToConformTo_.findEdgeNearestByType
    (
        surfHit.hitPoint(),
        edgeSearchDistCoeffSqr*targetCellSizeSqr,
        edHits,
        featuresHit
    );

    // Gather edge locations but do not add them to newEdgeLocations inside the
    // loop as they will prevent nearby edge locations of different types being
    // conformed to.

    DynamicList<point> currentEdgeLocations;

    forAll(edHits, i)
    {
        const pointIndexHit& edHit(edHits[i]);

        label featureHit = featuresHit[i];

        if (edHit.hit())
        {
            if (!nearFeaturePt(edHit.hitPoint()))
            {
                if
                (
                    magSqr(edHit.hitPoint() - surfHit.hitPoint())
                  < surfacePtReplaceDistCoeffSqr*targetCellSizeSqr
                )
                {
                    // If the point is within a given distance of a feature
                    // edge, give control to edge control points instead, this
                    // will prevent "pits" forming.

                    keepSurfacePoint = false;
                }

                if
                (
                    !nearFeatureEdgeLocation
                    (
                        edHit.hitPoint(),
                        newEdgeLocations,
                        existingEdgeLocations,
                        edgeLocationTree
                    )
                )
                {
                    // Do not place edge control points too close to a feature
                    // point or existing edge control points

                    featureEdgeHits.append(edHit);
                    featureEdgeFeaturesHit.append(featureHit);

                    currentEdgeLocations.append(edHit.hitPoint());
                }
            }
        }
    }

    newEdgeLocations.append(currentEdgeLocations);

    if (keepSurfacePoint)
    {
        surfaceHits.append(surfHit);

        hitSurfaces.append(hitSurface);
    }
}


void Foam::conformalVoronoiMesh::storeSurfaceConformation()
{
    Info<< nl << "    Storing surface conformation" << endl;

    surfaceConformationVertices_.setSize
    (
        number_of_vertices() - startOfSurfacePointPairs_
    );

    label surfPtI = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->index() >= startOfSurfacePointPairs_)
        {
            if (!vit->pairPoint())
            {
                FatalErrorIn("storeSurfaceConformation()")
                    << "Trying to store a vertex that is not a surface point"
                    << exit(FatalError);
            }

            surfaceConformationVertices_[surfPtI] = Vb(vit->point());

            surfaceConformationVertices_[surfPtI].index() =
            vit->index() - startOfSurfacePointPairs_;

            surfaceConformationVertices_[surfPtI].type() =
            vit->type() - startOfSurfacePointPairs_;

            surfPtI++;
        }
    }

    Info<< "    Stored " << surfaceConformationVertices_.size()
        << " vertices" << endl;
}


void Foam::conformalVoronoiMesh::reinsertSurfaceConformation()
{
    Info<< nl << "    Reinserting stored surface conformation" << endl;

    startOfSurfacePointPairs_ = number_of_vertices();

    forAll(surfaceConformationVertices_, v)
    {
        insertVb(surfaceConformationVertices_[v], startOfSurfacePointPairs_);
    }

    Info<< "    Reinserted " << number_of_vertices() - startOfSurfacePointPairs_
        << " vertices" << endl;
}


void Foam::conformalVoronoiMesh::calcDualMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    timeCheck();

    setVertexSizeAndAlignment();

    timeCheck();

    // // ~~~~~~~~~~~ removing short edges by indexing dual vertices ~~~~~~~~~~~~~~

    // for
    // (
    //     Triangulation::Finite_cells_iterator cit = finite_cells_begin();
    //     cit != finite_cells_end();
    //     ++cit
    // )
    // {
    //     cit->cellIndex() = -1;
    // }

    // points.setSize(number_of_cells());

    // // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // label dualVertI = 0;

    // // Scanning by number of short (dual) edges (nSE) attached to the
    // // circumcentre of each Delaunay tet.  A Delaunay tet may only have four
    // // dual edges emanating from its circumcentre, assigning positions and
    // // indices to those with 4 short edges attached first, then >= 3, then >= 2
    // // etc.
    // for (label nSE = 4; nSE >= 0; nSE--)
    // {
    //     Info<< nl << "Scanning for dual vertices with >= "
    //         << nSE
    //         << " short edges attached." << endl;

    //     for
    //     (
    //         Triangulation::Finite_cells_iterator cit = finite_cells_begin();
    //         cit != finite_cells_end();
    //         ++cit
    //     )
    //     {
    //         // If the Delaunay tet has an index already then it has either
    //         // evaluated itself and taken action or has had its index dictated
    //         // by a neighbouring tet with more short edges attached.

    //         if (cit->cellIndex() == -1)
    //         {
    //             point dualVertex = topoint(dual(cit));

    //             label shortEdges = 0;

    //             List<bool> edgeIsShort(4, false);

    //             List<bool> neighbourAlreadyIndexed(4, false);

    //             // Loop over the four facets of the Delaunay tet
    //             for (label f = 0; f < 4; f++)
    //             {
    //                 // Check that at least one of the vertices of the facet is
    //                 // an internal or boundary point
    //                 if
    //                 (
    //                     cit->vertex(vertex_triple_index(f, 0))->
    //                         internalOrBoundaryPoint()
    //                  || cit->vertex(vertex_triple_index(f, 1))->
    //                         internalOrBoundaryPoint()
    //                  || cit->vertex(vertex_triple_index(f, 2))->
    //                         internalOrBoundaryPoint()
    //                 )
    //                 {
    //                     point neighDualVertex;

    //                     label cNI = cit->neighbor(f)->cellIndex();

    //                     if (cNI == -1)
    //                     {
    //                         neighDualVertex = topoint(dual(cit->neighbor(f)));
    //                     }
    //                     else
    //                     {
    //                         neighDualVertex = points[cNI];
    //                     }

    //                     if
    //                     (
    //                         magSqr(dualVertex - neighDualVertex)
    //                       < sqr
    //                         (
    //                             minimumEdgeLength
    //                             (
    //                                 0.5*(dualVertex + neighDualVertex)
    //                             )
    //                         )
    //                     )
    //                     {
    //                         edgeIsShort[f] = true;

    //                         if (cNI > -1)
    //                         {
    //                             neighbourAlreadyIndexed[f] = true;
    //                         }

    //                         shortEdges++;
    //                     }
    //                 }
    //             }

    //             if (nSE == 0 && shortEdges == 0)
    //             {
    //                 // Final iteration and no short edges are found, index
    //                 // remaining dual vertices.

    //                 if
    //                 (
    //                     cit->vertex(0)->internalOrBoundaryPoint()
    //                  || cit->vertex(1)->internalOrBoundaryPoint()
    //                  || cit->vertex(2)->internalOrBoundaryPoint()
    //                  || cit->vertex(3)->internalOrBoundaryPoint()
    //                 )
    //                 {
    //                     cit->cellIndex() = dualVertI;
    //                     points[dualVertI] = dualVertex;
    //                     dualVertI++;
    //                 }
    //             }
    //             else if
    //             (
    //                 shortEdges >= nSE
    //             )
    //             {
    //                 // Info<< neighbourAlreadyIndexed << ' '
    //                 //     << edgeIsShort << endl;

    //                 label numUnindexedNeighbours = 1;

    //                 for (label f = 0; f < 4; f++)
    //                 {
    //                     if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
    //                     {
    //                         dualVertex += topoint(dual(cit->neighbor(f)));

    //                         numUnindexedNeighbours++;
    //                     }
    //                 }

    //                 dualVertex /= numUnindexedNeighbours;

    //                 label nearestExistingIndex = -1;

    //                 point nearestIndexedNeighbourPos = vector::zero;

    //                 scalar minDistSqrToNearestIndexedNeighbour = VGREAT;

    //                 for (label f = 0; f < 4; f++)
    //                 {
    //                     if (edgeIsShort[f] && neighbourAlreadyIndexed[f])
    //                     {
    //                         label cNI = cit->neighbor(f)->cellIndex();

    //                         point indexedNeighbourPos = points[cNI];

    //                         if
    //                         (
    //                             magSqr(indexedNeighbourPos - dualVertex)
    //                           < minDistSqrToNearestIndexedNeighbour
    //                         )
    //                         {
    //                             nearestExistingIndex = cNI;

    //                             nearestIndexedNeighbourPos =
    //                             indexedNeighbourPos;

    //                             minDistSqrToNearestIndexedNeighbour =
    //                             magSqr(indexedNeighbourPos - dualVertex);
    //                         }
    //                     }
    //                 }

    //                 if
    //                 (
    //                     nearestExistingIndex > -1
    //                  && minDistSqrToNearestIndexedNeighbour
    //                   < sqr
    //                     (
    //                         minimumEdgeLength
    //                         (
    //                             0.5*(nearestIndexedNeighbourPos + dualVertex)
    //                         )
    //                     )
    //                 )
    //                 {
    //                     points[nearestExistingIndex] =
    //                     0.5*(dualVertex + nearestIndexedNeighbourPos);

    //                     for (label f = 0; f < 4; f++)
    //                     {
    //                         if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
    //                         {
    //                             cit->neighbor(f)->cellIndex() =
    //                             nearestExistingIndex;
    //                         }
    //                     }

    //                     cit->cellIndex() = nearestExistingIndex;
    //                 }
    //                 else
    //                 {
    //                     for (label f = 0; f < 4; f++)
    //                     {
    //                         if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
    //                         {
    //                             cit->neighbor(f)->cellIndex() = dualVertI;
    //                         }
    //                     }

    //                     cit->cellIndex() = dualVertI;

    //                     points[dualVertI] = dualVertex;

    //                     dualVertI++;
    //                 }
    //             }
    //         }
    //     }
    // }

    // points.setSize(dualVertI);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~ dual cell indexing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // assigns an index to the Delaunay vertices which will be the dual cell
    // index used for owner neighbour assignment.

    // The indices of the points are reset which destroys the point-pair
    // matching, so the type of each vertex are reset to avoid any ambiguity.

    label dualCelli = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            vit->type() = Vb::ptInternalPoint;
            vit->index() = dualCelli;
            dualCelli++;
        }
        else
        {
            vit->type() = Vb::ptFarPoint;
            vit->index() = -1;
        }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~ dual face filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Loop over all dual faces and merge points to remove faces that
    // are not wanted.

    // Indexing Delaunay cells, which are Dual vertices

    label dualVertI = 0;

    points.setSize(number_of_cells());

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
        {
            cit->cellIndex() = dualVertI;
            points[dualVertI] = topoint(dual(cit));
            dualVertI++;
        }
        else
        {
            cit->cellIndex() = -1;
        }
    }

    points.setSize(dualVertI);

    Info<< nl << "    Merging close points" << endl;

    label nCollapsedFaces = 0;

    label nPtsMerged = 0;

    do
    {
        Map<label> dualPtIndexMap;

        nPtsMerged = mergeCloseDualVertices(points, dualPtIndexMap);

        Info<< "        Merged " << nPtsMerged << " points" << endl;

        reindexDualVertices(dualPtIndexMap);

    } while (nPtsMerged > 0);

    Info<< nl << "    Collapsing unnecessary faces" << endl;

    do
    {
        Map<label> dualPtIndexMap;

        nCollapsedFaces = 0;

        nCollapsedFaces = collapseFaces(points, dualPtIndexMap);

        reindexDualVertices(dualPtIndexMap);

        Info<< "        Collapsed " << nCollapsedFaces << " faces" << endl;

    } while (nCollapsedFaces > 0);

    // ~~~~~~~~~~~~ dual face and owner neighbour construction ~~~~~~~~~~~~~~~~~

    patchNames = geometryToConformTo_.patchNames();

    patchNames.setSize(patchNames.size() + 1);

    patchNames[patchNames.size() - 1] = "cvMesh_defaultPatch";

    label nPatches = patchNames.size();

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_edges());

    owner.setSize(number_of_edges());

    neighbour.setSize(number_of_edges());

    label dualFaceI = 0;

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_handle c = eit->first;
        Vertex_handle vA = c->vertex(eit->second);
        Vertex_handle vB = c->vertex(eit->third);

        if
        (
            vA->internalOrBoundaryPoint()
         || vB->internalOrBoundaryPoint()
        )
        {
            face newDualFace = buildDualFace(eit);

            //bool keepFace = assessFace(newDualFace, vA, vB, points);

            // if (newDualFace.size() >= 3 && keepFace)
            // {

            if (newDualFace.size() >= 3)
            {

                label dcA = vA->index();

                if (!vA->internalOrBoundaryPoint())
                {
                    dcA = -1;
                }

                label dcB = vB->index();

                if (!vB->internalOrBoundaryPoint())
                {
                    dcB = -1;
                }

                label dcOwn = -1;
                label dcNei = -1;

                if (dcA == -1 && dcB == -1)
                {
                    FatalErrorIn("calcDualMesh")
                        << "Attempting to create a face joining "
                        << "two external dual cells "
                        << exit(FatalError);
                }
                else if (dcA == -1 || dcB == -1)
                {
                    // boundary face, find which is the owner

                    if (dcA == -1)
                    {
                        dcOwn = dcB;

                        // reverse face order to correctly orientate normal
                        reverse(newDualFace);
                    }
                    else
                    {
                        dcOwn = dcA;
                    }

                    // Find which patch this face is on by finding the
                    // intersection with the surface of the Delaunay edge
                    // generating the face and identify the region of the
                    // intersection.

                    point ptA = topoint(vA->point());

                    point ptB = topoint(vB->point());

                    label patchIndex = geometryToConformTo_.findPatch(ptA, ptB);

                    if (patchIndex == -1)
                    {
                        patchIndex = patchNames.size() - 1;

                        WarningIn("Foam::conformalVoronoiMesh::calcDualMesh")
                            << "Dual face found between Dv pair " << nl
                            << "    ptA" << nl
                            << "    ptB" << nl
                            << "    that is not on a surface patch. Adding to "
                            << patchNames[patchIndex]
                            << endl;
                    }

                    patchFaces[patchIndex].append(newDualFace);
                    patchOwners[patchIndex].append(dcOwn);
                }
                else
                {
                    // internal face, find the lower cell to be the owner

                    if (dcB > dcA)
                    {
                        dcOwn = dcA;
                        dcNei = dcB;
                    }
                    else
                    {
                        dcOwn = dcB;
                        dcNei = dcA;

                        // reverse face order to correctly orientate normal
                        reverse(newDualFace);
                    }

                    faces[dualFaceI] = newDualFace;

                    owner[dualFaceI] = dcOwn;

                    neighbour[dualFaceI] = dcNei;

                    dualFaceI++;
                }
            }
        }
    }

    label nInternalFaces = dualFaceI;

    faces.setSize(nInternalFaces);
    owner.setSize(nInternalFaces);
    neighbour.setSize(nInternalFaces);

    timeCheck();

    sortFaces(faces, owner, neighbour);

    timeCheck();

    addPatches
    (
        nInternalFaces,
        faces,
        owner,
        patchNames,
        patchSizes,
        patchStarts,
        patchFaces,
        patchOwners,
        false
    );

    removeUnusedPoints(faces, points);

    timeCheck();

    // // Write out faces to be removed as a list of labels to be used in
    // // faceSet

    // DynamicList<label> facesToBeRemoved;

    // labelList nEdgeHistogram(12, 0);

    // forAll(faces, fI)
    // {
    //     const face& f = faces[fI];

    //     if (!assessFace(f, targetCellSize(f.centre(points)), points))
    //     {
    //         facesToBeRemoved.append(fI);

    //         nEdgeHistogram[f.size()]++;
    //     }
    // }

    // fileName fName = "facesToBeRemoved";

    // OFstream str(fName);

    // str << facesToBeRemoved;

    // Info<< nEdgeHistogram << endl;
}


void Foam::conformalVoronoiMesh::calcTetMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    labelList vertexMap(number_of_vertices());

    label vertI = 0;

    points.setSize(number_of_vertices());

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() || vit->pairPoint())
        {
            vertexMap[vit->index()] = vertI;
            points[vertI] =  topoint(vit->point());
            vertI++;
        }
    }

    points.setSize(vertI);

    label cellI = 0;

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
        {
             cit->cellIndex() = cellI++;
        }
        else
        {
            cit->cellIndex() = -1;
        }
    }

    patchNames = geometryToConformTo_.patchNames();

    patchNames.setSize(patchNames.size() + 1);

    patchNames[patchNames.size() - 1] = "cvMesh_defaultPatch";

    label nPatches = patchNames.size();

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_facets());

    owner.setSize(number_of_facets());

    neighbour.setSize(number_of_facets());

    label faceI = 0;

    labelList verticesOnTriFace(3, -1);

    face newFace(verticesOnTriFace);

    for
    (
        Triangulation::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const int oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        label c1I = c1->cellIndex();
        label c2I = c2->cellIndex();

        label ownerCell = -1;
        label neighbourCell = -1;

        if (c1I == -1 && c2I == -1)
        {
            // Both tets are outside, skip
            continue;
        }

        for (label i = 0; i < 3; i++)
        {
            verticesOnTriFace[i] = vertexMap
            [
                c1->vertex(vertex_triple_index(oppositeVertex, i))->index()
            ];
        }

        newFace = face(verticesOnTriFace);

        if (c1I == -1 || c2I == -1)
        {
            // Boundary face...
            if (c1I == -1)
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
                    << "    " << newFace.centre(points) << nl
                    << "    did not find a surface patch. Adding to "
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

    addPatches
    (
        nInternalFaces,
        faces,
        owner,
        patchNames,
        patchSizes,
        patchStarts,
        patchFaces,
        patchOwners,
        false
    );
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
            verticesOnFace.append(cc1I);
        }

        cc1++;

        cc2++;

    } while (cc1 != ccStart);

    return face(verticesOnFace);
}


bool Foam::conformalVoronoiMesh::assessFace
(
    const face& f,
    const Vertex_handle& vA,
    const Vertex_handle& vB,
    const pointField& pts
) const
{
    if (f.size() < 3)
    {
        // Invalid face, fewer than three points

        return false;
    }
    else if (f.size() == 3)
    {
        // Triangle face, handle specially
    }
    else
    {
        // Polygonal face
    }

    scalar averageCellSize = averageAnyCellSize(vA, vB);

    return assessFace(f, averageCellSize, pts);
}


bool Foam::conformalVoronoiMesh::assessFace
(
    const face& f,
    scalar targetFaceSize,
    const pointField& pts
) const
{
    scalar smallFaceAreaCoeff = sqr(1e-5);
    scalar highAspectRatioFaceAreaCoeff = 0.1;
    scalar aspectRatioLimit = 2.0;
    scalar targetArea = sqr(targetFaceSize);

    const edgeList& eds = f.edges();

    scalar perimeter = 0.0;

    forAll(eds, i)
    {
        perimeter += eds[i].mag(pts);

        vector edVec = eds[i].vec(pts);
    };

    scalar area = f.mag(pts);

    scalar equivalentSqrPerimeter = 4.0*sqrt(area);

    scalar aspectRatio = perimeter/max(equivalentSqrPerimeter, VSMALL);

    bool keepFace = true;

    if (area < smallFaceAreaCoeff*targetArea)
    {
        keepFace = false;
    }

    if
    (
        aspectRatio > aspectRatioLimit
     && area < highAspectRatioFaceAreaCoeff*targetArea)
    {
        keepFace = false;
    }

    // if (!keepFace)
    // {
    //     Info<< nl << "Area " << area << nl
    //         << "targetFaceSize " << targetFaceSize << nl
    //         << "Area ratio "
    //         << area/max(sqr(targetFaceSize), VSMALL) << nl
    //         << "aspectRatio " << aspectRatio << nl
    //         << endl;

    //     forAll(f, i)
    //     {
    //         meshTools::writeOBJ(Info, pts[f[i]]);
    //     }

    //     Info<< nl;
    // }

    return keepFace;
}


Foam::label Foam::conformalVoronoiMesh::mergeCloseDualVertices
(
    const pointField& pts,
    Map<label>& dualPtIndexMap
)
{
    label nPtsMerged = 0;

    scalar closenessTolerance = 1e-4;

    for
    (
        Triangulation::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const int oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        label& c1I = c1->cellIndex();
        label& c2I = c2->cellIndex();

        if (dualPtIndexMap.found(c1I) || dualPtIndexMap.found(c2I))
        {
            // One of the points of this edge has already been
            // merged this sweep, leave for next sweep

            continue;
        }

        if(c1I != -1 && c2I != -1 && (c1I != c2I))
        {
            if
            (
                magSqr(pts[c1I] - pts[c2I])
              < sqr(averageAnyCellSize(fit)*closenessTolerance)
            )
            {
                dualPtIndexMap.insert(c1I, c1I);
                dualPtIndexMap.insert(c2I, c1I);

                nPtsMerged++;
            }
        }
    }

    return nPtsMerged;
}


Foam::label Foam::conformalVoronoiMesh::collapseFaces
(
    pointField& pts,
    Map<label>& dualPtIndexMap
)
{
    label nCollapsedFaces = 0;

    scalar smallEdgeLengthCoeff = 1e-3;
    scalar smallFaceAreaCoeff = sqr(smallEdgeLengthCoeff);
    scalar collapseToEdgeCoeff = 0.02;
    scalar longestEdgeLengthRatio = 0.35;

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_circulator ccStart = incident_cells(*eit);
        Cell_circulator cc = ccStart;

        do
        {
            if (dualPtIndexMap.found(cc->cellIndex()))
            {
                // One of the points of this face has already been
                // collapsed this sweep, leave for next sweep
                continue;
            }

        } while (++cc != ccStart);

        Cell_handle c = eit->first;
        Vertex_handle vA = c->vertex(eit->second);
        Vertex_handle vB = c->vertex(eit->third);

        if
        (
            vA->internalOrBoundaryPoint()
         || vB->internalOrBoundaryPoint()
        )
        {
            scalar targetFaceSize = averageAnyCellSize(vA, vB);
            scalar targetArea = sqr(targetFaceSize);

            face dualFace = buildDualFace(eit);

            if (dualFace.size() < 3)
            {
                // This face has been collapsed already
                continue;
            }

            scalar area = dualFace.mag(pts);

            if (area < smallFaceAreaCoeff*targetArea)
            {
                // Collapse the dual face

                // Determine if the face should be collapsed to a line or a
                // point

                const edgeList& eds = dualFace.edges();

                label longestEdgeI = -1;

                scalar longestEdgeLength = -SMALL;

                scalar perimeter = 0.0;

                forAll(eds, edI)
                {
                    scalar edgeLength = eds[edI].mag(pts);

                    perimeter += edgeLength;

                    if (edgeLength > longestEdgeLength)
                    {
                        longestEdgeI = edI;

                        longestEdgeLength = edgeLength;
                    }
                }

                if
                (
                    longestEdgeLength > collapseToEdgeCoeff*targetFaceSize
                 && longestEdgeLength/perimeter > longestEdgeLengthRatio
                )
                {
                    // Collapse to edge

                    // Start at either end of the longest edge and consume the
                    // rest of the points of the face

                    const edge& longestEd = eds[longestEdgeI];

                    label longestEdStartPtI = longestEd.start();
                    label longestEdEndPtI = longestEd.end();

                    label revEdI = longestEdgeI;
                    label fwdEdI = longestEdgeI;

                    point revPt = pts[longestEdStartPtI];
                    point fwdPt = pts[longestEdEndPtI];

                    dualPtIndexMap.insert(longestEdStartPtI, longestEdStartPtI);
                    dualPtIndexMap.insert(longestEdEndPtI, longestEdEndPtI);

                    // Circulate around the face

                    // Info<< nl << "# Before " << dualFace << nl
                    //     << "# area " << area << nl
                    //     << "# " << longestEdStartPtI << " " << longestEdEndPtI << nl
                    //     << endl;

                    // forAll(dualFace, fPtI)
                    // {
                    //     meshTools::writeOBJ(Info, pts[dualFace[fPtI]]);
                    // }

                    for (label fcI = 1; fcI <= label(eds.size()/2); fcI++)
                    {
                        revEdI = eds.rcIndex(revEdI);
                        fwdEdI = eds.fcIndex(fwdEdI);

                        const edge& revEd = eds[revEdI];
                        const edge& fwdEd = eds[fwdEdI];

                        if (fcI < label(eds.size()/2))
                        {
                            revPt += pts[revEd.start()];
                            fwdPt += pts[fwdEd.end()];

                            dualPtIndexMap.insert
                            (
                                revEd.start(),
                                longestEdStartPtI
                            );

                            dualPtIndexMap.insert
                            (
                                fwdEd.end(),
                                longestEdEndPtI
                            );
                        }
                        else
                        {
                            // Final circulation

                            if
                            (
                                eds.size() % 2 == 1
                                && revEd.start() == fwdEd.end()
                            )
                            {
                                // Odd number of edges, give final point to
                                // the edge direction that has the shorter
                                // final edge

                                if (fwdEd.mag(pts) < revEd.mag(pts))
                                {
                                    fwdPt += pts[fwdEd.end()];

                                    dualPtIndexMap.insert
                                    (
                                        fwdEd.end(),
                                        longestEdEndPtI
                                    );

                                    revPt /= fcI;
                                    fwdPt /= (fcI + 1);
                                }
                                else
                                {
                                    revPt += pts[revEd.start()];

                                    dualPtIndexMap.insert
                                    (
                                        revEd.start(),
                                        longestEdStartPtI
                                    );

                                    revPt /= (fcI + 1);
                                    fwdPt /= fcI;
                                }
                            }
                            else if
                            (
                                eds.size() % 2 == 0
                             && revEd.start() == fwdEd.start()
                             && revEd.end() == fwdEd.end()
                            )
                            {
                                // Even number of edges

                                revPt /= fcI;
                                fwdPt /= fcI;
                            }
                            else
                            {
                                FatalErrorIn("Foam::conformalVoronoiMesh::collapseFace")
                                    << "Face circulation failed for face "
                                    << dualFace << nl
                                    << exit(FatalError);
                            }
                        }
                    }

                    // Info<< "# dualPtIndexMap " << dualPtIndexMap << endl;

                    // Move the position of the accumulated points
                    pts[longestEdStartPtI] = revPt;
                    pts[longestEdEndPtI] = fwdPt;

                    // {
                    //     face checkDualFace = buildDualFace(eit);

                    //     Info<< "# After " << checkDualFace << endl;
                    // }

                    // Info<< "# Collapsed" << endl;

                    // meshTools::writeOBJ(Info, revPt);
                    // meshTools::writeOBJ(Info, fwdPt);

                    nCollapsedFaces++;
                }
                else if
                (
                    longestEdgeLength <= collapseToEdgeCoeff*targetFaceSize
                )
                {
                    // Collapse to point

                    // Cell_circulator ccStart = incident_cells(*eit);
                    // Cell_circulator cc1 = ccStart;
                    // Cell_circulator cc2 = cc1;

                    // // Advance the second circulator so that it always stays on the next
                    // // cell around the edge;
                    // cc2++;

                    // label nPts = 0;

                    // point resultantPt = vector::zero;

                    // label ccStartI = cc1->cellIndex();

                    // do
                    // {
                    //     label& cc1I = cc1->cellIndex();
                    //     label& cc2I = cc2->cellIndex();

                    //     if (cc1I < 0 || cc2I < 0)
                    //     {
                    //         FatalErrorIn("Foam::conformalVoronoiMesh::collapseFace")
                    //             << "Dual face uses circumcenter defined by a "
                    //             << "Delaunay tetrahedron with no internal "
                    //             << "or boundary points.  Defining Delaunay edge ends: "
                    //             << topoint(vA->point()) << " "
                    //             << topoint(vB->point()) << nl
                    //             << exit(FatalError);
                    //     }

                    //     if (cc1I != cc2I)
                    //     {
                    //         resultantPt += pts[cc1I];
                    //         nPts++;
                    //     }

                    //     cc1I = ccStartI;
                    //     cc2I = ccStartI;
                    //     cc1++;
                    //     cc2++;

                    // } while (cc1 != ccStart);

                    // resultantPt /= nPts;

                    // pts[ccStartI] = resultantPt;

                    point resultantPt = vector::zero;

                    label collapseToPtI = dualFace[0];

                    forAll(dualFace, fPtI)
                    {
                        label ptI = dualFace[fPtI];

                        resultantPt += pts[ptI];

                        dualPtIndexMap.insert(ptI, collapseToPtI);
                    }

                    resultantPt /= dualFace.size();

                    pts[collapseToPtI] = resultantPt;

                    nCollapsedFaces++;
                }
            }
        }
    }

    return nCollapsedFaces;
}


void Foam::conformalVoronoiMesh::reindexDualFace
(
    const Triangulation::Finite_edges_iterator& eit,
    const Map<label>& dualPtIndexMap
)
{
    Cell_circulator ccStart = incident_cells(*eit);
    Cell_circulator cc = ccStart;

    do
    {
        if (dualPtIndexMap.found(cc->cellIndex()))
        {
            cc->cellIndex() = dualPtIndexMap[cc->cellIndex()];
        }

        cc++;

    } while (cc != ccStart);
}


void Foam::conformalVoronoiMesh::reindexDualVertices
(
    const Map<label>& dualPtIndexMap
)
{
    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
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

    // Two stage sort:
    // 1) sort by owner

    Info<< nl
        << "    Sorting faces, owner and neighbour into upper triangular order"
        << endl;

    labelList oldToNew;

    sortedOrder(owner, oldToNew);

    oldToNew = invert(oldToNew.size(), oldToNew);

    inplaceReorder(oldToNew, faces);
    inplaceReorder(oldToNew, owner);
    inplaceReorder(oldToNew, neighbour);

    // 2) in each block of owners sort by neighbour

    // Reset map.  Elements that are not sorted will retain their -1
    // value, which will mean that they are ignored by inplaceReorder

    oldToNew = -1;

    label ownerBlockStart = 0;

    for (label o = 1; o < owner.size(); o++)
    {
        label blockLength = -1;

        if (owner[o] > owner[o-1])
        {
            blockLength = o - ownerBlockStart;
        }
        else if (o == owner.size() - 1)
        {
            // If the last element is not a jump in owner, then it
            // needs to trigger a sort of the last block, but with a
            // block length that is one element longer so that it
            // sorts itself.

            // If it is a jump in owner, then it will form a block of
            // length one, and so will not need sorted.

            blockLength = o - ownerBlockStart + 1;
        }

        if (blockLength >= 1)
        {
            labelList blockIndices =
                identity(blockLength) + ownerBlockStart;

            SubList<label> neighbourBlock
            (
                neighbour,
                blockLength,
                ownerBlockStart
            );

            sortedOrder(neighbourBlock, blockIndices);

            blockIndices = invert(blockIndices.size(), blockIndices);

            forAll(blockIndices, b)
            {
                oldToNew[ownerBlockStart + b] =
                blockIndices[b] + ownerBlockStart;
            }

            ownerBlockStart = o;
        }
    }

    // owner does not need re-sorted
    inplaceReorder(oldToNew, faces);
    inplaceReorder(oldToNew, neighbour);
}


void Foam::conformalVoronoiMesh::addPatches
(
    const label nInternalFaces,
    faceList& faces,
    labelList& owner,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    List<DynamicList<face> >& patchFaces,
    List<DynamicList<label> >& patchOwners,
    bool includeEmptyPatches
) const
{
    label nTotalPatches = patchNames.size();

    label nValidPatches = 0;

    PackedBoolList validPatch(nTotalPatches, false);

    wordList allPatchNames = patchNames;

    patchSizes.setSize(nTotalPatches);
    patchStarts.setSize(nTotalPatches);

    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        // Check if the patch has any faces.  Never create an empty
        // default patch.

        if
        (
            patchFaces[p].size()
         || (includeEmptyPatches && (p != nTotalPatches - 1))
        )
        {
            patchNames[nValidPatches] = allPatchNames[p];
            patchSizes[nValidPatches] = patchFaces[p].size();
            patchStarts[nValidPatches] = nInternalFaces + nBoundaryFaces;

            nBoundaryFaces += patchSizes[p];

            nValidPatches++;

            validPatch[p] = 1;
        }
        else
        {
            // Warn if a patch is empty and includeEmptyPatches is
            // false, unless it is the default patch.

            if (p != nTotalPatches - 1)
            {
                WarningIn("void addPatches")
                    << "Patch " << patchNames[p]
                    << " has no faces, not creating." << endl;
            }
        }
    }

    patchNames.setSize(nValidPatches);
    patchSizes.setSize(nValidPatches);
    patchStarts.setSize(nValidPatches);

    faces.setSize(nInternalFaces + nBoundaryFaces);
    owner.setSize(nInternalFaces + nBoundaryFaces);

    label faceI = nInternalFaces;

    forAll(patchFaces, p)
    {
        if (validPatch[p])
        {
            forAll(patchFaces[p], f)
            {
                faces[faceI] = patchFaces[p][f];
                owner[faceI] = patchOwners[p][f];

                faceI++;
            }
        }
    }
}


void Foam::conformalVoronoiMesh::removeUnusedPoints
(
    faceList& faces,
    pointField& pts
) const
{
    Info<< nl << "    Removing unused points after filtering" << endl;

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

    // Move all of the used faces to the start of the pointField and
    // truncate it

    forAll(ptUsed, ptUI)
    {
        if (ptUsed[ptUI] == true)
        {
            oldToNew[ptUI] = pointI++;
        }
    }

    inplaceReorder(oldToNew, pts);

    Info<< "        Removing "
        << pts.size() - pointI
        << " unused points" << endl;

    pts.setSize(pointI);

    // Renumber the faces to use the new point numbers

    forAll(faces, fI)
    {
        inplaceRenumber(oldToNew, faces[fI]);
    }
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
