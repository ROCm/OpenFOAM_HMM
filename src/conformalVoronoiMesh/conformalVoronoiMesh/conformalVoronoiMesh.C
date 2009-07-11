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
#include "uint.H"
#include "ulong.H"

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
    timeCheck();

    createFeaturePoints();
    timeCheck();

    insertInitialPoints();
    timeCheck();

    conformToSurface();
    timeCheck();

    writePoints("allPoints.obj", false);
    timeCheck();

    writeMesh();
    timeCheck();

    writeTargetCellSize();
    timeCheck();
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
            Foam::cos(i*mathematicalConstant::twoPi/s),
            Foam::sin(i*mathematicalConstant::twoPi/s),
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

    allGeometry_[hitSurface].getNormal
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
            << "Parallel normals detected in spoke search."
            << nl << exit(FatalError);
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

    if (fName != fileName::null)
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

    if (fName != fileName::null)
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
        180*(mathematicalConstant::pi + acos(mag(nA & nB)))
       /mathematicalConstant::pi;

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

    insertPoint(point(-bigSpan, -bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan, -bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan,  bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan,  bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan, -bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan, -bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan,  bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan,  bigSpan , bigSpan), Vb::FAR_POINT);

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

    writePoints("featureVertices.obj", false);
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

    label nVert = startOfInternalPoints_;

    Info<< nl << "Inserting initial points" << endl;

    std::vector<Point> initialPoints = initialPointsMethod_->initialPoints();

    Info<< "    " << initialPoints.size() << " points to insert..." << endl;

    // using the range insert (faster than inserting points one by one)
    insert(initialPoints.begin(), initialPoints.end());

    Info<< "    " << number_of_vertices() - startOfInternalPoints_
        << " points inserted" << endl;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->uninitialised())
        {
            vit->index() = nVert++;
        }
    }

    writePoints("initialPoints.obj", true);
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


void Foam::conformalVoronoiMesh::addSurfaceAndEdgeHits
(
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

    bool keepSurfacePoint = true;

    if (nearFeaturePt(surfHit.hitPoint()))
    {
        keepSurfacePoint = false;
    }

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
            if
            (
                magSqr(edHit.hitPoint() - surfHit.hitPoint())
              < surfacePtReplaceDistCoeffSqr*targetCellSizeSqr
            )
            {
                // If the point is within a given distance of a feature edge,
                // give control to edge control points instead, this will
                // prevent "pits" forming.

                keepSurfacePoint = false;
            }

            if (!nearFeaturePt(edHit.hitPoint()))
            {
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
    // ~~~~~~~~~~~ removing short edges by indexing dual vertices ~~~~~~~~~~~~~~

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        cit->cellIndex() = -1;
    }

    points.setSize(number_of_cells());

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Looking up minEdgeLenSqr with a dummy point, in future will be available
    // as a local value to be looked up in-place.

    scalar minEdgeLenSqr = sqr(minimumEdgeLength(vector::zero));

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label dualVerti = 0;

    // Scanning by number of short (dual) edges (nSE) attached to the
    // circumcentre of each Delaunay tet.  A Delaunay tet may only have four
    // dual edges emanating from its circumcentre, assigning positions and
    // indices to those with 4 short edges attached first, then >= 3, then >= 2
    // etc.
    for (label nSE = 4; nSE >= 0; nSE--)
    {
        Info<< nl << "Scanning for dual vertices with >= "
            << nSE
            << " short edges attached." << endl;

        for
        (
            Triangulation::Finite_cells_iterator cit = finite_cells_begin();
            cit != finite_cells_end();
            ++cit
        )
        {
            // If the Delaunay tet has an index already then it has either
            // evaluated itself and taken action or has had its index dictated
            // by a neighbouring tet with more short edges attached.

            if (cit->cellIndex() == -1)
            {
                point dualVertex = topoint(dual(cit));

                label shortEdges = 0;

                List<bool> edgeIsShort(4, false);

                List<bool> neighbourAlreadyIndexed(4, false);

                // Loop over the four facets of the Delaunay tet
                for (label f = 0; f < 4; f++)
                {
                    // Check that at least one of the vertices of the facet is
                    // an internal or boundary point
                    if
                    (
                        cit->vertex(vertex_triple_index(f, 0))->
                        internalOrBoundaryPoint()
                        || cit->vertex(vertex_triple_index(f, 1))->
                        internalOrBoundaryPoint()
                        || cit->vertex(vertex_triple_index(f, 2))->
                        internalOrBoundaryPoint()
                    )
                    {
                        point neighDualVertex;

                        label cNI = cit->neighbor(f)->cellIndex();

                        if (cNI == -1)
                        {
                            neighDualVertex = topoint(dual(cit->neighbor(f)));
                        }
                        else
                        {
                            neighDualVertex = points[cNI];
                        }

                        if
                        (
                            magSqr(dualVertex - neighDualVertex) < minEdgeLenSqr
                        )
                        {
                            edgeIsShort[f] = true;

                            if (cNI > -1)
                            {
                                neighbourAlreadyIndexed[f] = true;
                            }

                            shortEdges++;
                        }
                    }
                }

                if (nSE == 0 && shortEdges == 0)
                {
                    // Final iteration and no short edges are found, index
                    // remaining dual vertices.

                    if
                    (
                        cit->vertex(0)->internalOrBoundaryPoint()
                     || cit->vertex(1)->internalOrBoundaryPoint()
                     || cit->vertex(2)->internalOrBoundaryPoint()
                     || cit->vertex(3)->internalOrBoundaryPoint()
                    )
                    {
                        cit->cellIndex() = dualVerti;
                        points[dualVerti] = dualVertex;
                        dualVerti++;
                    }
                }
                else if
                (
                    shortEdges >= nSE
                )
                {
                    // Info<< neighbourAlreadyIndexed << ' '
                    //     << edgeIsShort << endl;

                    label numUnindexedNeighbours = 1;

                    for (label f = 0; f < 4; f++)
                    {
                        if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                        {
                            dualVertex += topoint(dual(cit->neighbor(f)));

                            numUnindexedNeighbours++;
                        }
                    }

                    dualVertex /= numUnindexedNeighbours;

                    label nearestExistingIndex = -1;

                    point nearestIndexedNeighbourPos = vector::zero;

                    scalar minDistSqrToNearestIndexedNeighbour = VGREAT;

                    for (label f = 0; f < 4; f++)
                    {
                        if (edgeIsShort[f] && neighbourAlreadyIndexed[f])
                        {
                            label cNI = cit->neighbor(f)->cellIndex();

                            point indexedNeighbourPos = points[cNI];

                            if
                            (
                                magSqr(indexedNeighbourPos - dualVertex)
                              < minDistSqrToNearestIndexedNeighbour
                            )
                            {
                                nearestExistingIndex = cNI;

                                nearestIndexedNeighbourPos =
                                indexedNeighbourPos;

                                minDistSqrToNearestIndexedNeighbour =
                                magSqr(indexedNeighbourPos - dualVertex);
                            }
                        }
                    }

                    if
                    (
                        nearestExistingIndex > -1
                     && minDistSqrToNearestIndexedNeighbour < minEdgeLenSqr
                    )
                    {
                        points[nearestExistingIndex] =
                        0.5*(dualVertex + nearestIndexedNeighbourPos);

                        for (label f = 0; f < 4; f++)
                        {
                            if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                            {
                                cit->neighbor(f)->cellIndex() =
                                nearestExistingIndex;
                            }
                        }

                        cit->cellIndex() = nearestExistingIndex;
                    }
                    else
                    {
                        for (label f = 0; f < 4; f++)
                        {
                            if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                            {
                                cit->neighbor(f)->cellIndex() = dualVerti;
                            }
                        }

                        cit->cellIndex() = dualVerti;

                        points[dualVerti] = dualVertex;

                        dualVerti++;
                    }
                }
            }
        }
    }

    points.setSize(dualVerti);

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
            vit->type() = Vb::INTERNAL_POINT;
            vit->index() = dualCelli;
            dualCelli++;
        }
        else
        {
            vit->type() = Vb::FAR_POINT;
            vit->index() = -1;
        }
    }

    // ~~~~~~~~~~~~ dual face and owner neighbour construction ~~~~~~~~~~~~~~~~~

    //label nPatches = qSurf_.patches().size() + 1;

    //label defaultPatchIndex = qSurf_.patches().size();

    label nPatches = 1;

    label defaultPatchIndex = 0;

    patchNames.setSize(nPatches);

    //const geometricSurfacePatchList& surfacePatches = qSurf_.patches();

    // forAll(surfacePatches, sP)
    // {
    //     patchNames[sP] = surfacePatches[sP].name();
    // }

    patchNames[defaultPatchIndex] = "cvMesh_defaultPatch";

    patchSizes.setSize(nPatches);

    patchStarts.setSize(nPatches);

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_edges());

    owner.setSize(number_of_edges());

    neighbour.setSize(number_of_edges());

    label dualFacei = 0;

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
                    FatalErrorIn("Foam::conformalVoronoiMesh::calcDualMesh")
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

            if (verticesOnFace.size() >= 3)
            {
                face newDualFace(verticesOnFace);

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

                    //pointIndexHit pHit = qSurf_.tree().findLineAny(ptA, ptB);

                    //label patchIndex = qSurf_[pHit.index()].region();

                    label patchIndex = defaultPatchIndex;

                    if (patchIndex == -1)
                    {
                        patchIndex = defaultPatchIndex;

                        WarningIn("Foam::conformalVoronoiMesh::calcDualMesh")
                            << "Dual face found that is not on a surface "
                            << "patch. Adding to "
                            << patchNames[defaultPatchIndex]
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

                    faces[dualFacei] = newDualFace;

                    owner[dualFacei] = dcOwn;

                    neighbour[dualFacei] = dcNei;

                    dualFacei++;
                }
            }
            // else
            // {
            //     Info<< verticesOnFace.size()
            //         << " size face not created." << endl;
            // }
        }
    }

    label nInternalFaces = dualFacei;

    faces.setSize(nInternalFaces);

    owner.setSize(nInternalFaces);

    neighbour.setSize(nInternalFaces);

    // ~~~~~~~~ sort owner, reordinging neighbour and faces to match ~~~~~~~~~~~
    // two stage sort for upper triangular order:  sort by owner first, then for
    // each block of owners sort by neighbour

    labelList sortingIndices;

    // Stage 1

    {
        SortableList<label> sortedOwner(owner);

        sortingIndices = sortedOwner.indices();
    }

    {
        labelList copyOwner(owner.size());

        forAll(sortingIndices, sI)
        {
            copyOwner[sI] = owner[sortingIndices[sI]];
        }

        owner = copyOwner;
    }

    {
        labelList copyNeighbour(neighbour.size());

        forAll(sortingIndices, sI)
        {
            copyNeighbour[sI] = neighbour[sortingIndices[sI]];
        }

        neighbour = copyNeighbour;
    }

    {
        faceList copyFaces(faces.size());

        forAll(sortingIndices, sI)
        {
            copyFaces[sI] = faces[sortingIndices[sI]];
        }

        faces = copyFaces;
    }

    // Stage 2

    sortingIndices = -1;

    DynamicList<label> ownerCellJumps;

    // Force first owner entry to be a jump
    ownerCellJumps.append(0);

    for (label o = 1; o < owner.size(); o++)
    {
        if (owner[o] > owner[o-1])
        {
            ownerCellJumps.append(o);
        }
    }

    forAll(ownerCellJumps, oCJ)
    {
        label start = ownerCellJumps[oCJ];

        label length;

        if (oCJ == ownerCellJumps.size() - 1)
        {
            length = owner.size() - start;
        }
        else
        {
            length = ownerCellJumps[oCJ + 1] - start;
        }

        SubList<label> neighbourBlock(neighbour, length, start);

        SortableList<label> sortedNeighbourBlock(neighbourBlock);

        forAll(sortedNeighbourBlock, sNB)
        {
            sortingIndices[start + sNB] =
            sortedNeighbourBlock.indices()[sNB] + start;
        }
    }

    // Perform sort

    {
        labelList copyOwner(owner.size());

        forAll(sortingIndices, sI)
        {
            copyOwner[sI] = owner[sortingIndices[sI]];
        }

        owner = copyOwner;
    }

    {
        labelList copyNeighbour(neighbour.size());

        forAll(sortingIndices, sI)
        {
            copyNeighbour[sI] = neighbour[sortingIndices[sI]];
        }

        neighbour = copyNeighbour;
    }

    {
        faceList copyFaces(faces.size());

        forAll(sortingIndices, sI)
        {
            copyFaces[sI] = faces[sortingIndices[sI]];
        }

        faces = copyFaces;
    }

    // ~~~~~~~~ add patch information ~~~~~~~~~~~

    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        patchSizes[p] = patchFaces[p].size();

        patchStarts[p] = nInternalFaces + nBoundaryFaces;

        nBoundaryFaces += patchSizes[p];
    }

    faces.setSize(nInternalFaces + nBoundaryFaces);

    owner.setSize(nInternalFaces + nBoundaryFaces);

    forAll(patchFaces, p)
    {
        forAll(patchFaces[p], f)
        {
            faces[dualFacei] = patchFaces[p][f];

            owner[dualFacei] = patchOwners[p][f];

            dualFacei++;
        }
    }
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::conformToSurface()
{
    Info<< nl << "Conforming to surfaces" << endl;

    startOfSurfacePointPairs_ = number_of_vertices();

    // Initialise containers to store the edge conformation locations
    DynamicList<point> newEdgeLocations;

    pointField existingEdgeLocations(0);

    autoPtr<indexedOctree<treeDataPoint> > edgeLocationTree;

    // Initialise the edgeLocationTree
    buildEdgeLocationTree(edgeLocationTree, existingEdgeLocations);

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

        Info<< nl <<"    Initial conformation " << nl
            << "    number_of_vertices " << number_of_vertices() << nl
            << "    surfaceHits.size() " << surfaceHits.size() << nl
            << "    featureEdgeHits.size() " << featureEdgeHits.size()
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
    }

    label iterationNo = 0;

    label maxIterations = 10;
    Info << "    MAX ITERATIONS HARD CODED TO "<< maxIterations << endl;

    // Set totalHits to a positive value to enter the while loop on the first
    // iteration
    label totalHits = 1;

    while (totalHits > 0 && iterationNo < maxIterations)
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

        Info<< nl <<"    iterationNo " << iterationNo << nl
            << "    number_of_vertices " << number_of_vertices() << nl
            << "    surfaceHits.size() " << surfaceHits.size() << nl
            << "    featureEdgeHits.size() " << featureEdgeHits.size()
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
    }
}


void Foam::conformalVoronoiMesh::move()
{
    timeCheck();

    scalar relaxation = relaxationModel_->relaxation();

    Info<< nl << "    Relaxation = " << relaxation << endl;

    pointField dualVertices(number_of_cells());

    pointField displacementAccumulator(startOfSurfacePointPairs_, vector::zero);

    List<bool> pointToBeRetained(startOfSurfacePointPairs_, true);

    DynamicList<point> newPointsToInsert;

    label dualVerti = 0;

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
            cit->cellIndex() = dualVerti;

            dualVertices[dualVerti] = topoint(dual(cit));

            dualVerti++;
        }
    }

    dualVertices.setSize(dualVerti);

    timeCheck();

    Info<< nl << "    Calculating target cell alignment and size" << endl;

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

            vit->alignment() = requiredAlignment(pt);

            vit->targetCellSize() = targetCellSize(pt);
        }
    }

    timeCheck();

    Info<< nl << "    Looping over all dual faces" << endl;

    vectorField cartesianDirections(3);

    cartesianDirections[0] = vector(0,0,1);
    cartesianDirections[1] = vector(0,1,0);
    cartesianDirections[2] = vector(1,0,0);

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
            Cell_circulator ccStart = incident_cells(*eit);
            Cell_circulator cc = ccStart;

            DynamicList<label> verticesOnFace;

            do
            {
                if (!is_infinite(cc))
                {
                    if (cc->cellIndex() < 0)
                    {
                        FatalErrorIn("conformalVoronoiMesh::move")
                            << "Dual face uses circumcenter defined by a "
                            << " Delaunay tetrahedron with no internal "
                            << "or boundary points."
                            << exit(FatalError);
                    }

                    verticesOnFace.append(cc->cellIndex());
                }
            } while (++cc != ccStart);

            verticesOnFace.shrink();

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

                scalar targetCellSize =
                0.5*(vA->targetCellSize() + vB->targetCellSize());

                scalar targetFaceArea = sqr(targetCellSize);

                if
                (
                    alignmentDotProd
                    > cvMeshControls().cosAlignmentAcceptanceAngle()
                )
                {
                    alignmentDir *= 0.5*targetCellSize;

                    vector delta = alignmentDir - 0.5*rAB;
                }
            }
        }
    }
}


// ************************************************************************* //
