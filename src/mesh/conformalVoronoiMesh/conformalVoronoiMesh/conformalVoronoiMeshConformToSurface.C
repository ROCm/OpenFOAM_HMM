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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

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
      % cvMeshControls().surfaceConformationRebuildFrequency() == 0
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
    vector& displacement,
    label callCount
) const
{
    callCount++;

    // Do not allow infinite recursion
    if (callCount > 7)
    {
        return;
    }

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

        limitDisplacement(vit, displacement, callCount);
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


// ************************************************************************* //
