/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "conformationSurfaces.H"
#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conformationSurfaces::conformationSurfaces
(
    const conformalVoronoiMesh& cvMesh,
    const searchableSurfaces& allGeometry,
    const dictionary& surfaceConformationDict
)
:
    cvMesh_(cvMesh),
    allGeometry_(allGeometry),
    features_(),
    locationInMesh_(surfaceConformationDict.lookup("locationInMesh")),
    surfaces_(),
    allGeometryToSurfaces_(),
    baffleSurfaces_(),
    patchNames_(0),
    patchOffsets_(),
    bounds_(),
    globalBounds_(),
    spanMag_(),
    spanMagSqr_(),
    processorDomains_(),
    referenceVolumeTypes_(0)
{
    const dictionary& surfacesDict
    (
        surfaceConformationDict.subDict("geometryToConformTo")
    );

    const dictionary& additionalFeaturesDict
    (
        surfaceConformationDict.subDict("additionalFeatures")
    );

    Info<< nl << "Reading geometryToConformTo" << endl;

    surfaces_.setSize(surfacesDict.size(), -1);

    allGeometryToSurfaces_.setSize(allGeometry_.size(), -1);

    baffleSurfaces_.setSize(surfacesDict.size(), false);

    // Features may be attached to host surfaces or independent
    features_.setSize(surfacesDict.size() + additionalFeaturesDict.size());

    label featureI = 0;

    patchOffsets_.setSize(surfacesDict.size(), -1);

    label surfI = 0;

    forAllConstIter(dictionary, surfacesDict, iter)
    {
        word surfaceName = iter().keyword();

        surfaces_[surfI] = allGeometry_.findSurfaceID(surfaceName);

        allGeometryToSurfaces_[surfaces_[surfI]] = surfI;

        if (surfaces_[surfI] < 0)
        {
            FatalErrorIn("Foam::conformationSurfaces::conformationSurfaces")
                << "No surface " << iter().keyword() << " found. "
                << "Valid geometry is " << nl << allGeometry_.names()
                << exit(FatalError);
        }

        Info<< nl << "    " << iter().keyword() << endl;

        patchOffsets_[surfI] = patchNames_.size();

        patchNames_.append(allGeometry.regionNames()[surfaces_[surfI]]);

        const dictionary& surfaceSubDict(surfacesDict.subDict(surfaceName));

        baffleSurfaces_[surfI] = Switch
        (
            surfaceSubDict.lookupOrDefault("baffleSurface", false)
        );

        word featureMethod = surfaceSubDict.lookupOrDefault
        (
            "featureMethod",
            word("none")
        );

        if (featureMethod == "featureEdgeMesh")
        {
            fileName feMeshName(surfaceSubDict.lookup("featureEdgeMesh"));

            Info<< "    features: " << feMeshName<< endl;

            features_.set
            (
                featureI++,
                new extendedFeatureEdgeMesh
                (
                    IOobject
                    (
                        feMeshName,
                        cvMesh_.time().findInstance
                        (
                            "featureEdgeMesh",
                            feMeshName
                        ),
                        "featureEdgeMesh",
                        cvMesh_.time(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );
        }
        else if (featureMethod == "extractFeatures")
        {
            notImplemented
            (
                "conformationSurfaces::conformationSurfaces, "
                "else if (featureMethod == \"extractFeatures\")"
            );
        }
        else if (featureMethod == "none")
        {
            // fileName feMeshName(surfaceName + "_noFeatures");

            // features_.set
            // (
            //     featureI++,
            //     new extendedFeatureEdgeMesh
            //     (
            //         IOobject
            //         (
            //             feMeshName,
            //             cvMesh_.time().constant(),
            //             "featureEdgeMesh",
            //             cvMesh_.time(),
            //             IOobject::NO_READ,
            //             IOobject::NO_WRITE,
            //             false
            //         )
            //     )
            // );
        }
        else
        {
            FatalErrorIn("Foam::conformationSurfaces::conformationSurfaces")
                << "No valid featureMethod found for surface " << surfaceName
                << nl << "Use \"featureEdgeMesh\" or \"extractFeatures\"."
                << exit(FatalError);
        }

        surfI++;
    }

    if (!additionalFeaturesDict.empty())
    {
        Info<< nl << "Reading additionalFeatures" << endl;
    }

    forAllConstIter(dictionary, additionalFeaturesDict, iter)
    {
        word featureName = iter().keyword();

        Info<< nl << "    " << iter().keyword() << endl;

        const dictionary& featureSubDict
        (
            additionalFeaturesDict.subDict(featureName)
        );

        word featureMethod = featureSubDict.lookupOrDefault
        (
            "featureMethod",
            word("none")
        );

        if (featureMethod == "featureEdgeMesh")
        {
            fileName feMeshName(featureSubDict.lookup("featureEdgeMesh"));

            Info<< "    features: " << feMeshName << endl;

            features_.set
            (
                featureI++,
                new extendedFeatureEdgeMesh
                (
                    IOobject
                    (
                        feMeshName,
                        cvMesh_.time().findInstance
                        (
                            "featureEdgeMesh",
                            feMeshName
                        ),
                        "featureEdgeMesh",
                        cvMesh_.time(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );
        }
    }

    // Remove unnecessary space from the features list
    features_.setSize(featureI);

    globalBounds_ = treeBoundBox
    (
        searchableSurfacesQueries::bounds(allGeometry_, surfaces_)
    );

    // Extend the global bounds to stop the bound box sitting on the surfaces
    // to be conformed to
    globalBounds_ = globalBounds_.extend(cvMesh_.rndGen(), 1e-4);

    if (Pstream::parRun())
    {
        label procLimit = 8;

        if (Pstream::nProcs() != procLimit)
        {
            FatalErrorIn
            (
                "Foam::conformationSurfaces::conformationSurfaces"
                "("
                    "const conformalVoronoiMesh& cvMesh, "
                    "const searchableSurfaces& allGeometry, "
                    "const dictionary& surfaceConformationDict"
                ")"
            )
                << "Hard coded to " << procLimit << " << procs"
                << exit(FatalError);
        }

        bounds_ = globalBounds_.subBbox(direction(Pstream::myProcNo()));

        processorDomains_.setSize(Pstream::nProcs());

        processorDomains_[Pstream::myProcNo()] = treeBoundBoxList(1, bounds_);

        Pstream::gatherList(processorDomains_);

        Pstream::scatterList(processorDomains_);
    }
    else
    {
        bounds_ = globalBounds_;

        processorDomains_.setSize(0);
    }

    if (Pstream::parRun() && cvMesh_.cvMeshControls().objOutput())
    {
        Info<< "global bounds " << globalBounds_ << endl;

        Pout<< "processor bounds " << bounds_ << endl;

        OFstream str
        (
            cvMesh_.time().path()
           /"proc_" + name(Pstream::myProcNo()) + "_bounds.obj"
        );

        Pout<< "Writing " << str.name() << endl;

        pointField bbPoints(bounds_.points());

        forAll(bbPoints, i)
        {
            meshTools::writeOBJ(str, bbPoints[i]);
        }

        forAll(treeBoundBox::faces, i)
        {
            const face& f = treeBoundBox::faces[i];

            str << "f"
                << ' ' << f[0] + 1
                << ' ' << f[1] + 1
                << ' ' << f[2] + 1
                << ' ' << f[3] + 1
                << nl;
        }

        Info<< "processor domains " << processorDomains_ << endl;
    }

    spanMag_ = bounds_.mag();

    spanMagSqr_ = sqr(spanMag_);

    // Look at all surfaces at determine whether the locationInMesh point is
    // inside or outside each, to establish a signature for the domain to be
    // meshed.

    referenceVolumeTypes_.setSize
    (
        surfaces_.size(),
        searchableSurface::UNKNOWN
    );

    forAll(surfaces_, s)
    {
        const searchableSurface& surface(allGeometry_[surfaces_[s]]);

        if (surface.hasVolumeType())
        {
            pointField pts(1, locationInMesh_);

            List<searchableSurface::volumeType> vTypes
            (
                pts.size(),
                searchableSurface::UNKNOWN
            );

            surface.getVolumeType(pts, vTypes);

            referenceVolumeTypes_[s] = vTypes[0];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformationSurfaces::~conformationSurfaces()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::conformationSurfaces::overlaps(const treeBoundBox& bb) const
{
    forAll(surfaces_, s)
    {
        if (allGeometry_[surfaces_[s]].overlaps(bb))
        {
            return true;
        }
    }

    return false;
}


bool Foam::conformationSurfaces::positionOnThisProc(const point& pt) const
{
    // This is likely to give problems when a point is on the boundary between
    // two processors.

    return bounds_.contains(pt);
}


Foam::Field<bool> Foam::conformationSurfaces::inside
(
    const pointField& samplePts
) const
{
    return wellInside(samplePts, scalarField(samplePts.size(), 0.0));
}


bool Foam::conformationSurfaces::inside
(
    const point& samplePt
) const
{
    return wellInside(pointField(1, samplePt), scalarField(1, 0))[0];
}


Foam::Field<bool> Foam::conformationSurfaces::outside
(
    const pointField& samplePts
) const
{
    return wellOutside(samplePts, scalarField(samplePts.size(), 0.0));
}


bool Foam::conformationSurfaces::outside
(
    const point& samplePt
) const
{
    return !inside(samplePt);
}


Foam::Field<bool> Foam::conformationSurfaces::wellInOutSide
(
    const pointField& samplePts,
    const scalarField& testDistSqr,
    bool testForInside
) const
{
    List<List<searchableSurface::volumeType> > surfaceVolumeTests
    (
        surfaces_.size(),
        List<searchableSurface::volumeType>
        (
            samplePts.size(),
            searchableSurface::UNKNOWN
        )
    );

    // Get lists for the volumeTypes for each sample wrt each surface
    forAll(surfaces_, s)
    {
        const searchableSurface& surface(allGeometry_[surfaces_[s]]);

        if (surface.hasVolumeType())
        {
            surface.getVolumeType(samplePts, surfaceVolumeTests[s]);
        }
    }

    // Compare the volumeType result for each point wrt to each surface with the
    // reference value and if the points are inside the surface by a given
    // distanceSquared

    // Assume that the point is wellInside until demonstrated otherwise.
    Field<bool> inOutSidePoint(samplePts.size(), testForInside);

    //Check if the points are inside the surface by the given distance squared

    labelList hitSurfaces;

    List<pointIndexHit> hitInfo;

    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        surfaces_,
        samplePts,
        testDistSqr,
        hitSurfaces,
        hitInfo
    );

    forAll(samplePts, i)
    {
        const pointIndexHit& pHit = hitInfo[i];

        if (pHit.hit())
        {
            // If the point is within range of the surface, then it can't be
            // well (in|out)side
            inOutSidePoint[i] = false;

            continue;
        }

        forAll(surfaces_, s)
        {
            // If one of the pattern tests is failed, then the point cannot be
            // inside, therefore, if this is a testForInside = true call, the
            // result is false.  If this is a testForInside = false call, then
            // the result is true.
            if (surfaceVolumeTests[s][i] != referenceVolumeTypes_[s])
            {
                inOutSidePoint[i] = !testForInside;

                break;
            }
        }
    }

    return inOutSidePoint;
}


Foam::Field<bool> Foam::conformationSurfaces::wellInside
(
    const pointField& samplePts,
    const scalarField& testDistSqr
) const
{
    return wellInOutSide(samplePts, testDistSqr, true);
}


bool Foam::conformationSurfaces::wellInside
(
    const point& samplePt,
    scalar testDistSqr
) const
{
    return wellInside(pointField(1, samplePt), scalarField(1, testDistSqr))[0];
}


Foam::Field<bool> Foam::conformationSurfaces::wellOutside
(
    const pointField& samplePts,
    const scalarField& testDistSqr
) const
{
    return wellInOutSide(samplePts, testDistSqr, false);
}


bool Foam::conformationSurfaces::wellOutside
(
    const point& samplePt,
    scalar testDistSqr
) const
{
    return wellOutside(pointField(1, samplePt), scalarField(1, testDistSqr))[0];
}


bool Foam::conformationSurfaces::findSurfaceAnyIntersection
(
    const point& start,
    const point& end
) const
{
    labelList hitSurfaces;
    List<pointIndexHit> hitInfo;

    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        surfaces_,
        pointField(1, start),
        pointField(1, end),
        hitSurfaces,
        hitInfo
    );

    return hitInfo[0].hit();
}


void Foam::conformationSurfaces::findSurfaceAnyIntersection
(
    const point& start,
    const point& end,
    pointIndexHit& surfHit,
    label& hitSurface
) const
{
    labelList hitSurfaces;
    List<pointIndexHit> hitInfo;

    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        surfaces_,
        pointField(1, start),
        pointField(1, end),
        hitSurfaces,
        hitInfo
    );

    surfHit = hitInfo[0];

    if (surfHit.hit())
    {
        // hitSurfaces has returned the index of the entry in surfaces_ that was
        // found, not the index of the surface in allGeometry_, translating this
        // to allGeometry_

        hitSurface = surfaces_[hitSurfaces[0]];
    }
}


void Foam::conformationSurfaces::findSurfaceNearestIntersection
(
    const point& start,
    const point& end,
    pointIndexHit& surfHit,
    label& hitSurface
) const
{
    labelList hitSurfacesStart;
    List<pointIndexHit> hitInfoStart;
    labelList hitSurfacesEnd;
    List<pointIndexHit> hitInfoEnd;

    searchableSurfacesQueries::findNearestIntersection
    (
        allGeometry_,
        surfaces_,
        pointField(1, start),
        pointField(1, end),
        hitSurfacesStart,
        hitInfoStart,
        hitSurfacesEnd,
        hitInfoEnd
    );

    surfHit = hitInfoStart[0];

    if (surfHit.hit())
    {
        // hitSurfaces has returned the index of the entry in surfaces_ that was
        // found, not the index of the surface in allGeometry_, translating this
        // to allGeometry_

        hitSurface = surfaces_[hitSurfacesStart[0]];
    }
}


void Foam::conformationSurfaces::findSurfaceNearest
(
    const point& sample,
    scalar nearestDistSqr,
    pointIndexHit& surfHit,
    label& hitSurface
) const
{
    labelList hitSurfaces;
    List<pointIndexHit> surfaceHits;

    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        surfaces_,
        pointField(1, sample),
        scalarField(1, nearestDistSqr),
        hitSurfaces,
        surfaceHits
    );

    surfHit = surfaceHits[0];

    if (surfHit.hit())
    {
        // hitSurfaces has returned the index of the entry in surfaces_ that was
        // found, not the index of the surface in allGeometry_, translating this
        // to allGeometry_

        hitSurface = surfaces_[hitSurfaces[0]];
    }
}


void Foam::conformationSurfaces::findSurfaceNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& surfaceHits,
    labelList& hitSurfaces
) const
{
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        surfaces_,
        samples,
        nearestDistSqr,
        hitSurfaces,
        surfaceHits
    );

    forAll(surfaceHits, i)
    {
        if (surfaceHits[i].hit())
        {
            // hitSurfaces has returned the index of the entry in surfaces_ that
            // was found, not the index of the surface in allGeometry_,
            // translating this to the surface in allGeometry_.

            hitSurfaces[i] = surfaces_[hitSurfaces[i]];
        }
    }
}


void Foam::conformationSurfaces::findEdgeNearest
(
    const point& sample,
    scalar nearestDistSqr,
    pointIndexHit& edgeHit,
    label& featureHit
) const
{
    pointField samples(1, sample);
    scalarField nearestDistsSqr(1, nearestDistSqr);

    List<pointIndexHit> edgeHits;
    labelList featuresHit;

    findEdgeNearest
    (
        samples,
        nearestDistsSqr,
        edgeHits,
        featuresHit
    );

    edgeHit = edgeHits[0];
    featureHit = featuresHit[0];
}


void Foam::conformationSurfaces::findEdgeNearest
(
    const pointField& samples,
    const scalarField& nearestDistsSqr,
    List<pointIndexHit>& edgeHits,
    labelList& featuresHit
) const
{
    // Initialise
    featuresHit.setSize(samples.size());
    featuresHit = -1;
    edgeHits.setSize(samples.size());

    // Work arrays
    scalarField minDistSqr(nearestDistsSqr);
    List<pointIndexHit> hitInfo(samples.size());

    forAll(features_, testI)
    {
        features_[testI].nearestFeatureEdge
        (
            samples,
            minDistSqr,
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, pointI)
        {
            if (hitInfo[pointI].hit())
            {
                minDistSqr[pointI] = magSqr
                (
                    hitInfo[pointI].hitPoint()
                  - samples[pointI]
                );
                edgeHits[pointI] = hitInfo[pointI];
                featuresHit[pointI] = testI;
            }
        }
    }
}


void Foam::conformationSurfaces::findEdgeNearestByType
(
    const point& sample,
    scalar nearestDistSqr,
    List<pointIndexHit>& edgeHits,
    List<label>& featuresHit
) const
{
    // Initialise
    featuresHit.setSize(extendedFeatureEdgeMesh::nEdgeTypes);
    featuresHit = -1;
    edgeHits.setSize(extendedFeatureEdgeMesh::nEdgeTypes);

    // Work arrays
    scalarField minDistSqr(extendedFeatureEdgeMesh::nEdgeTypes, nearestDistSqr);
    List<pointIndexHit> hitInfo(extendedFeatureEdgeMesh::nEdgeTypes);

    forAll(features_, testI)
    {
        features_[testI].nearestFeatureEdgeByType
        (
            sample,
            minDistSqr,
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, typeI)
        {
            if (hitInfo[typeI].hit())
            {
                minDistSqr[typeI] = magSqr(hitInfo[typeI].hitPoint() - sample);
                edgeHits[typeI] = hitInfo[typeI];
                featuresHit[typeI] = testI;
            }
        }
    }
}


void Foam::conformationSurfaces::writeFeatureObj(const fileName& prefix) const
{
    OFstream ftStr(cvMesh_.time().path()/prefix + "_allFeatures.obj");

    Pout<< nl << "Writing all features to " << ftStr.name() << endl;

    label verti = 0;

    forAll(features_, i)
    {
        const extendedFeatureEdgeMesh& fEM(features_[i]);
        const pointField pts(fEM.points());
        const edgeList eds(fEM.edges());

        ftStr << "g " << fEM.name() << endl;

        forAll(eds, j)
        {
            const edge& e = eds[j];

            meshTools::writeOBJ(ftStr, pts[e[0]]); verti++;
            meshTools::writeOBJ(ftStr, pts[e[1]]); verti++;
            ftStr << "l " << verti-1 << ' ' << verti << endl;
        }
    }
}


Foam::label Foam::conformationSurfaces::findPatch
(
    const point& ptA,
    const point& ptB
) const
{
    pointIndexHit surfHit;
    label hitSurface;

    findSurfaceAnyIntersection(ptA, ptB, surfHit, hitSurface);

    if (!surfHit.hit())
    {
        return -1;
    }

    labelList surfLocalRegion;

    allGeometry_[hitSurface].getRegion
    (
        List<pointIndexHit>(1, surfHit),
        surfLocalRegion
    );

    return
        surfLocalRegion[0] + patchOffsets_[allGeometryToSurfaces_[hitSurface]];
}


Foam::label Foam::conformationSurfaces::findPatch(const point& pt) const
{
    pointIndexHit surfHit;
    label hitSurface;

    findSurfaceNearest
    (
        pt,
        spanMagSqr_,
        surfHit,
        hitSurface
    );

    if (!surfHit.hit())
    {
        return -1;
    }

    labelList surfLocalRegion;

    allGeometry_[hitSurface].getRegion
    (
        List<pointIndexHit>(1, surfHit),
        surfLocalRegion
    );

    return
        surfLocalRegion[0] + patchOffsets_[allGeometryToSurfaces_[hitSurface]];
}


// ************************************************************************* //
