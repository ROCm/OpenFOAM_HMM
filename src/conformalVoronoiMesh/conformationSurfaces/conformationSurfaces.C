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
    baffleSurfaces_(),
    bounds_()
{
    const dictionary& surfacesDict
    (
        surfaceConformationDict.subDict("geometryToConformTo")
    );

    Info<< nl << "Reading geometryToConformTo" << endl;

    surfaces_.setSize(surfacesDict.size());

    baffleSurfaces_.setSize(surfacesDict.size());

    features_.setSize(surfacesDict.size());

    label surfI = 0;

    forAllConstIter(dictionary, surfacesDict, iter)
    {
        word surfaceName = iter().keyword();

        surfaces_[surfI] = allGeometry_.findSurfaceID(surfaceName);

        if (surfaces_[surfI] < 0)
        {
            FatalErrorIn("Foam::conformationSurfaces::conformationSurfaces")
                << "No surface " << iter().keyword() << " found. "
                << "Valid geometry is " << nl << allGeometry_.names()
                << exit(FatalError);
        }

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

            features_.set
            (
                surfI,
                new featureEdgeMesh
                (
                    IOobject
                    (
                        feMeshName,
                        cvMesh_.time().findInstance
                        (
                            "featureEdgeMesh", feMeshName
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
            fileName feMeshName(surfaceName + "_noFeatures");

            features_.set
            (
                surfI,
                new featureEdgeMesh
                (
                    IOobject
                    (
                        feMeshName,
                        cvMesh_.time().constant(),
                        "featureEdgeMesh",
                        cvMesh_.time(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            );
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

    bounds_ = searchableSurfacesQueries::bounds(allGeometry_, surfaces_);

    writeFeatureObj("cvMesh");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformationSurfaces::~conformationSurfaces()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Field<bool> Foam::conformationSurfaces::inside
(
    const pointField& samplePts
) const
{
    return wellInside(samplePts, 0.0);
}


Foam::Field<bool> Foam::conformationSurfaces::outside
(
    const pointField& samplePts
) const
{
    return wellOutside(samplePts, 0.0);
}


Foam::Field<bool> Foam::conformationSurfaces::wellInside
(
    const pointField& samplePts,
    const scalar dist2
) const
{
    // Look at all surfaces at determine whether the locationInMesh point is
    // inside or outside each, to establish a signature for the domain to be
    // meshed.

    List<searchableSurface::volumeType> referenceVolumeTypes
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

            referenceVolumeTypes[s] = vTypes[0];
        }
    }

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

    Field<bool> insidePoints(samplePts.size(), true);

    //Check if the points are inside the surface by the given distance squared

    scalarField testDistSqr(insidePoints.size(), dist2);

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
            insidePoints[i] = false;

            continue;
        }

        forAll(surfaces_, s)
        {
            if (surfaceVolumeTests[s][i] != referenceVolumeTypes[s])
            {
                insidePoints[i] = false;

                break;
            }
        }
    }

    return insidePoints;
}


Foam::Field<bool> Foam::conformationSurfaces::wellOutside
(
    const pointField& samplePts,
    const scalar dist2
) const
{
    notImplemented("Field<bool> Foam::conformationSurfaces::wellOutside");

    return Field<bool>(samplePts.size(), true);
}


bool Foam::conformationSurfaces::findSurfaceAnyIntersection
(
    point start,
    point end
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
    point start,
    point end,
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
    point start,
    point end,
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
    featuresHit.setSize(featureEdgeMesh::nEdgeTypes);
    featuresHit = -1;
    edgeHits.setSize(featureEdgeMesh::nEdgeTypes);

    // Work arrays
    scalarField minDistSqr(featureEdgeMesh::nEdgeTypes, nearestDistSqr);
    List<pointIndexHit> hitInfo(featureEdgeMesh::nEdgeTypes);

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


void Foam::conformationSurfaces::writeFeatureObj
(
    const fileName& prefix
) const
{
    OFstream ftStr(prefix + "_allFeatures.obj");
    Pout<< nl << "Writing all features to " << ftStr.name() << endl;

    label verti = 0;

    forAll(features_, i)
    {
        const featureEdgeMesh& fEM(features_[i]);
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


// ************************************************************************* //
