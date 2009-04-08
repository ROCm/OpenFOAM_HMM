/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    features_
    (
        IOobject
        (
            "features",
            cvMesh_.time().constant(),
            "featureEdgeMesh",
            cvMesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointField(0),
        edgeList(0)
    ),
    locationInMesh_(surfaceConformationDict.lookup("locationInMesh")),
    surfaces_(0),
    bounds_()
{
    const dictionary& surfacesDict
    (
        surfaceConformationDict.subDict("geometryToConformTo")
    );

    surfaces_.setSize(surfacesDict.size());

    label surfI = 0;

    forAllConstIter(dictionary, surfacesDict, iter)
    {
        surfaces_[surfI] = allGeometry_.findSurfaceID(iter().keyword());

        if (surfaces_[surfI] < 0)
        {
            FatalErrorIn("Foam::conformationSurfaces::conformationSurfaces")
                << "No surface " << iter().keyword() << " found. "
                << "Valid geometry is " << nl << allGeometry_.names()
                << exit(FatalError);
        }

        surfI++;
    }

    bounds_ = searchableSurfacesQueries::bounds(allGeometry_, surfaces_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformationSurfaces::~conformationSurfaces()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



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



// ************************************************************************* //
