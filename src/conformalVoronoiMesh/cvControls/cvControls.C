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

#include "cvControls.H"
#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cvControls::cvControls
(
    const conformalVoronoiMesh& cvMesh,
    const IOdictionary& cvMeshDict
)
:
    cvMesh_(cvMesh),
    cvMeshDict_(cvMeshDict)
{
    // General parameters

    const boundBox& bb = cvMesh_.geometryToConformTo().bounds();

    span_ =
        max(mag(bb.max().x()), mag(bb.min().x()))
      + max(mag(bb.max().y()), mag(bb.min().y()))
      + max(mag(bb.max().z()), mag(bb.min().z()));

    // Surface conformation controls

    const dictionary& surfDict(cvMeshDict_.subDict("surfaceConformation"));

    pointPairDistanceCoeff_ = readScalar
    (
        surfDict.lookup("pointPairDistanceCoeff")
    );

    mixedFeaturePointPPDistanceCoeff_ = readScalar
    (
        surfDict.lookup("mixedFeaturePointPPDistanceCoeff")
    );

    featurePointExclusionDistanceCoeff_ = readScalar
    (
        surfDict.lookup("featurePointExclusionDistanceCoeff")
    );

    surfaceSearchDistanceCoeff_ = readScalar
    (
        surfDict.lookup("surfaceSearchDistanceCoeff")
    );

    maxSurfaceProtrusionCoeff_ = readScalar
    (
        surfDict.lookup("maxSurfaceProtrusionCoeff")
    );

    maxQuadAngle_= readScalar(surfDict.lookup("maxQuadAngle"));

    // Motion control controls

    const dictionary& motionDict(cvMeshDict_.subDict("motionControl"));

    const dictionary& insertionDict
    (
        motionDict.subDict("pointInsertionCriteria")
    );

    cellCentreInsertionDistCoeff_ = readScalar
    (
        insertionDict.lookup("cellCentreDistCoeff")
    );

    faceAreaRatioCoeff_ = readScalar
    (
        insertionDict.lookup("faceAreaRatioCoeff")
    );

    alignmentAcceptanceAngle_ = readScalar
    (
        insertionDict.lookup("alignmentAcceptanceAngle")
    );

    const dictionary& removalDict
    (
        motionDict.subDict("pointRemovalCriteria")
    );

    cellCentreRemovalDistCoeff_ = readScalar
    (
        removalDict.lookup("cellCentreDistCoeff")
    );

    // polyMesh filtering controls

    const dictionary& filteringDict(cvMeshDict_.subDict("polyMeshFiltering"));

    minimumEdgeLengthCoeff_ = readScalar
    (
        filteringDict.lookup("minimumEdgeLengthCoeff")
    );

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cvControls::~cvControls()
{}


// ************************************************************************* //
