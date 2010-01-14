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

    spanSqr_ = sqr(span_);

    timeChecks_ = true;

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

    featureEdgeExclusionDistanceCoeff_ = readScalar
    (
        surfDict.lookup("featureEdgeExclusionDistanceCoeff")
    );

    surfaceSearchDistanceCoeff_ = readScalar
    (
        surfDict.lookup("surfaceSearchDistanceCoeff")
    );

    maxSurfaceProtrusionCoeff_ = readScalar
    (
        surfDict.lookup("maxSurfaceProtrusionCoeff")
    );

    maxQuadAngle_ = readScalar(surfDict.lookup("maxQuadAngle"));

    surfaceConformationRebuildFrequency_ = max
    (
        1,
        readLabel(surfDict.lookup("surfaceConformationRebuildFrequency"))
    );

    // Motion control controls

    const dictionary& motionDict(cvMeshDict_.subDict("motionControl"));

    objOutput_ = Switch(motionDict.lookupOrDefault<Switch>("objOutput", false));

    alignmentSearchSpokes_ = readLabel
    (
        motionDict.lookup("alignmentSearchSpokes")
    );

    cosAlignmentAcceptanceAngle_ = cos
    (
        degToRad(readScalar(motionDict.lookup("alignmentAcceptanceAngle")))
    );

    // Point removal criteria

    const dictionary& insertionDict
    (
        motionDict.subDict("pointInsertionCriteria")
    );

    insertionDistCoeff_ = readScalar
    (
        insertionDict.lookup("cellCentreDistCoeff")
    );

    faceAreaRatioCoeff_ = readScalar
    (
        insertionDict.lookup("faceAreaRatioCoeff")
    );

    cosInsertionAcceptanceAngle_ = cos
    (
        degToRad(readScalar(insertionDict.lookup("acceptanceAngle")))
    );

    // Point removal criteria

    const dictionary& removalDict
    (
        motionDict.subDict("pointRemovalCriteria")
    );

    removalDistCoeff_ = readScalar
    (
        removalDict.lookup("cellCentreDistCoeff")
    );

    // polyMesh filtering controls

    const dictionary& filteringDict(cvMeshDict_.subDict("polyMeshFiltering"));

    filterSizeCoeff_ = readScalar
    (
        filteringDict.lookup("filterSizeCoeff")
    );

    mergeClosenessCoeff_ = readScalar
    (
        filteringDict.lookup("mergeClosenessCoeff")
    );

    continueFilteringOnBadInitialPolyMesh_ = Switch
    (
        filteringDict.lookupOrDefault<Switch>
        (
            "continueFilteringOnBadInitialPolyMesh",
            false
        )
    );

    surfaceStepFaceAngle_ = readScalar
    (
        filteringDict.lookup("surfaceStepFaceAngle")
    );

    edgeCollapseGuardFraction_ = readScalar
    (
        filteringDict.lookup("edgeCollapseGuardFraction")
    );

    maxCollapseFaceToPointSideLengthCoeff_ = readScalar
    (
        filteringDict.lookup("maxCollapseFaceToPointSideLengthCoeff")
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cvControls::~cvControls()
{}


// ************************************************************************* //
