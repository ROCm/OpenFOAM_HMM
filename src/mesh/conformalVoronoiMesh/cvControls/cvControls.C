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

    // Controls for coarse surface conformation

    const dictionary& coarseDict
    (
        surfDict.subDict("coarseConformationControls")
    );

    const dictionary& coarseInitialDict
    (
        coarseDict.subDict("initial")
    );

    const dictionary& coarseIterationDict
    (
        coarseDict.subDict("iteration")
    );

    edgeSearchDistCoeffSqr_coarse_initial_ = sqr
    (
        readScalar
        (
            coarseInitialDict.lookup("edgeSearchDistCoeff")
        )
    );

    edgeSearchDistCoeffSqr_coarse_iteration_ = sqr
    (
        readScalar
        (
            coarseIterationDict.lookup("edgeSearchDistCoeff")
        )
    );

    surfacePtReplaceDistCoeffSqr_coarse_initial_ = sqr
    (
        readScalar
        (
            coarseInitialDict.lookup("surfacePtReplaceDistCoeff")
        )
    );

    surfacePtReplaceDistCoeffSqr_coarse_iteration_ = sqr
    (
        readScalar
        (
            coarseIterationDict.lookup("surfacePtReplaceDistCoeff")
        )
    );

    maxConformationIterations_coarse_ = readLabel
    (
        coarseDict.lookup("maxIterations")
    );

    iterationToIntialHitRatioLimit_coarse_ = readScalar
    (
        coarseDict.lookup("iterationToIntialHitRatioLimit")
    );

    // Controls for fine surface conformation

    const dictionary& fineDict
    (
        surfDict.subDict("fineConformationControls")
    );

    const dictionary& fineInitialDict
    (
        fineDict.subDict("initial")
    );

    const dictionary& fineIterationDict
    (
        fineDict.subDict("iteration")
    );

    edgeSearchDistCoeffSqr_fine_initial_ = sqr
    (
        readScalar
        (
            fineInitialDict.lookup("edgeSearchDistCoeff")
        )
    );

    edgeSearchDistCoeffSqr_fine_iteration_ = sqr
    (
        readScalar
        (
            fineIterationDict.lookup("edgeSearchDistCoeff")
        )
    );

    surfacePtReplaceDistCoeffSqr_fine_initial_ = sqr
    (
        readScalar
        (
            fineInitialDict.lookup("surfacePtReplaceDistCoeff")
        )
    );

    surfacePtReplaceDistCoeffSqr_fine_iteration_ = sqr
    (
        readScalar
        (
            fineIterationDict.lookup("surfacePtReplaceDistCoeff")
        )
    );

    maxConformationIterations_fine_ = readLabel
    (
        fineDict.lookup("maxIterations")
    );

    iterationToIntialHitRatioLimit_fine_ = readScalar
    (
        fineDict.lookup("iterationToIntialHitRatioLimit")
    );

    // Motion control controls

    const dictionary& motionDict(cvMeshDict_.subDict("motionControl"));

    objOutput_ = Switch(motionDict.lookupOrDefault<Switch>("objOutput", false));

    timeChecks_ = Switch
    (
        motionDict.lookupOrDefault<Switch>("timeChecks", false)
    );

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

    filterErrorReductionCoeff_ = readScalar
    (
        filteringDict.lookup("filterErrorReductionCoeff")
    );

    filterCountSkipThreshold_ = readLabel
    (
        filteringDict.lookup("filterCountSkipThreshold")
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::cvControls::edgeSearchDistCoeffSqrInitial
(
    int reconfMode
) const
{
    if (reconfMode == conformalVoronoiMesh::rmCoarse)
    {
        return edgeSearchDistCoeffSqr_coarse_initial_;
    }

    else if (reconfMode == conformalVoronoiMesh::rmFine)
    {
        return edgeSearchDistCoeffSqr_fine_initial_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::cvControls::edgeSearchDistCoeffSqrInitial"
            "("
                "int reconfMode"
            ") const"
        )   << "Unknown reconformationMode " << reconfMode
            << exit(FatalError);
    }

    return 0;
}


Foam::scalar Foam::cvControls::edgeSearchDistCoeffSqrIteration
(
    int reconfMode
) const
{
    if (reconfMode == conformalVoronoiMesh::rmCoarse)
    {
        return edgeSearchDistCoeffSqr_coarse_iteration_;
    }

    else if (reconfMode == conformalVoronoiMesh::rmFine)
    {
        return edgeSearchDistCoeffSqr_fine_iteration_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::cvControls::edgeSearchDistCoeffSqrIteration"
            "("
                "int reconfMode"
            ") const"
        )   << "Unknown reconformationMode " << reconfMode
            << exit(FatalError);
    }

    return 0;
}


Foam::scalar Foam::cvControls::surfacePtReplaceDistCoeffSqrInitial
(
    int reconfMode
) const
{
    if (reconfMode == conformalVoronoiMesh::rmCoarse)
    {
        return surfacePtReplaceDistCoeffSqr_coarse_initial_;
    }

    else if (reconfMode == conformalVoronoiMesh::rmFine)
    {
        return surfacePtReplaceDistCoeffSqr_fine_initial_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::cvControls::surfacePtReplaceDistCoeffSqrInitial"
            "("
                "int reconfMode"
            ") const"
        )   << "Unknown reconformationMode " << reconfMode
            << exit(FatalError);
    }

    return 0;
}


Foam::scalar Foam::cvControls::surfacePtReplaceDistCoeffSqrIteration
(
    int reconfMode
) const
{
    if (reconfMode == conformalVoronoiMesh::rmCoarse)
    {
        return surfacePtReplaceDistCoeffSqr_coarse_iteration_;
    }

    else if (reconfMode == conformalVoronoiMesh::rmFine)
    {
        return surfacePtReplaceDistCoeffSqr_fine_iteration_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::cvControls::surfacePtReplaceDistCoeffSqrIteration"
            "("
                "int reconfMode"
            ") const"
        )   << "Unknown reconformationMode " << reconfMode
            << exit(FatalError);
    }

    return 0;
}


Foam::label Foam::cvControls::maxConformationIterations
(
    int reconfMode
) const
{
    if (reconfMode == conformalVoronoiMesh::rmCoarse)
    {
        return maxConformationIterations_coarse_;
    }

    else if (reconfMode == conformalVoronoiMesh::rmFine)
    {
        return maxConformationIterations_fine_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::cvControls::maxConformationIterations"
            "("
                "int reconfMode"
            ") const"
        )   << "Unknown reconformationMode " << reconfMode
            << exit(FatalError);
    }

    return 0;
}


Foam::scalar Foam::cvControls::iterationToIntialHitRatioLimit
(
    int reconfMode
) const
{
    if (reconfMode == conformalVoronoiMesh::rmCoarse)
    {
        return iterationToIntialHitRatioLimit_coarse_;
    }

    else if (reconfMode == conformalVoronoiMesh::rmFine)
    {
        return iterationToIntialHitRatioLimit_fine_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::cvControls::iterationToIntialHitRatioLimit"
            "("
                "int reconfMode"
            ") const"
        )   << "Unknown reconformationMode " << reconfMode
            << exit(FatalError);
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cvControls::~cvControls()
{}


// ************************************************************************* //
