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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::IOdictionary& Foam::cvControls::cvMeshDict() const
{
    return cvMeshDict_;
}


inline Foam::scalar Foam::cvControls::span() const
{
    return span_;
}


inline Foam::scalar Foam::cvControls::spanSqr() const
{
    return spanSqr_;
}


inline bool Foam::cvControls::timeChecks() const
{
    return timeChecks_;
}


inline Foam::scalar Foam::cvControls::pointPairDistanceCoeff() const
{
    return pointPairDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::mixedFeaturePointPPDistanceCoeff() const
{
    return mixedFeaturePointPPDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::featurePointExclusionDistanceCoeff() const
{
    return featurePointExclusionDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::featureEdgeExclusionDistanceCoeff() const
{
    return featureEdgeExclusionDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::surfaceSearchDistanceCoeff() const
{
    return surfaceSearchDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::maxSurfaceProtrusionCoeff() const
{
    return maxSurfaceProtrusionCoeff_;
}


inline Foam::scalar Foam::cvControls::maxQuadAngle() const
{
    return maxQuadAngle_;
}


inline Foam::label Foam::cvControls::surfaceConformationRebuildFrequency() const
{
    return surfaceConformationRebuildFrequency_;
}


inline bool Foam::cvControls::objOutput() const
{
    return objOutput_;
}


inline Foam::label Foam::cvControls::alignmentSearchSpokes() const
{
    return alignmentSearchSpokes_;
}


inline Foam::scalar Foam::cvControls::cosAlignmentAcceptanceAngle() const
{
    return cosAlignmentAcceptanceAngle_;
}


inline Foam::scalar Foam::cvControls::insertionDistCoeff() const
{
    return insertionDistCoeff_;
}


inline Foam::scalar Foam::cvControls::faceAreaRatioCoeff() const
{
    return faceAreaRatioCoeff_;
}


inline Foam::scalar Foam::cvControls::cosInsertionAcceptanceAngle() const
{
    return cosInsertionAcceptanceAngle_;
}


inline Foam::scalar Foam::cvControls::removalDistCoeff() const
{
    return removalDistCoeff_;
}


inline Foam::scalar Foam::cvControls::filterSizeCoeff() const
{
    return filterSizeCoeff_;
}


inline Foam::scalar Foam::cvControls::mergeClosenessCoeff() const
{
    return mergeClosenessCoeff_;
}


inline Foam::scalar Foam::cvControls::surfaceStepFaceAngle() const
{
    return surfaceStepFaceAngle_;
}


inline Foam::scalar Foam::cvControls::edgeCollapseGuardFraction() const
{
    return edgeCollapseGuardFraction_;
}


inline Foam::scalar
Foam::cvControls::maxCollapseFaceToPointSideLengthCoeff() const
{
    return maxCollapseFaceToPointSideLengthCoeff_;
}

// ************************************************************************* //
