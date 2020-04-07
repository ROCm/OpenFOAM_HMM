/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "partialFaceAreaWeightAMI.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::setNextFaces
(
    label& startSeedi,
    label& srcFacei,
    label& tgtFacei,
    const bitSet& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces,
    const bool errorOnNotFound
) const
{
    faceAreaWeightAMI<SourcePatch, TargetPatch>::setNextFaces
    (
        startSeedi,
        srcFacei,
        tgtFacei,
        mapFlag,
        seedFaces,
        visitedFaces,
        false // no error on not found
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::
partialFaceAreaWeightAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch
)
:
    faceAreaWeightAMI<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        triMode,
        reverseTarget,
        requireMatch,
        true // false // Not performing restart on low weights - valid for partial match
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
bool Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::conformal() const
{
    return false;
}


template<class SourcePatch, class TargetPatch>
bool Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    pointListList& srcCentroids,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    scalarList& srcMagSf,
    scalarList& tgtMagSf,
    autoPtr<mapDistribute>& srcMapPtr,
    autoPtr<mapDistribute>& tgtMapPtr,
    label srcFacei,
    label tgtFacei
)
{
     bool ok =
        faceAreaWeightAMI<SourcePatch, TargetPatch>::calculate
        (
            srcAddress,
            srcWeights,
            srcCentroids,
            tgtAddress,
            tgtWeights,
            srcMagSf,
            tgtMagSf,
            srcMapPtr,
            tgtMapPtr,
            srcFacei,
            tgtFacei
        );

    if (ok)
    {
        if (this->distributed())
        {
            scalarList newTgtMagSf(std::move(tgtMagSf));

            // Assign default sizes. Override selected values with
            // calculated values. This is to support ACMI
            // where some of the target faces are never used (so never get sent
            // over and hence never assigned to)
            tgtMagSf = this->tgtPatch0_.magFaceAreas();

            for (const labelList& smap : this->extendedTgtMapPtr_->subMap())
            {
                UIndirectList<scalar>(tgtMagSf, smap) =
                    UIndirectList<scalar>(newTgtMagSf, smap);
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
