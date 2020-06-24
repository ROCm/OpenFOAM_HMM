/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(partialFaceAreaWeightAMI, 0);
    addToRunTimeSelectionTable
    (
        AMIInterpolation,
        partialFaceAreaWeightAMI,
        dict
    );
    addToRunTimeSelectionTable
    (
        AMIInterpolation,
        partialFaceAreaWeightAMI,
        component
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::partialFaceAreaWeightAMI::setNextFaces
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
    return faceAreaWeightAMI::setNextFaces
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

Foam::partialFaceAreaWeightAMI::partialFaceAreaWeightAMI
(
    const dictionary& dict,
    const bool reverseTarget
)
:
    faceAreaWeightAMI(dict, reverseTarget)
{}


Foam::partialFaceAreaWeightAMI::partialFaceAreaWeightAMI
(
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection,
    const faceAreaIntersect::triangulationMode triMode,
    const bool restartUncoveredSourceFace
)
:
    faceAreaWeightAMI
    (
        requireMatch,
        reverseTarget,
        lowWeightCorrection,
        triMode,
        restartUncoveredSourceFace
    )
{}


Foam::partialFaceAreaWeightAMI::partialFaceAreaWeightAMI
(
    const partialFaceAreaWeightAMI& ami
)
:
    faceAreaWeightAMI(ami)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::partialFaceAreaWeightAMI::conformal() const
{
    return false;
}


bool Foam::partialFaceAreaWeightAMI::calculate
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr
)
{
    if (faceAreaWeightAMI::calculate(srcPatch, tgtPatch, surfPtr))
    {
        if (distributed())
        {
            scalarList newTgtMagSf(std::move(tgtMagSf_));

            // Assign default sizes. Override selected values with calculated
            // values. This is to support ACMI where some of the target faces
            // are never used (so never get sent over and hence never assigned
            // to)
            tgtMagSf_ = tgtPatch0().magFaceAreas();

            for (const labelList& smap : this->extendedTgtMapPtr_->subMap())
            {
                UIndirectList<scalar>(tgtMagSf_, smap) =
                    UIndirectList<scalar>(newTgtMagSf, smap);
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
