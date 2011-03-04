/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "NonInertialFrameForce.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonInertialFrameForce<CloudType>::NonInertialFrameForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict),
    WName_
    (
        this->coeffs().template lookupOrDefault<word>
        (
            "linearAccelerationName",
            "linearAcceleration"
        )
    ),
    W_(vector::zero),
    omegaMagName_
    (
        this->coeffs().template lookupOrDefault<word>
        (
            "angularVelocityMagName",
            "angularVelocityMag"
        )
    ),
    omegaMag_(0.0),
    omegaDotName_
    (
        this->coeffs().template lookupOrDefault<word>
        (
            "angularAccelerationName",
            "angularAcceleration"
        )
    ),
    omegaDot_(vector::zero),
    axisName_
    (
        this->coeffs().template lookupOrDefault<word>
        (
            "axisName",
            "axis"
        )
    ),
    axis_(vector::zero),
    hasAxis_(false),
    axisRefPointName_
    (
        this->coeffs().template lookupOrDefault<word>
        (
            "axisRefPointName",
            "axisRefPoint"
        )
    ),
    axisRefPoint_(vector::zero)
{}


template<class CloudType>
Foam::NonInertialFrameForce<CloudType>::NonInertialFrameForce
(
    const NonInertialFrameForce& niff
)
:
    ParticleForce<CloudType>(niff),
    WName_(niff.WName_),
    W_(niff.W_),
    omegaMagName_(niff.omegaMagName_),
    omegaMag_(niff.omegaMag_),
    omegaDotName_(niff.omegaDotName_),
    omegaDot_(niff.omegaDot_),
    axisName_(niff.axisName_),
    axis_(niff.axis_),
    hasAxis_(niff.hasAxis_),
    axisRefPointName_(niff.axisRefPointName_),
    axisRefPoint_(niff.axisRefPoint_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonInertialFrameForce<CloudType>::~NonInertialFrameForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::NonInertialFrameForce<CloudType>::cacheFields(const bool store)
{
    W_ = vector::zero;
    omegaMag_ = 0.0;
    omegaDot_ = vector::zero;
    axis_ = vector::zero;
    hasAxis_ = false;
    axisRefPoint_ = vector::zero;

    if (store)
    {
        if
        (
            this->mesh().template
                foundObject<uniformDimensionedScalarField>(omegaMagName_)
        )
        {
            uniformDimensionedScalarField omegaMag = this->mesh().template
                lookupObject<uniformDimensionedScalarField>(omegaMagName_);

            omegaMag_ = omegaMag.value();

            // If omegaMag is found, require that the axis and axisRefPoint is
            // found.
            uniformDimensionedVectorField a = this->mesh().template
                lookupObject<uniformDimensionedVectorField>(axisName_);

            axis_ = a.value();

            hasAxis_ = true;

            scalar axisMag = mag(axis_);

            if (mag(axis_) < SMALL)
            {
                FatalErrorIn
                (
                    "void Foam::NonInertialFrameForce<CloudType>::"
                    "cacheFields(const bool store)"
                )   << axisName_ << " " << axis_ << " too small."
                    << abort(FatalError);
            }

            axis_ /= axisMag;

            uniformDimensionedVectorField axisRefPoint = this->mesh().template
                lookupObject<uniformDimensionedVectorField>(axisRefPointName_);

            axisRefPoint_ = axisRefPoint.value();

            // Only look for omegaDot is omegaMag is found, optional.
            if
            (
                this->mesh().template
                    foundObject<uniformDimensionedVectorField>(omegaDotName_)
            )
            {
                uniformDimensionedVectorField omegaDot = this->mesh().template
                    lookupObject<uniformDimensionedVectorField>(omegaDotName_);

                omegaDot_ = omegaDot.value();
            }
        }

        if
        (
            this->mesh().template
                foundObject<uniformDimensionedVectorField>(WName_)
        )
        {
            uniformDimensionedVectorField W = this->mesh().template
                lookupObject<uniformDimensionedVectorField>(WName_);

            W_ = W.value();
        }
        else if (!hasAxis_)
        {
            WarningIn
            (
                "void Foam::NonInertialFrameForce<CloudType>::"
                "cacheFields(const bool store)"
            )  << "No " << typeName << " variables found." << endl;
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::NonInertialFrameForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    value.Su() += -mass*W_;

    if (hasAxis_)
    {
        const vector pRel = p.position() - axisRefPoint_;

        const vector r = pRel - axis_*(axis_ & pRel);

        vector omega = axis_*omegaMag_;

        value.Su() +=
            mass
           *(
               (r ^ omegaDot_)
             + 2.0*(p.U() ^ omega)
             + (omega ^ (r ^ omega))
            );
    }

    return value;
}


// ************************************************************************* //
