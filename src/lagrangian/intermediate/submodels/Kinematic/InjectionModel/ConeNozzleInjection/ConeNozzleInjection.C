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

#include "ConeNozzleInjection.H"
#include "DataEntry.H"
#include "mathematicalConstants.H"
#include "distributionModel.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::ConeNozzleInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    injectionMethod_(imPoint),
    outerNozzleDiameter_
    (
        readScalar(this->coeffDict().lookup("outerNozzleDiameter"))
    ),
    innerNozzleDiameter_
    (
        readScalar(this->coeffDict().lookup("innerNozzleDiameter"))
    ),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    position_(this->coeffDict().lookup("position")),
    injectorCell_(-1),
    tetFaceI_(-1),
    tetPtI_(-1),
    direction_(this->coeffDict().lookup("direction")),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    volumeFlowRate_
    (
        DataEntry<scalar>::New
        (
            "volumeFlowRate",
            this->coeffDict()
        )
    ),
    Cd_
    (
        DataEntry<scalar>::New
        (
            "Cd",
            this->coeffDict()
        )
    ),
    thetaInner_
    (
        DataEntry<scalar>::New
        (
            "thetaInner",
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        DataEntry<scalar>::New
        (
            "thetaOuter",
            this->coeffDict()
        )
    ),
    sizeDistribution_
    (
        distributionModels::distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
    tanVec1_(vector::zero),
    tanVec2_(vector::zero),
    normal_(vector::zero)
{
    if (innerNozzleDiameter_ >= outerNozzleDiameter_)
    {
        FatalErrorIn
        (
            "Foam::ConeNozzleInjection<CloudType>::ConeNozzleInjection"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )<< "innerNozzleDiameter >= outerNozzleDiameter" << nl
         << exit(FatalError);
    }

    word injectionMethodType = this->coeffDict().lookup("injectionMethod");

    if (injectionMethodType == "disc")
    {
        injectionMethod_ = imDisc;
    }
    else if (injectionMethodType == "point")
    {
        injectionMethod_ = imPoint;

        // Set/cache the injector cell
        this->findCellAtPosition
        (
            injectorCell_,
            tetFaceI_,
            tetPtI_,
            position_,
            false
        );
    }
    else
    {
        FatalErrorIn
        (
            "Foam::InjectionModel<CloudType>::InjectionModel"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )<< "injectionMethod must be either 'point' or 'disc'" << nl
         << exit(FatalError);
    }

    cachedRandom& rndGen = this->owner().rndGen();

    // Normalise direction vector
    direction_ /= mag(direction_);

    // Determine direction vectors tangential to direction
    vector tangent = vector::zero;
    scalar magTangent = 0.0;

    while(magTangent < SMALL)
    {
        vector v = rndGen.sample01<vector>();

        tangent = v - (v & direction_)*direction_;
        magTangent = mag(tangent);
    }

    tanVec1_ = tangent/magTangent;
    tanVec2_ = direction_^tanVec1_;

    // Set total volume to inject
    this->volumeTotal_ = volumeFlowRate_().integrate(0.0, duration_);
}


template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::ConeNozzleInjection
(
    const ConeNozzleInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    injectionMethod_(im.injectionMethod_),
    outerNozzleDiameter_(im.outerNozzleDiameter_),
    innerNozzleDiameter_(im.innerNozzleDiameter_),
    duration_(im.duration_),
    position_(im.position_),
    injectorCell_(im.injectorCell_),
    direction_(im.direction_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    volumeFlowRate_(im.volumeFlowRate_().clone().ptr()),
    Cd_(im.Cd_().clone().ptr()),
    thetaInner_(im.thetaInner_().clone().ptr()),
    thetaOuter_(im.thetaOuter_().clone().ptr()),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    tanVec1_(im.tanVec1_),
    tanVec2_(im.tanVec1_),
    normal_(im.normal_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::~ConeNozzleInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ConeNozzleInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::ConeNozzleInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return floor((time1 - time0)*parcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ConeNozzleInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return volumeFlowRate_().integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
)
{
    cachedRandom& rndGen = this->owner().rndGen();

    scalar beta = mathematical::twoPi*rndGen.sample01<scalar>();
    normal_ = tanVec1_*cos(beta) + tanVec2_*sin(beta);

    switch (injectionMethod_)
    {
        case imPoint:
        {
            position = position_;
            cellOwner = injectorCell_;
            tetFaceI = tetFaceI_;
            tetPtI = tetPtI_;

            break;
        }
        case imDisc:
        {
            scalar frac = rndGen.sample01<scalar>();
            scalar dr = outerNozzleDiameter_ - innerNozzleDiameter_;
            scalar r = 0.5*(innerNozzleDiameter_ + frac*dr);
            position = position_ + r*normal_;

            this->findCellAtPosition
            (
                cellOwner,
                tetFaceI,
                tetPtI,
                position,
                false
            );
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void Foam::ConeNozzleInjection<CloudType>::setPositionAndCell"
                "("
                    "const label, "
                    "const label, "
                    "const scalar, "
                    "vector&, "
                    "label&"
                ")"
            )<< "Unknown injectionMethod type" << nl
             << exit(FatalError);
        }
    }
}


template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    cachedRandom& rndGen = this->owner().rndGen();

    // set particle velocity
    const scalar deg2Rad = mathematical::pi/180.0;

    scalar t = time - this->SOI_;
    scalar ti = thetaInner_().value(t);
    scalar to = thetaOuter_().value(t);
    scalar coneAngle = rndGen.sample01<scalar>()*(to - ti) + ti;

    coneAngle *= deg2Rad;
    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);

    vector normal = alpha*normal_;
    vector dirVec = dcorr*direction_;
    dirVec += normal;
    dirVec /= mag(dirVec);

    scalar Ao = 0.25*mathematical::pi*outerNozzleDiameter_*outerNozzleDiameter_;
    scalar Ai = 0.25*mathematical::pi*innerNozzleDiameter_*innerNozzleDiameter_;
    scalar massFlowRate =
        this->massTotal()*volumeFlowRate_().value(t)/this->volumeTotal();

    scalar Umag = massFlowRate/(parcel.rho()*Cd_().value(t)*(Ao - Ai));
    parcel.U() = Umag*dirVec;

    // set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::ConeNozzleInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ConeNozzleInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
