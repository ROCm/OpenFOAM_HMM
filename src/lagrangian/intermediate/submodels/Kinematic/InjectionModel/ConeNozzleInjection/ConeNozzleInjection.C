/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
#include "TimeFunction1.H"
#include "unitConversion.H"
#include "distributionModel.H"

using namespace Foam::constant;


template<class CloudType>
const Foam::Enum
<
    typename Foam::ConeNozzleInjection<CloudType>::injectionMethod
>
Foam::ConeNozzleInjection<CloudType>::injectionMethodNames
({
    { injectionMethod::imPoint, "point" },
    { injectionMethod::imDisc, "disc" },
    { injectionMethod::imMovingPoint, "movingPoint" },
});

template<class CloudType>
const Foam::Enum
<
    typename Foam::ConeNozzleInjection<CloudType>::flowType
>
Foam::ConeNozzleInjection<CloudType>::flowTypeNames
({
    { flowType::ftConstantVelocity, "constantVelocity" },
    { flowType::ftPressureDrivenVelocity, "pressureDrivenVelocity" },
    { flowType::ftFlowRateAndDischarge, "flowRateAndDischarge" },
});


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setInjectionMethod()
{
    switch (injectionMethod_)
    {
        case injectionMethod::imPoint:
        case injectionMethod::imDisc:
        {
            this->coeffDict().readEntry("position", position_);
            break;
        }
        case injectionMethod::imMovingPoint:
        {
            positionVsTime_.reset(this->coeffDict());
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled injection method "
                << injectionMethodNames[injectionMethod_]
                << exit(FatalError);
        }
    }
}


template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setFlowType()
{
    switch (flowType_)
    {
        case flowType::ftConstantVelocity:
        {
            this->coeffDict().readEntry("UMag", UMag_);
            break;
        }
        case flowType::ftPressureDrivenVelocity:
        {
            Pinj_.reset(this->coeffDict());
            break;
        }
        case flowType::ftFlowRateAndDischarge:
        {
            Cd_.reset(this->coeffDict());
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled flow type "
                << flowTypeNames[flowType_]
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::ConeNozzleInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    injectionMethod_
    (
        injectionMethodNames.get("injectionMethod", this->coeffDict())
    ),
    flowType_(flowTypeNames.get("flowType", this->coeffDict())),
    outerDiameter_(this->coeffDict().getScalar("outerDiameter")),
    innerDiameter_(this->coeffDict().getScalar("innerDiameter")),
    duration_(this->coeffDict().getScalar("duration")),
    positionVsTime_(owner.db().time(), "position"),
    position_(Zero),
    injectorCell_(-1),
    tetFacei_(-1),
    tetPti_(-1),
    direction_(this->coeffDict().lookup("direction")),
    parcelsPerSecond_(this->coeffDict().getScalar("parcelsPerSecond")),
    flowRateProfile_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()
        )
    ),
    thetaInner_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "thetaInner",
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "thetaOuter",
            this->coeffDict()
        )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
    tanVec1_(Zero),
    tanVec2_(Zero),
    normal_(Zero),

    UMag_(0.0),
    Cd_(owner.db().time(), "Cd"),
    Pinj_(owner.db().time(), "Pinj")
{
    if (innerDiameter_ >= outerDiameter_)
    {
        FatalErrorInFunction
            << "Inner diameter must be less than the outer diameter:" << nl
            << "    innerDiameter: " << innerDiameter_ << nl
            << "    outerDiameter: " << outerDiameter_
            << exit(FatalError);
    }

    duration_ = owner.db().time().userTimeToTime(duration_);

    setInjectionMethod();

    setFlowType();

    Random& rndGen = this->owner().rndGen();

    // Normalise direction vector
    direction_.normalise();

    // Determine direction vectors tangential to direction
    vector tangent = Zero;
    scalar magTangent = 0.0;

    while(magTangent < SMALL)
    {
        vector v = rndGen.globalSample01<vector>();

        tangent = v - (v & direction_)*direction_;
        magTangent = mag(tangent);
    }

    tanVec1_ = tangent/magTangent;
    tanVec2_ = direction_^tanVec1_;

    // Set total volume to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_);

    updateMesh();
}


template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::ConeNozzleInjection
(
    const ConeNozzleInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    injectionMethod_(im.injectionMethod_),
    flowType_(im.flowType_),
    outerDiameter_(im.outerDiameter_),
    innerDiameter_(im.innerDiameter_),
    duration_(im.duration_),
    positionVsTime_(im.positionVsTime_),
    position_(im.position_),
    injectorCell_(im.injectorCell_),
    tetFacei_(im.tetFacei_),
    tetPti_(im.tetPti_),
    direction_(im.direction_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    thetaInner_(im.thetaInner_),
    thetaOuter_(im.thetaOuter_),
    sizeDistribution_(im.sizeDistribution_.clone()),
    tanVec1_(im.tanVec1_),
    tanVec2_(im.tanVec2_),
    normal_(im.normal_),
    UMag_(im.UMag_),
    Cd_(im.Cd_),
    Pinj_(im.Pinj_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeNozzleInjection<CloudType>::~ConeNozzleInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::updateMesh()
{
    // Set/cache the injector cell info for static methods

    switch (injectionMethod_)
    {
        case injectionMethod::imPoint:
        {
            this->findCellAtPosition
            (
                injectorCell_,
                tetFacei_,
                tetPti_,
                position_
            );
            break;
        }
        default:
        {
            // Do nothing for the other methods
        }
    }
}


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

    return 0;
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
        return flowRateProfile_.integrate(time0, time1);
    }

    return 0.0;
}


template<class CloudType>
void Foam::ConeNozzleInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    Random& rndGen = this->owner().rndGen();

    scalar beta = mathematical::twoPi*rndGen.globalSample01<scalar>();
    normal_ = tanVec1_*cos(beta) + tanVec2_*sin(beta);

    switch (injectionMethod_)
    {
        case injectionMethod::imPoint:
        {
            position = position_;
            cellOwner = injectorCell_;
            tetFacei = tetFacei_;
            tetPti = tetPti_;

            break;
        }
        case injectionMethod::imMovingPoint:
        {
            position = positionVsTime_.value(time - this->SOI_);

            this->findCellAtPosition
            (
                cellOwner,
                tetFacei,
                tetPti,
                position
            );

            break;
        }
        case injectionMethod::imDisc:
        {
            scalar frac = rndGen.globalSample01<scalar>();
            scalar dr = outerDiameter_ - innerDiameter_;
            scalar r = 0.5*(innerDiameter_ + frac*dr);

            position = position_ + r*normal_;

            this->findCellAtPosition
            (
                cellOwner,
                tetFacei,
                tetPti,
                position
            );
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled injection method "
                << injectionMethodNames[injectionMethod_]
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
    Random& rndGen = this->owner().rndGen();

    // Set particle velocity
    scalar t = time - this->SOI_;
    scalar ti = thetaInner_.value(t);
    scalar to = thetaOuter_.value(t);
    scalar coneAngle = degToRad(rndGen.sample01<scalar>()*(to - ti) + ti);

    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);

    vector normal = alpha*normal_;
    vector dirVec = dcorr*direction_;
    dirVec += normal;
    dirVec.normalise();

    switch (flowType_)
    {
        case flowType::ftConstantVelocity:
        {
            parcel.U() = UMag_*dirVec;
            break;
        }
        case flowType::ftPressureDrivenVelocity:
        {
            scalar pAmbient = this->owner().pAmbient();
            scalar rho = parcel.rho();
            scalar UMag = ::sqrt(2.0*(Pinj_.value(t) - pAmbient)/rho);
            parcel.U() = UMag*dirVec;
            break;
        }
        case flowType::ftFlowRateAndDischarge:
        {
            scalar Ao = 0.25*mathematical::pi*outerDiameter_*outerDiameter_;
            scalar Ai = 0.25*mathematical::pi*innerDiameter_*innerDiameter_;
            scalar massFlowRate =
                this->massTotal()
               *flowRateProfile_.value(t)
               /this->volumeTotal();

            scalar Umag = massFlowRate/(parcel.rho()*Cd_.value(t)*(Ao - Ai));
            parcel.U() = Umag*dirVec;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled injection method "
                << flowTypeNames[flowType_]
                << exit(FatalError);
        }
    }

    // Set particle diameter
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
