/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "ConeInjection.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjection<CloudType>::ConeInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    SOI_(readScalar(coeffDict_.lookup("SOI"))),
    duration_(readScalar(coeffDict_.lookup("duration"))),
    position_(coeffDict_.lookup("position")),
    direction_(coeffDict_.lookup("direction")),
    parcelsPerSecond_(readScalar(coeffDict_.lookup("parcelsPerSecond"))),
    volumeFlowRate_
    (
        DataEntry<scalar>::New
        (
            "volumeFlowRate",
            coeffDict_
        )
    ),
    Umag_
    (
        DataEntry<scalar>::New
        (
            "Umag",
            coeffDict_
        )
    ),
    thetaInner_
    (
        DataEntry<scalar>::New
        (
            "thetaInner",
            coeffDict_
        )
    ),
    thetaOuter_
    (
        DataEntry<scalar>::New
        (
            "thetaOuter",
            coeffDict_
        )
    ),
    parcelPDF_
    (
        pdf::New
        (
            coeffDict_.subDict("parcelPDF"),
            owner.rndGen()
        )
    ),
    tanVec1_(vector::zero),
    tanVec2_(vector::zero)
//    nParticlesPerParcel_(0.0)
{
    // Normalise direction vector
    direction_ /= mag(direction_);

    // Determine direction vectors tangential to direction
    vector tangent = vector::zero;
    scalar magTangent = 0.0;

    Random rnd(label(0));

    while (magTangent < SMALL)
    {
        vector v = rnd.vector01();

        tangent = v - (v & direction_)*direction_;
        magTangent = mag(tangent);
    }

    tanVec1_ = tangent/magTangent;
    tanVec2_ = direction_^tanVec1_;

    Info<< "tanVec1_ = " << tanVec1_ << endl;
    Info<< "tanVec2_ = " << tanVec2_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjection<CloudType>::~ConeInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ConeInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::timeStart() const
{
    return SOI_;
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::timeEnd() const
{
    return SOI_ + duration_;
}


template<class CloudType>
Foam::vector Foam::ConeInjection<CloudType>::position
(
    const label,
    const scalar,
    const polyMeshInfo& meshInfo,
    Random&
) const
{
    vector pos = position_;
    if (meshInfo.caseIs2d())
    {
        if (meshInfo.caseIs2dWedge())
        {
            pos.component(meshInfo.emptyComponent()) = 0.0;
        }
        else if (meshInfo.caseIs2dSlab())
        {
            pos.component(meshInfo.emptyComponent()) =
                meshInfo.centrePoint().component(meshInfo.emptyComponent());
        }
        else
        {
            FatalErrorIn
            (
                "Foam::vector Foam::ConeInjection<CloudType>::position"
            )   << "Could not determine 2-D case geometry" << nl
                << abort(FatalError);
        }
    }

    return pos;
}


template<class CloudType>
Foam::label Foam::ConeInjection<CloudType>::nParcelsToInject
(
    const label,
    const scalar time0,
    const scalar time1
) const
{
    scalar t0 = time0 - SOI_;
    scalar t1 = time1 - SOI_;

    return round((t1 - t0)*parcelsPerSecond_);
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::volume
(
    const scalar time0,
    const scalar time1,
    const polyMeshInfo&
) const
{
    scalar t0 = time0 - SOI_;
    scalar t1 = time1 - SOI_;

    return volumeFlowRate_().integrate(t0, t1);
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::volumeFraction
(
    const scalar time0,
    const scalar time1
) const
{
    scalar t0 = time0 - SOI_;
    scalar t1 = time1 - SOI_;

    return
        volumeFlowRate_().integrate(t0, t1)
       /volumeFlowRate_().integrate(SOI_, SOI_ + duration_);
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::d0
(
    const label,
    const scalar
) const
{
    return parcelPDF_().sample();
}


template<class CloudType>
Foam::vector Foam::ConeInjection<CloudType>::velocity
(
    const label,
    const scalar time,
    const polyMeshInfo& meshInfo,
    Random& rndGen
) const
{
    const scalar deg2Rad = mathematicalConstant::pi/180.0;
    scalar t = time - SOI_;
    scalar ti = thetaInner_().value(t);
    scalar to = thetaOuter_().value(t);
    scalar coneAngle = deg2Rad*(rndGen.scalar01()*(to - ti) + ti);

    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);
    scalar beta = 2.0*mathematicalConstant::pi*rndGen.scalar01();

    vector normal = alpha*(tanVec1_*cos(beta) + tanVec2_*sin(beta));
    vector dirVec = dcorr*direction_;
    dirVec += normal;
    if (meshInfo.caseIs2dSlab())
    {
        dirVec.component(meshInfo.emptyComponent()) = 0.0;
    }

    dirVec /= mag(dirVec);

    return Umag_().value(t)*dirVec;
}


// ************************************************************************* //
