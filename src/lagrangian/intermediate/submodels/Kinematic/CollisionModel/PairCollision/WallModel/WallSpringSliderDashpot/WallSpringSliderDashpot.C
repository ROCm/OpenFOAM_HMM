/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "WallSpringSliderDashpot.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
void Foam::WallSpringSliderDashpot<CloudType>::findMinMaxProperties
(
    scalar& RMin,
    scalar& rhoMax,
    scalar& UMagMax
) const
{
    RMin = VGREAT;
    rhoMax = -VGREAT;
    UMagMax = -VGREAT;

    forAllConstIter(typename CloudType, this->owner(), iter)
    {
        const typename CloudType::parcelType& p = iter();

        // Finding minimum diameter to avoid excessive arithmetic
        RMin = min(p.d(), RMin);

        rhoMax = max(p.rho(), rhoMax);

        UMagMax = max
        (
            mag(p.U()) + mag(p.omega())*p.r(),
            UMagMax
        );
    }

    // Transform the minimum diameter into minimum radius
    //     rMin = dMin/2
    // then rMin into minimum R,
    //     1/RMin = 1/rMin + 1/rMin,
    //     RMin = rMin/2 = dMin/4

    RMin /= 4.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::WallSpringSliderDashpot<CloudType>::WallSpringSliderDashpot
(
    const dictionary& dict,
    CloudType& cloud
)
:
    WallModel<CloudType>(dict, cloud, typeName),
    E_(dimensionedScalar(this->coeffDict().lookup("youngsModulus")).value()),
    nu_(dimensionedScalar(this->coeffDict().lookup("poissonsRatio")).value()),
    alpha_(dimensionedScalar(this->coeffDict().lookup("alpha")).value()),
    b_(dimensionedScalar(this->coeffDict().lookup("b")).value()),
    mu_(dimensionedScalar(this->coeffDict().lookup("mu")).value()),
    collisionResolutionSteps_
    (
        readScalar
        (
            this->coeffDict().lookup("collisionResolutionSteps")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::WallSpringSliderDashpot<CloudType>::~WallSpringSliderDashpot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::WallSpringSliderDashpot<CloudType>::controlsTimestep() const
{
    return true;
}


template<class CloudType>
Foam::label Foam::WallSpringSliderDashpot<CloudType>::nSubCycles() const
{
    // if (!(this->owner().size()))
    // {
    //     return 1;
    // }

    // scalar RMin;
    // scalar rhoMax;
    // scalar UMagMax;

    // findMinMaxProperties(RMin, rhoMax, UMagMax);

    // // Note:  pi^(7/5)*(5/4)^(2/5) = 5.429675
    // scalar minCollisionDeltaT =
    //     5.429675
    //    *RMin
    //    *pow(rhoMax/(Estar_*sqrt(UMagMax) + VSMALL), 0.4)
    //    /collisionResolutionSteps_;

    // return ceil(this->owner().time().deltaTValue()/minCollisionDeltaT);

    return 1;
}


template<class CloudType>
void Foam::WallSpringSliderDashpot<CloudType>::evaluateWall
(
    typename CloudType::parcelType& p,
    const List<point>& flatSites,
    const List<point>& sharpSites
) const
{
    scalar pNu = this->owner().constProps().poissonsRatio();

    scalar pE = this->owner().constProps().youngsModulus();

    scalar Estar = 1/((1 - sqr(pNu))/pE + (1 - sqr(nu_))/E_);

    scalar kN = (4.0/3.0)*sqrt(p.r())*Estar;

    forAll(flatSites, siteI)
    {
        vector r_PW = p.position() - flatSites[siteI];

        scalar normalOverlapMag = p.r() - mag(r_PW);

        vector rHat_PW = r_PW/(mag(r_PW) + VSMALL);

        scalar etaN = alpha_*sqrt(p.mass()*kN)*pow025(normalOverlapMag);

        vector fN_PW =
            rHat_PW
           *(kN*pow(normalOverlapMag, b_) - etaN*(p.U() & rHat_PW));

        p.f() += fN_PW;
    }

    // Treating sharp sites like flat sites

    forAll(sharpSites, siteI)
    {
        vector r_PW = p.position() - sharpSites[siteI];

        scalar normalOverlapMag = p.r() - mag(r_PW);

        vector rHat_PW = r_PW/(mag(r_PW) + VSMALL);

        scalar etaN = alpha_*sqrt(p.mass()*kN)*pow025(normalOverlapMag);

        vector fN_PW =
            rHat_PW
           *(kN*pow(normalOverlapMag, b_) - etaN*(p.U() & rHat_PW));

        p.f() += fN_PW;
    }
}


// ************************************************************************* //
