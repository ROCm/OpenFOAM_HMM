/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "WallSpringSliderDashpot.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
void Foam::WallSpringSliderDashpot<CloudType>::findMinMaxProperties
(
    scalar& rMin,
    scalar& rhoMax,
    scalar& UMagMax
) const
{
    rMin = VGREAT;
    rhoMax = -VGREAT;
    UMagMax = -VGREAT;

    forAllConstIter(typename CloudType, this->owner(), iter)
    {
        const typename CloudType::parcelType& p = iter();

        // Finding minimum diameter to avoid excessive arithmetic
        rMin = min(p.d(), rMin);

        rhoMax = max(p.rho(), rhoMax);

        UMagMax = max
        (
            mag(p.U()) + mag(p.omega())*p.d()/2,
            UMagMax
        );
    }

    // Transform the minimum diameter into minimum radius
    //     rMin = dMin/2

    rMin /= 2.0;
}

template <class CloudType>
void Foam::WallSpringSliderDashpot<CloudType>::evaluateWall
(
    typename CloudType::parcelType& p,
    const point& site,
    scalar pNu,
    scalar pE,
    scalar Estar,
    scalar kN
) const
{
    vector r_PW = p.position() - site;

    scalar normalOverlapMag = p.d()/2 - mag(r_PW);

    vector rHat_PW = r_PW/(mag(r_PW) + VSMALL);

    scalar etaN = alpha_*sqrt(p.mass()*kN)*pow025(normalOverlapMag);

    vector fN_PW =
    rHat_PW
    *(kN*pow(normalOverlapMag, b_) - etaN*(p.U() & rHat_PW));

    p.f() += fN_PW;
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
    if (!(this->owner().size()))
    {
        return 1;
    }

    scalar rMin;
    scalar rhoMax;
    scalar UMagMax;

    findMinMaxProperties(rMin, rhoMax, UMagMax);

    scalar pNu = this->owner().constProps().poissonsRatio();

    scalar pE = this->owner().constProps().youngsModulus();

    scalar Estar = 1/((1 - sqr(pNu))/pE + (1 - sqr(nu_))/E_);

    // Note:  pi^(7/5)*(5/4)^(2/5) = 5.429675
    scalar minCollisionDeltaT =
        5.429675
       *rMin
       *pow(rhoMax/(Estar*sqrt(UMagMax) + VSMALL), 0.4)
       /collisionResolutionSteps_;

    return ceil(this->owner().time().deltaTValue()/minCollisionDeltaT);
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

    scalar kN = (4.0/3.0)*sqrt(p.d()/2)*Estar;

    forAll(flatSites, siteI)
    {
        evaluateWall(p, flatSites[siteI], pNu, pE, Estar, kN);
    }

    forAll(sharpSites, siteI)
    {
        // Treating sharp sites like flat sites

        evaluateWall(p, sharpSites[siteI], pNu, pE, Estar, kN);
    }
}


// ************************************************************************* //
