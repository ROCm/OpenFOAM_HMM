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

#include "SpringSliderDashpot.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SpringSliderDashpot<CloudType>::SpringSliderDashpot
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PairFunction<CloudType>(dict, cloud, typeName),
    Estar_(),
    Gstar_(),
    sigma_(dimensionedScalar(this->coeffDict().lookup("sigma")).value()),
    alpha_(dimensionedScalar(this->coeffDict().lookup("alpha")).value()),
    b_(dimensionedScalar(this->coeffDict().lookup("b")).value()),
    mu_(dimensionedScalar(this->coeffDict().lookup("mu")).value())
{
    scalar E = dimensionedScalar(this->coeffDict().lookup("E")).value();

    Estar_ = E/(2.0*(1.0 - sqr(sigma_)));

    scalar G = E/(2.0*(1.0 + sigma_));

    Gstar_ = G/(2.0*(2.0 - sigma_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SpringSliderDashpot<CloudType>::~SpringSliderDashpot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SpringSliderDashpot<CloudType>::evaluatePair
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    vector r_AB = (pA.position() - pB.position());

    scalar normalOverlapMag = 0.5*(pA.d() + pB.d()) - mag(r_AB);

    if (normalOverlapMag > 0)
    {
        //Particles in collision

        vector rHat_AB = r_AB/(mag(r_AB) + VSMALL);

        vector U_AB = pA.U() - pB.U();

        // Effective radius
        scalar R = 0.5*pA.d()*pB.d()/(pA.d() + pB.d());

        // Effective mass
        scalar M = pA.mass()*pB.mass()/(pA.mass() + pB.mass());

        scalar kN = (4.0/3.0)*sqrt(R)*Estar_;

        scalar etaN = alpha_*sqrt(M*kN)*pow025(normalOverlapMag);

        // Normal force
        vector fN_AB =
            rHat_AB
           *(kN*pow(normalOverlapMag, b_) - etaN*(U_AB & rHat_AB));

        pA.f() += fN_AB;
        pB.f() += -fN_AB;

        vector USlip_AB =
            U_AB - (U_AB & rHat_AB)*rHat_AB
          + (pA.omega() ^ (pA.r()*-rHat_AB))
          - (pB.omega() ^ (pB.r()*rHat_AB));

        scalar deltaT = this->owner().mesh().time().deltaT().value();

        vector& tangentialOverlap_AB =
            pA.collisionRecords().matchRecord
            (
                pB.origProc(),
                pB.origId()
            ).collisionData();

        vector& tangentialOverlap_BA =
            pB.collisionRecords().matchRecord
            (
                pA.origProc(),
                pA.origId()
            ).collisionData();

        vector deltaTangentialOverlap_AB = USlip_AB * deltaT;

        tangentialOverlap_AB += deltaTangentialOverlap_AB;
        tangentialOverlap_BA += -deltaTangentialOverlap_AB;

        scalar kT = 8.0*sqrt(R*normalOverlapMag)*Gstar_;

        scalar& etaT = etaN;

        // Tangential force
        vector fT_AB =
            -min(kT*mag(tangentialOverlap_AB), mu_*mag(fN_AB))
           *tangentialOverlap_AB/mag(tangentialOverlap_AB)
          - etaT*USlip_AB;

        pA.f() += fT_AB;
        pB.f() += -fT_AB;

        pA.tau() += (pA.r()*-rHat_AB) ^ fT_AB;
        pB.tau() += (pB.r()*rHat_AB) ^ -fT_AB;
    }
}


// ************************************************************************* //
