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
    b_(dimensionedScalar(this->coeffDict().lookup("b")).value())
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

        scalar etaN = alpha_*sqrt(M*kN)*sqrt(sqrt(normalOverlapMag));

        // Normal force
        vector fN_AB =
            rHat_AB
           *(kN*pow(normalOverlapMag, b_) - etaN*(U_AB & rHat_AB));

        pA.f() += fN_AB;
        pB.f() += -fN_AB;

        vector Uslip_AB =
            U_AB
          - (U_AB & rHat_AB)*rHat_AB
          - (pA.omega() ^ (pA.r()*rHat_AB))
          - (pB.omega() ^ (pB.r()*rHat_AB));

        const scalar deltaT = this->owner().mesh().time().deltaT().value();

        // TODO retrieve tangentialOverlap from previous collision
        vector tangentialOverlap = vector::zero;

        tangentialOverlap += Uslip_AB * deltaT;

        // const scalar& etaT = etaN;

        // Tangential force
        // fT_AB =
    }
}

// + Add this force to the sum of forces for this particle, + If
// normalOverlap < 0 then there is no collision between this pair and
// any record of collision in the previous timestep and the
// accumulated value of tangentialOverlap are removed.

// + If normalOverlap > 0 then a check is made to see if these
//   particles were colliding in the previous step, if so, retrieve
//   the previous value of tangentialOverlap, if not, create a
//   collision record with tangentialOverlap = 0.

// + Calculate Delta(tangentialOverlap):
//       Delta(tangentialOverlap) = vSlip * dt
//   where dt is the current timestep and vSlip:
//       vSlip = vRel - (vRel & n)n - omega1 ^ r1*n - omega2 ^ r2*n
//   adding Delta(tangentialOverlap) to the current value of tangentialOverlap
//   for this collision pair.

// + Using the current value of tangentialOverlap for the pair,
//   calculate the tangential component of force on this particle, Ft:
//   Ft = -min(kT*mag(tangentialOverlap), mu*mag(Fn))
//      *tangentialOverlap/mag(tangentialOverlap) - etaT*vSlip
//   Where mu is the coefficient of friction (values f in table 1?),
//   kT is a function of normalOverlap, r1, r2, E1, E2, sigma1 and
//   sigma2, and etaT = etaN.

// + Add Ft and its torque to the particle, and the corresponding
//   parts to the other particle.
//   Corresponding torque
//       ((r1*-n) ^ Fn)
//   ^ is the cross product, the point of application of the
//   force relative to the particle's position (assumed to be its centre of
//   mass) is (r1*-n).
//   The other particle receives the negative of this force value and
//   calculates its torque contribution as ((r2*n) ^ -Fn).

// ************************************************************************* //
