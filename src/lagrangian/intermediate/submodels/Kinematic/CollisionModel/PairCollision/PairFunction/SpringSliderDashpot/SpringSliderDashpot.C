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

    scalar normalOverlap = 0.5*(pA.d() + pB.d()) - mag(r_AB);

    if (normalOverlap > 0)
    {
        //Particles in collision

        vector rHat_AB = r_AB/(mag(r_AB) + VSMALL);

        vector U_AB = pA.U() - pB.U();

        // Effective radius
        scalar R = 0.5*pA.d()*pB.d()/(pA.d() + pB.d());

        // Effective mass
        scalar M = pA.mass()*pB.mass()/(pA.mass() + pB.mass());

        scalar kN = (4.0/3.0)*sqrt(R)*Estar_;

        scalar etaN = alpha_*sqrt(M*kN)*sqrt(sqrt(normalOverlap));

        vector f_AB =
            rHat_AB
           *(kN*pow(normalOverlap, b_) - etaN*(U_AB & rHat_AB));

        pA.f() += f_AB;
        pB.f() += -f_AB;
    }
}


// ************************************************************************* //
