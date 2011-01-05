/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011 OpenCFD Ltd.
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

#include "NonSphereDrag.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonSphereDrag<CloudType>::NonSphereDrag
(
    const dictionary& dict,
    CloudType& owner
)
:
    DragModel<CloudType>(dict, owner, typeName),
    phi_(readScalar(this->coeffDict().lookup("phi"))),
    a_(exp(2.3288 - 6.4581*phi_ + 2.4486*sqr(phi_))),
    b_(0.0964 + 0.5565*phi_),
    c_(exp(4.9050 - 13.8944*phi_ + 18.4222*sqr(phi_) - 10.2599*pow3(phi_))),
    d_(exp(1.4681 + 12.2584*phi_ - 20.7322*sqr(phi_) + 15.8855*pow3(phi_)))
{
    if (phi_ <= 0 || phi_ > 1)
    {
        FatalErrorIn
        (
            "NonSphereDrag<CloudType>::NonSphereDrag"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )   << "Ratio of surface of sphere having same volume as particle to "
            << "actual surface area of particle (phi) must be greater than 0 "
            << "and less than or equal to 1" << exit(FatalError);
    }
}


template<class CloudType>
Foam::NonSphereDrag<CloudType>::NonSphereDrag
(
    const NonSphereDrag<CloudType>& dm
)
:
    DragModel<CloudType>(dm),
    phi_(dm.phi_),
    a_(dm.a_),
    b_(dm.b_),
    c_(dm.c_),
    d_(dm.d_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonSphereDrag<CloudType>::~NonSphereDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::NonSphereDrag<CloudType>::Cd(const scalar Re) const
{
    return 24.0/(Re + ROOTVSMALL)*(1.0 + a_*pow(Re, b_)) + Re*c_/(Re + d_);
}


// ************************************************************************* //
