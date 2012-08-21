/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "VirtualMassForce.H"
#include "fvcDdt.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VirtualMassForce<CloudType>::VirtualMassForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    UName_(this->coeffs().template lookupOrDefault<word>("U", "U")),
    Cvm_(readScalar(this->coeffs().lookup("Cvm"))),
    DUcDtPtr_(NULL)
{}


template<class CloudType>
Foam::VirtualMassForce<CloudType>::VirtualMassForce
(
    const VirtualMassForce& vmf
)
:
    ParticleForce<CloudType>(vmf),
    UName_(vmf.UName_),
    Cvm_(vmf.Cvm_),
    DUcDtPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VirtualMassForce<CloudType>::~VirtualMassForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VirtualMassForce<CloudType>::cacheFields(const bool store)
{
    if (store && !DUcDtPtr_)
    {
        const volVectorField& Uc = this->mesh().template
            lookupObject<volVectorField>(UName_);

        DUcDtPtr_ = new volVectorField(fvc::ddt(Uc) + (Uc & fvc::grad(Uc)));
    }
    else
    {
        if (DUcDtPtr_)
        {
            delete DUcDtPtr_;
            DUcDtPtr_ = NULL;
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::VirtualMassForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    value.Su() = mass*p.rhoc()/p.rho()*Cvm_*DUcDt()[p.cell()];

    return value;
}


template<class CloudType>
Foam::scalar Foam::VirtualMassForce<CloudType>::massAdd
(
    const typename CloudType::parcelType& p,
    const scalar mass
) const
{
    return mass*p.rhoc()/p.rho()*Cvm_;
}


// ************************************************************************* //
