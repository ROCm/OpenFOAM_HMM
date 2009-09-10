/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "mathConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
inline Foam::InteractingKinematicParcel<ParcelType>::
constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    dict_(parentDict.subDict("constantProperties")),
    rhoMin_(dimensionedScalar(dict_.lookup("rhoMin")).value()),
    rho0_(dimensionedScalar(dict_.lookup("rho0")).value()),
    minParticleMass_
    (
        dimensionedScalar(dict_.lookup("minParticleMass")).value()
    )
{}


template <class ParcelType>
inline Foam::InteractingKinematicParcel<ParcelType>::trackData::trackData
(
    InteractingKinematicCloud<ParcelType>& cloud,
    const constantProperties& constProps,
    const interpolation<scalar>& rhoInterp,
    const interpolation<vector>& UInterp,
    const interpolation<scalar>& muInterp,
    const vector& g,
    trackPart part
)
:
    Particle<ParcelType>::trackData(cloud),
    cloud_(cloud),
    constProps_(constProps),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    muInterp_(muInterp),
    g_(g),
    part_(part)
{}


template <class ParcelType>
inline Foam::InteractingKinematicParcel<ParcelType>::InteractingKinematicParcel
(
    InteractingKinematicCloud<ParcelType>& owner,
    const vector& position,
    const label cellI
)
:
    Particle<ParcelType>(owner, position, cellI),
    typeId_(owner.parcelTypeId()),
    nParticle_(0),
    d_(0.0),
    U_(vector::zero),
    f_(vector::zero),
    pi_(vector::zero),
    tau_(vector::zero),
    rho_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    collisionRecords_(),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{}


template <class ParcelType>
inline Foam::InteractingKinematicParcel<ParcelType>::InteractingKinematicParcel
(
    InteractingKinematicCloud<ParcelType>& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const vector& f0,
    const vector& pi0,
    const vector& tau0,
    const constantProperties& constProps
)
:
    Particle<ParcelType>(owner, position, cellI),
    typeId_(typeId),
    nParticle_(nParticle0),
    d_(d0),
    U_(U0),
    f_(f0),
    pi_(pi0),
    tau_(tau0),
    rho_(constProps.rho0()),
    tTurb_(0.0),
    UTurb_(vector::zero),
    collisionRecords_(),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template <class ParcelType>
inline const Foam::dictionary&
Foam::InteractingKinematicParcel<ParcelType>::constantProperties::dict() const
{
    return dict_;
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::constantProperties::rhoMin() const
{
    return rhoMin_;
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_;
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::
constantProperties::minParticleMass() const
{
    return minParticleMass_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template <class ParcelType>
inline Foam::InteractingKinematicCloud<ParcelType>&
Foam::InteractingKinematicParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


template <class ParcelType>
inline
const typename Foam::InteractingKinematicParcel<ParcelType>::constantProperties&
Foam::InteractingKinematicParcel<ParcelType>::trackData::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::InteractingKinematicParcel<ParcelType>::trackData::rhoInterp() const
{
    return rhoInterp_;
}


template <class ParcelType>
inline const Foam::interpolation<Foam::vector>&
Foam::InteractingKinematicParcel<ParcelType>::trackData::UInterp() const
{
    return UInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::InteractingKinematicParcel<ParcelType>::trackData::muInterp() const
{
    return muInterp_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::InteractingKinematicParcel<ParcelType>::trackData::g() const
{
    return g_;
}


template<class ParcelType>
inline
typename Foam::InteractingKinematicParcel<ParcelType>::trackData::trackPart
Foam::InteractingKinematicParcel<ParcelType>::trackData::part() const
{
    return part_;
}


template<class ParcelType>
inline
typename Foam::InteractingKinematicParcel<ParcelType>::trackData::trackPart&
Foam::InteractingKinematicParcel<ParcelType>::trackData::part()
{
    return part_;
}


// * * * * * * * InteractingKinematicParcel Member Functions  * * * * * * * //

template <class ParcelType>
inline Foam::label Foam::InteractingKinematicParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::d() const
{
    return d_;
}


template <class ParcelType>
inline const Foam::vector&
Foam::InteractingKinematicParcel<ParcelType>::U() const
{
    return U_;
}


template <class ParcelType>
inline const Foam::vector&
Foam::InteractingKinematicParcel<ParcelType>::f() const
{
    return f_;
}


template <class ParcelType>
inline const Foam::vector&
Foam::InteractingKinematicParcel<ParcelType>::pi() const
{
    return pi_;
}


template <class ParcelType>
inline const Foam::vector&
Foam::InteractingKinematicParcel<ParcelType>::tau() const
{
    return tau_;
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::rho() const
{
    return rho_;
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template <class ParcelType>
inline const Foam::vector&
Foam::InteractingKinematicParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}


template <class ParcelType>
inline const Foam::collisionRecordList&
Foam::InteractingKinematicParcel<ParcelType>::collisionRecords() const
{
    return collisionRecords_;
}


template <class ParcelType>
inline Foam::label Foam::InteractingKinematicParcel<ParcelType>::typeId()
{
    return typeId_;
}


template <class ParcelType>
inline Foam::scalar& Foam::InteractingKinematicParcel<ParcelType>::nParticle()
{
    return nParticle_;
}


template <class ParcelType>
inline Foam::scalar& Foam::InteractingKinematicParcel<ParcelType>::d()
{
    return d_;
}


template <class ParcelType>
inline Foam::vector& Foam::InteractingKinematicParcel<ParcelType>::U()
{
    return U_;
}


template <class ParcelType>
inline Foam::vector& Foam::InteractingKinematicParcel<ParcelType>::f()
{
    return f_;
}


template <class ParcelType>
inline Foam::vector& Foam::InteractingKinematicParcel<ParcelType>::pi()
{
    return pi_;
}


template <class ParcelType>
inline Foam::vector& Foam::InteractingKinematicParcel<ParcelType>::tau()
{
    return tau_;
}


template <class ParcelType>
inline Foam::scalar& Foam::InteractingKinematicParcel<ParcelType>::rho()
{
    return rho_;
}


template <class ParcelType>
inline Foam::scalar& Foam::InteractingKinematicParcel<ParcelType>::tTurb()
{
    return tTurb_;
}


template <class ParcelType>
inline Foam::vector& Foam::InteractingKinematicParcel<ParcelType>::UTurb()
{
    return UTurb_;
}


template <class ParcelType>
inline Foam::collisionRecordList&
Foam::InteractingKinematicParcel<ParcelType>::collisionRecords()
{
    return collisionRecords_;
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::wallImpactDistance
(
    const vector&
) const
{
    // To make sure that the particles jump cyclic, wedge and
    // processor patches in the correct location, this must be zero,
    // or the calculation of lambda must be revised

    return 0.0;
    // return 0.5*d_;
}


template <class ParcelType>
inline Foam::label
Foam::InteractingKinematicParcel<ParcelType>::faceInterpolation() const
{
    // Use volume-based interpolation if dealing with external faces
    if (this->cloud().internalFace(this->face()))
    {
        return this->face();
    }
    else
    {
        return -1;
    }
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::massCell
(
    const label cellI
) const
{
    return rhoc_*this->cloud().pMesh().cellVolumes()[cellI];
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}


template <class ParcelType>
inline Foam::vector
Foam::InteractingKinematicParcel<ParcelType>::omega() const
{
    return pi_/momentOfInertia();
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::r() const
{
    return 0.5*d_;
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::volume() const
{
    return volume(d_);
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::volume(const scalar d) const
{
    return math::pi/6.0*pow3(d);
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::areaP() const
{
    return areaP(d_);
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::areaP(const scalar d) const
{
    return 0.25*areaS(d);
}


template <class ParcelType>
inline Foam::scalar Foam::InteractingKinematicParcel<ParcelType>::areaS() const
{
    return areaS(d_);
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::areaS(const scalar d) const
{
    return math::pi*d*d;
}


template <class ParcelType>
inline Foam::scalar
Foam::InteractingKinematicParcel<ParcelType>::Re
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar muc
) const
{
    return rhoc*mag(U - Uc_)*d/muc;
}


// ************************************************************************* //
