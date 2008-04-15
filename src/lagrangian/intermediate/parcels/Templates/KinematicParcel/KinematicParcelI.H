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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
inline Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    rho0_(dimensionedScalar(dict.lookup("rho0")).value()),
    minParticleMass_(dimensionedScalar(dict.lookup("minParticleMass")).value())
{}


template <class ParcelType>
inline Foam::KinematicParcel<ParcelType>::trackData::trackData
(
    KinematicCloud<ParcelType>& cloud,
    const constantProperties& constProps,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& muInterp,
    const vector& g
)
:
    Particle<ParcelType>::trackData(cloud),
    cloud_(cloud),
    constProps_(constProps),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    muInterp_(muInterp),
    g_(g)
{}


template <class ParcelType>
inline Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    KinematicCloud<ParcelType>& owner,
    const label typeId,
    const vector& position,
    const label celli,
    const scalar d0,
    const vector& U0,
    const scalar nParticle0,
    const constantProperties& constProps
)
:
    Particle<ParcelType>(owner, position, celli),
    typeId_(typeId),
    d_(d0),
    U_(U0),
    Ur_(vector::zero),
    nParticle_(nParticle0),
    rho_(constProps.rho0()),
    tTurb_(0.0),
    UTurb_(vector::zero)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template <class ParcelType>
inline const Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_;
}


template <class ParcelType>
inline const Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::minParticleMass() const
{
    return minParticleMass_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template <class ParcelType>
inline Foam::KinematicCloud<ParcelType>&
Foam::KinematicParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


template <class ParcelType>
inline const typename Foam::KinematicParcel<ParcelType>::constantProperties&
Foam::KinematicParcel<ParcelType>::trackData::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::trackData::rhoInterp() const
{
    return rhoInterp_;
}


template <class ParcelType>
inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::KinematicParcel<ParcelType>::trackData::UInterp() const
{
    return UInterp_;
}


template<class ParcelType>
inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::trackData::muInterp() const
{
    return muInterp_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::trackData::g() const
{
    return g_;
}


// * * * * * * * * * * KinematicParcel Member Functions  * * * * * * * * * * //

template <class ParcelType>
inline const Foam::label& Foam::KinematicParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template <class ParcelType>
inline const Foam::scalar& Foam::KinematicParcel<ParcelType>::d() const
{
    return d_;
}


template <class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::d()
{
    return d_;
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::wallImpactDistance
(
    const vector&
) const
{
    return 0.5*d_;
}


template <class ParcelType>
inline Foam::label Foam::KinematicParcel<ParcelType>::faceInterpolation() const
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
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::U() const
{
    return U_;
}


template <class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::U()
{
    return U_;
}


template <class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::Ur() const
{
    return Ur_;
}


template <class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::Ur()
{
    return Ur_;
}


template <class ParcelType>
inline const Foam::scalar& Foam::KinematicParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}


template <class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::nParticle()
{
    return nParticle_;
}


template <class ParcelType>
inline const Foam::scalar& Foam::KinematicParcel<ParcelType>::rho() const
{
    return rho_;
}


template <class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::rho()
{
    return rho_;
}


template <class ParcelType>
inline const Foam::scalar& Foam::KinematicParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template <class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::tTurb()
{
    return tTurb_;
}


template <class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}


template <class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::UTurb()
{
    return UTurb_;
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::volume() const
{
    return mathematicalConstant::pi/6*pow(d_, 3);
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaP() const
{
    return 0.25*mathematicalConstant::pi*d_*d_;
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaS() const
{
    return mathematicalConstant::pi*d_*d_;
}


// ************************************************************************* //
