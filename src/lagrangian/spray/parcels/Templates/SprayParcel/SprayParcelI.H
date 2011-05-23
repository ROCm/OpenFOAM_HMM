/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ReactingParcel<ParcelType>::constantProperties
    (
        parentDict,
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("rhoMin")).value(),
        1000.0, // rho0 bypassed
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("minParticleMass")).value(),
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("T0")).value(),
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("TMin")).value(),
        4000.0, // cp0 bypassed
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("epsilon0")).value(),
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("f0")).value(),
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("Pr")).value(),
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("pMin")).value(),
        false, // constant volume
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("Tvap")).value(),
        dimensionedScalar((parentDict.subDict("constantProperties")).lookup("Tbp")).value()
    )
{

}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::trackData::trackData
(
    SprayCloud<ParcelType>& cloud,
    const constantProperties& constProps,
    const interpolation<scalar>& rhoInterp,
    const interpolation<vector>& UInterp,
    const interpolation<scalar>& muInterp,
    const interpolation<scalar>& TInterp,
    const interpolation<scalar>& CpInterp,
    const interpolation<scalar>& pInterp,
    const vector& g
)
:
    ReactingParcel<ParcelType>::trackData
    (
        cloud,
        constProps,
        rhoInterp,
        UInterp,
        muInterp,
        TInterp,
        CpInterp,
        pInterp,
        g
    ),
    cloud_(cloud),
    constProps_(constProps)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::SprayParcel
(
    SprayCloud<ParcelType>& owner,
    const vector& position,
    const label cellI
)
:
    ReactingParcel<ParcelType>(owner, position, cellI),
    d0_(this->d()),
    position0_(position),
    liquidCore_(owner.atomization().initLiquidCore()),
    KHindex_(0.0),
    y_(0.0),
    yDot_(0.0),
    tc_(0.0),
    ms_(0.0),
    injector_(1.0),
    tMom_(GREAT),
    user_(0.0)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::SprayParcel
(
    SprayCloud<ParcelType>& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const scalarField& Y0,
    const scalar& dInit,
    const scalar& liquidCore,
    const scalar& KHindex,
    const scalar& y,
    const scalar& yDot,
    const scalar& tc,
    const scalar& ms,
    const scalar& injector,
    const scalar& tMom,
    const scalar& user,
    const constantProperties& constProps
)
:
    ReactingParcel<ParcelType>
    (
        owner,
        position,
        cellI,
        typeId,
        nParticle0,
        d0,
        U0,
        Y0,
        constProps
    ),
    d0_(dInit),
    position0_(position),
    liquidCore_(liquidCore),
    KHindex_(KHindex),
    y_(y),
    yDot_(yDot),
    tc_(tc),
    ms_(ms),
    injector_(injector),
    tMom_(tMom),
    user_(user)
{}

// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::SprayCloud<ParcelType>&
Foam::SprayParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


template<class ParcelType>
inline const typename Foam::SprayParcel<ParcelType>::constantProperties&
Foam::SprayParcel<ParcelType>::trackData::constProps() const
{
    return constProps_;
}



// * * * * * * * * * * SprayParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::d0() const
{
    return d0_;
}

template<class ParcelType>
inline const Foam::vector& Foam::SprayParcel<ParcelType>::position0() const
{
    return position0_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::liquidCore() const
{
    return liquidCore_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::KHindex() const
{
    return KHindex_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::y() const
{
    return y_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::yDot() const
{
    return yDot_;
}


template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::tc() const
{
    return tc_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::ms() const
{
    return ms_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::injector() const
{
    return injector_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::tMom() const
{
    return tMom_;
}


template<class ParcelType>
inline const Foam::scalar& Foam::SprayParcel<ParcelType>::user() const
{
    return user_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::d0()
{
    return d0_;
}

template<class ParcelType>
inline Foam::vector& Foam::SprayParcel<ParcelType>::position0()
{
    return position0_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::liquidCore()
{
    return liquidCore_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::KHindex()
{
    return KHindex_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::y()
{
    return y_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::yDot()
{
    return yDot_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::tc()
{
    return tc_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::ms()
{
    return ms_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::injector()
{
    return injector_;
}

template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::tMom()
{
    return tMom_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::user()
{
    return user_;
}
