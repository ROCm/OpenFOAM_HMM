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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& dict
)
:
    ReactingParcel<ParcelType>::constantProperties(dict),
    Ldevol_(dimensionedScalar(dict.lookup("Ldevol")).value())
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::trackData::trackData
(
    ReactingMultiphaseCloud<ParcelType>& cloud,
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
inline Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    ReactingMultiphaseCloud<ParcelType>& owner,
    const label typeId,
    const vector& position,
    const label celli,
    const scalar d0,
    const vector& U0,
    const scalar nParticle0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const scalarField& YMixture0,
    const constantProperties& constProps
)
:
    ReactingParcel<ParcelType>
    (
        owner,
        typeId,
        position,
        celli,
        d0,
        U0,
        nParticle0,
        YMixture0,
        constProps
    ),
    YGas_(YGas0),
    YLiquid_(YLiquid0),
    YSolid_(YSolid0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::Ldevol() const
{
    return Ldevol_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ReactingMultiphaseCloud<ParcelType>&
Foam::ReactingMultiphaseParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


template<class ParcelType>
inline const typename Foam::ReactingMultiphaseParcel<ParcelType>::
constantProperties&
Foam::ReactingMultiphaseParcel<ParcelType>::trackData::constProps() const
{
    return constProps_;
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YGas() const
{
    return YGas_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YGas()
{
    return YGas_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YLiquid() const
{
    return YLiquid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YLiquid()
{
    return YLiquid_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YSolid() const
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YSolid()
{
    return YSolid_;
}


// ************************************************************************* //
