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

template <class ParcelType>
inline Foam::DsmcParcel<ParcelType>::constantProperties::constantProperties()
:
    mass_(0),
    d_(0)
{}


template <class ParcelType>
inline Foam::DsmcParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    mass_(readScalar(dict.lookup("mass"))),
    d_(readScalar(dict.lookup("diameter")))
{}


template <class ParcelType>
inline Foam::DsmcParcel<ParcelType>::trackData::trackData
(
    DsmcCloud<ParcelType>& cloud
)
:
    Particle<ParcelType>::trackData(cloud),
    cloud_(cloud)
{}


template <class ParcelType>
inline Foam::DsmcParcel<ParcelType>::DsmcParcel
(
    DsmcCloud<ParcelType>& owner,
    const vector& position,
    const vector& U,
    const label celli,
    const label typeId
)
:
    Particle<ParcelType>(owner, position, celli),
    U_(U),
    typeId_(typeId)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template <class ParcelType>
inline Foam::scalar
Foam::DsmcParcel<ParcelType>::constantProperties::mass() const
{
    return mass_;
}


template <class ParcelType>
inline Foam::scalar
Foam::DsmcParcel<ParcelType>::constantProperties::d() const
{
    return d_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template <class ParcelType>
inline Foam::DsmcCloud<ParcelType>&
Foam::DsmcParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


// * * * * * * * * * * DsmcParcel Member Functions  * * * * * * * * * * //

template <class ParcelType>
inline Foam::label Foam::DsmcParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template <class ParcelType>
inline Foam::scalar Foam::DsmcParcel<ParcelType>::wallImpactDistance
(
    const vector&
) const
{
    Info<< "PLACEHOLDER FOR Foam::DsmcParcel<ParcelType>::wallImpactDistance"
        << endl;

    return 0;
    // return 0.5*d_;
}


template <class ParcelType>
inline const Foam::vector& Foam::DsmcParcel<ParcelType>::U() const
{
    return U_;
}


template <class ParcelType>
inline Foam::vector& Foam::DsmcParcel<ParcelType>::U()
{
    return U_;
}


// ************************************************************************* //
