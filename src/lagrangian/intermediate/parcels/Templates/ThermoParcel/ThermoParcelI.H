/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties()
:
    KinematicParcel<ParcelType>::constantProperties(),
    T0_(0.0),
    TMin_(0.0),
    Cp0_(0.0),
    epsilon0_(0.0),
    f0_(0.0),
    Pr_(0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    KinematicParcel<ParcelType>::constantProperties(cp),
    T0_(cp.T0_),
    TMin_(cp.TMin_),
    Cp0_(cp.Cp0_),
    epsilon0_(cp.epsilon0_),
    f0_(cp.f0_),
    Pr_(cp.Pr_)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    KinematicParcel<ParcelType>::constantProperties(parentDict),
    T0_(readScalar(this->dict().lookup("T0"))),
    TMin_(readScalar(this->dict().lookup("TMin"))),
    Cp0_(readScalar(this->dict().lookup("Cp0"))),
    epsilon0_(readScalar(this->dict().lookup("epsilon0"))),
    f0_(readScalar(this->dict().lookup("f0"))),
    Pr_(readScalar(this->dict().lookup("Pr")))
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::trackData::trackData
(
    ThermoCloud<ParcelType>& cloud,
    typename KinematicParcel<ParcelType>::trackData::trackPart part
)
:
    KinematicParcel<ParcelType>::trackData(cloud, part),
    cloud_(cloud),
    Cp_
    (
        IOobject
        (
            "Cp",
            cloud.db().time().timeName(),
            cloud.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.thermo().thermo().Cp()
    ),
    TInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.T()
        )
    ),
    CpInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            Cp_
        )
    ),
    GInterp_(NULL)
{
    if (cloud.radiation())
    {
        GInterp_.reset
        (
            interpolation<scalar>::New
            (
                cloud.solution().interpolationSchemes(),
                cloud.mesh().objectRegistry::template
                    lookupObject<volScalarField>("G")
            ).ptr()
        );
    }
}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    ThermoCloud<ParcelType>& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    KinematicParcel<ParcelType>(owner, position, cellI, tetFaceI, tetPtI),
    T_(0.0),
    Cp_(0.0),
    Tc_(0.0),
    Cpc_(0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    ThermoCloud<ParcelType>& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const constantProperties& constProps
)
:
    KinematicParcel<ParcelType>
    (
        owner,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    T_(constProps.T0()),
    Cp_(constProps.Cp0()),
    Tc_(0.0),
    Cpc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::T0() const
{
    return T0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::TMin() const
{
    return TMin_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::Cp0() const
{
    return Cp0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::epsilon0() const
{
    return epsilon0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::f0() const
{
    return f0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::Pr() const
{
    return Pr_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ThermoCloud<ParcelType>&
Foam::ThermoParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::ThermoParcel<ParcelType>::trackData::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackData::TInterp() const
{
    return TInterp_();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackData::CpInterp() const
{
    return CpInterp_();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackData::GInterp() const
{
    if (!GInterp_.valid())
    {
        FatalErrorIn
        (
            "inline const Foam::interpolation<Foam::scalar>&"
            "Foam::ThermoParcel<ParcelType>::trackData::GInterp() const"
        )   << "Radiation G interpolation object not set"
            << abort(FatalError);
    }

    return GInterp_();
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::hs() const
{
    return Cp_*(T_ - 298.15);
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::T()
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::Cp()
{
    return Cp_;
}


// ************************************************************************* //
