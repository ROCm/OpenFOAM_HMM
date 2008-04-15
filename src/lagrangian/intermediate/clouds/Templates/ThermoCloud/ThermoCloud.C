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

#include "ThermoCloud.H"
#include "HeatTransferModel.H"

#include "interpolationCellPoint.H"
#include "ThermoParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoCloud<ParcelType>::ThermoCloud
(
    const word& cloudType,
    const volPointInterpolation& vpi,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    basicThermo& thermo
)
:
    KinematicCloud<ParcelType>
    (
        cloudType,
        vpi,
        rho,
        U,
        thermo.mu(),
        g
    ),
    thermoCloud(),
    constProps_(this->particleProperties()),
    carrierThermo_(thermo),
    heatTransferModel_
    (
        HeatTransferModel<ThermoCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    radiation_(this->particleProperties().lookup("radiation")),
    hTrans_
    (
        IOobject
        (
            this->cloudName() + "hTrans",
            this->runTime().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar("zero", dimensionSet(1, 2, -2, 0, 0), 0.0)
    ),
    hCoeff_
    (
        IOobject
        (
            this->cloudName() + "hCoeff",
            this->runTime().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar("zero", dimensionSet(1, 2, -3, -1, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoCloud<ParcelType>::~ThermoCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::resetSourceTerms()
{
    KinematicCloud<ParcelType>::resetSourceTerms();
    hTrans_.field() = 0.0;
    hCoeff_.field() = 0.0;
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::evolve()
{
    const volScalarField& T = carrierThermo_.T();
    const volScalarField cp = carrierThermo_.Cp();

    interpolationCellPoint<scalar> rhoInterp(this->vpi(), this->rho());
    interpolationCellPoint<vector> UInterp(this->vpi(), this->U());
    interpolationCellPoint<scalar> muInterp(this->vpi(), this->mu());
    interpolationCellPoint<scalar> TInterp(this->vpi(), T);
    interpolationCellPoint<scalar> cpInterp(this->vpi(), cp);

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterp,
        UInterp,
        muInterp,
        TInterp,
        cpInterp,
        this->g().value()
    );

    inject(td);

    move(td);
}


template<class ParcelType>
template<class TrackingData>
void Foam::ThermoCloud<ParcelType>::move
(
    TrackingData& td
)
{
    if (this->coupled())
    {
        resetSourceTerms();
    }
    Cloud<ParcelType>::move(td);
}


template<class ParcelType>
template<class TrackingData>
void Foam::ThermoCloud<ParcelType>::inject
(
    TrackingData& td
)
{
    KinematicCloud<ParcelType>::inject(td);
/*
    scalar time = this->runTime().value();

    scalar pRho = td.constProps().rho0();

    // Number of parcels to introduce during this timestep
    const label nParcels = this->injection().nParcelsToInject
    (
        this->nInjections(),
        this->time0(),
        time
    );

    // Return if no parcels are required
    if (!nParcels)
    {
        this->postInjectCheck();
        return;
    }

    // Volume of particles to introduce during this timestep
    scalar pVolume = this->injection().volume
    (
         this->time0(),
         time,
         this->meshInfo()
    );

    // Volume fraction to introduce during this timestep
    scalar pVolumeFraction =
        this->injection().volumeFraction(this->time0(), time);

    // Duration of injection period during this timestep
    scalar deltaT = min
    (
        this->runTime().deltaT().value(),
        min
        (
            time - this->injection().timeStart(),
            this->injection().timeEnd() - this->time0()
        )
    );

    // Pad injection time if injection starts during this timestep
    scalar padTime = max
    (
        0.0,
        this->injection().timeStart() - this->time0()
    );

    // Introduce new parcels linearly with time
    for (label iParcel=0; iParcel<nParcels; iParcel++)
    {
        // Calculate the pseudo time of injection for parcel 'iParcel'
        scalar timeInj = this->time0() + padTime + deltaT*iParcel/nParcels;

        // Determine injected parcel properties
        vector pPosition = this->injection().position
        (
            iParcel,
            timeInj,
            this->meshInfo(),
            this->rndGen()
        );

        // Diameter of parcels
        scalar pDiameter = this->injection().d0(iParcel, timeInj);

        // Number of particles per parcel
        scalar pNumberOfParticles = this->setNumberOfParticles
        (
            nParcels,
            pDiameter,
            pVolumeFraction,
            pRho,
            pVolume
        );

        // Velocity of parcels
        vector pU = this->injection().velocity(iParcel, timeInj);

        // Determine the injection cell
        label pCell = -1;
        this->setInjectorCellAndPosition(pCell, pPosition);

        if (pCell >= 0)
        {
            // construct the parcel that is to be injected
            ParcelType* pPtr = new ParcelType
            (
                td.cloud(),
                this->parcelTypeId(),
                pPosition,
                pCell,
                pDiameter,
                pU,
                pNumberOfParticles,
                td.constProps()
            );

            scalar dt = time - timeInj;

            pPtr->stepFraction() = (this->runTime().deltaT().value() - dt)
                /this->runTime().deltaT().value();

            this->injectParcel(td, pPtr);
         }
    }

    this->postInjectCheck();

    if (debug)
    {
        this->dumpParticlePositions();
    }
*/
}


template<class ParcelType>
template<class TrackingData>
void Foam::ThermoCloud<ParcelType>::injectParcel
(
    TrackingData& td,
    ParcelType* p
)
{
    KinematicCloud<ParcelType>::injectParcel(td, p);
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::postInjectCheck()
{
    KinematicCloud<ParcelType>::postInjectCheck();
}


// ************************************************************************* //
