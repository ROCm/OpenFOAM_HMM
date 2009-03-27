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

#include "ThermoParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::updateCellQuantities
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    KinematicParcel<ParcelType>::updateCellQuantities(td, dt, cellI);

    Tc_ = td.TInterp().interpolate(this->position(), cellI);
    if (Tc_ < SMALL)
    {
        WarningIn
        (
            "void Foam::ThermoParcel<ParcelType>::updateCellQuantities"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Temperature < " << SMALL << " in cell " << cellI << nl << endl;
    }

    cpc_ = td.cpInterp().interpolate(this->position(), cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector U0 = this->U_;
    const scalar rho0 = this->rho_;
    const scalar T0 = this->T_;
    const scalar cp0 = this->cp_;
    const scalar mass0 = this->mass();


    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    vector dUTrans = vector::zero;
    scalar dhTrans = 0.0;


    // Heat transfer
    // ~~~~~~~~~~~~~

    // No additional enthalpy sources
    scalar Sh = 0.0;

    // Calculate new particle velocity
    scalar htc = 0.0;
    scalar T1 =
        calcHeatTransfer
        (
            td,
            dt,
            cellI,
            d0,
            U0,
            rho0,
            T0,
            cp0,
            Sh,
            htc,
            dhTrans
        );


    // Motion
    // ~~~~~~

    // No additional forces
    vector Fx = vector::zero;

    // Calculate new particle velocity
    scalar Cud = 0.0;
    vector U1 =
        calcVelocity(td, dt, cellI, d0, U0, rho0, mass0, Fx, Cud, dUTrans);


    //  Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Coefficient to be applied in carrier phase momentum coupling
        td.cloud().UCoeff()[cellI] += np0*mass0*Cud;

        // Update enthalpy transfer
        td.cloud().hTrans()[cellI] += np0*dhTrans;

        // Coefficient to be applied in carrier phase enthalpy coupling
        td.cloud().hCoeff()[cellI] += np0*htc*this->areaS();
    }

    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this->U_ = U1;
    T_ = T1;
}


template<class ParcelType>
template <class TrackData>
Foam::scalar Foam::ThermoParcel<ParcelType>::calcHeatTransfer
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar d,
    const vector& U,
    const scalar rho,
    const scalar T,
    const scalar cp,
    const scalar Sh,
    scalar& htc,
    scalar& dhTrans
)
{
    if (!td.cloud().heatTransfer().active())
    {
        htc = 0.0;
        return T;
    }

    // Calc heat transfer coefficient
    htc = td.cloud().heatTransfer().h
    (
        d,
        U - this->Uc_,
        this->rhoc_,
        rho,
        cpc_,
        cp,
        this->muc_
    );


    // Determine new particle temperature
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Particle area
    scalar Ap = this->areaS(d);

    // Determine ap and bp coefficients
    scalar ap = Tc_ + Sh/(htc*Ap + ROOTVSMALL);
    scalar bp = 6.0*htc/(rho*d*cp);
    if (td.cloud().radiation())
    {
        // Carrier phase incident radiation field
        // - The G field is not interpolated to the parcel position
        //   Instead, the cell centre value is applied directly
        const scalarField& G = td.cloud().mesh().objectRegistry
            ::lookupObject<volScalarField>("G");

        // Helper variables
        const scalar sigma = radiation::sigmaSB.value();
        const scalar epsilon = td.constProps().epsilon0();
        const scalar D = epsilon*sigma*pow3(T)/(htc + ROOTVSMALL) + 1.0;
        ap += 0.25*epsilon*G[cellI]/(htc + ROOTVSMALL);
        ap /= D;
        bp *= D;
    }

    // Integrate to find the new parcel temperature
    IntegrationScheme<scalar>::integrationResult Tres =
        td.cloud().TIntegrator().integrate(T, dt, ap, bp);

    // Enthalpy transfer
    // - Using average particle temperature
    dhTrans = dt*Ap*htc*(Tres.average() - Tc_);

    return Tres.value();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ThermoParcelIO.C"

// ************************************************************************* //

