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
void Foam::ThermoParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    KinematicParcel<ParcelType>::setCellValues(td, dt, cellI);

    cpc_ = td.cpInterp().interpolate(this->position(), cellI);

    Tc_ = td.TInterp().interpolate(this->position(), cellI);

    if (Tc_ < td.constProps().TMin())
    {
        WarningIn
        (
            "void Foam::ThermoParcel<ParcelType>::setCellValues"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Limiting temperature in cell " << cellI << " to "
            << td.constProps().TMin() <<  nl << endl;

        Tc_ = td.constProps().TMin();
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    this->Uc_ += td.cloud().UTrans()[cellI]/this->massCell(cellI);

    scalar cpMean = td.cpInterp().psi()[cellI];
    Tc_ += td.cloud().hsTrans()[cellI]/(cpMean*this->massCell(cellI));
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

    // Initial enthalpy state
    scalar H0 = cp0*T0;


    // Heat transfer
    // ~~~~~~~~~~~~~

    // No additional enthalpy sources
    scalar Sh = 0.0;

    // Calculate new particle velocity
    scalar T1 = calcHeatTransfer(td, dt, cellI, d0, U0, rho0, T0, cp0, Sh);

    // Calculate new enthalpy state
    scalar H1 = cp0*T1;


    // Motion
    // ~~~~~~

    // No additional forces
    vector Fx = vector::zero;

    // Calculate new particle velocity
    vector U1 = calcVelocity(td, dt, cellI, d0, U0, rho0, mass0, Fx);


    //  Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*mass0*(U0 - U1);

        // Update sensible enthalpy transfer
        td.cloud().hsTrans()[cellI] += np0*mass0*(H0 - H1);
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
    const scalar Sh
)
{
    if (!td.cloud().heatTransfer().active())
    {
        return T;
    }

    // Calc heat transfer coefficient
    scalar htc = td.cloud().heatTransfer().h
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

    // Determine ap and bp coefficients
    scalar ap = Tc_ + Sh/(htc*this->areaS(d) + ROOTVSMALL);
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

    return Tres.value();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ThermoParcelIO.C"

// ************************************************************************* //
