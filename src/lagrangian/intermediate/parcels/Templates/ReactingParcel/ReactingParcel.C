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

#include "ReactingParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::updateCellQuantities
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    ThermoParcel<ParcelType>::updateCellQuantities(td, dt, celli);

    pc_ = td.pInterp().interpolate(this->position(), celli);
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcCoupled
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const vector U0 = this->U_;
    const scalar mass0 = this->mass();
    const scalar cp0 = this->cp_;
    const scalar np0 = this->nParticle_;
    const scalar T0 = this->T_;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Enthalpy transfer from the particle to the carrier phase
    scalar dhTrans = 0.0;

    // Mass transfer from particle to carrier phase
    // - components exist in particle already
    scalarList dMassMT(td.cloud().gases().size(), 0.0);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, Cud, dUTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    scalar T1 = calcHeatTransfer(td, dt, celli, htc, dhTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~
    // Calculate phase change
    // ~~~~~~~~~~~~~~~~~~~~~~
    scalarField X = td.cliud().composition().X(0, YMixture_);
    calcPhaseChange(td, dt, T, X, dMassMT);

    // New total mass
    const scalar mass1 = mass0 - sum(dMassMT);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Transfer mass lost from particle to carrier mass source
    forAll(dMassMT, i)
    {
        td.cloud().rhoTrans(i)[celli] += np0*dMassMT[i];
    }

    // Update momentum transfer
    td.cloud().UTrans()[celli] += np0*dUTrans;

    // Accumulate coefficient to be applied in carrier phase momentum coupling
    td.cloud().UCoeff()[celli] += np0*mass0*Cud;

    // Update enthalpy transfer
    td.cloud().hTrans()[celli] += np0*dhTrans;

    // Accumulate coefficient to be applied in carrier phase enthalpy coupling
    td.cloud().hCoeff()[celli] += np0*htc*this->areaS();


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;

        // Absorb particle(s) into carrier phase
        forAll(dMassMT, i)
        {
            td.cloud().rhoTrans(i)[celli] += np0*dMassMT[i];
        }
        td.cloud().UTrans()[celli] += np0*mass1*U1;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
    {
        this->U_ = U1;
        this->T_ = T1;
        //        this->cp_ = ??? // TODO:

        // Update particle density or diameter
        if (td.constProps().constantVolume())
        {
            this->rho_ = mass1/this->volume();
        }
        else
        {
            this->d_ = cbrt(mass1/this->rho_*6.0/mathematicalConstant::pi);
        }
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcUncoupled
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();
    const scalar cp0 = this->cp_;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Enthalpy transfer from the particle to the carrier phase
    scalar dhTrans = 0.0;

    // Mass transfer from particle to carrier phase
    // - components exist in particle already
    scalarList dMassMT(td.cloud().gases().size(), 0.0);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, Cud, dUTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    scalar T1 = calcHeatTransfer(td, dt, celli, htc, dhTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~
    // Calculate phase change
    // ~~~~~~~~~~~~~~~~~~~~~~
    scalarField X = td.cloud().composition().X(0, YMixture_);
    scalar dMassPC = calcPhaseChange(td, dt, T, X, dMassMT);
    T1 -= td.constProps().Lvap()*dMassPC/(0.5*mass0*cp0);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialise enthalpy retention to zero
    scalar dhRet = 0.0;

    // New total mass
    const scalar mass1 = mass0 - sum(dMassMT);

    // New specific heat capacity
    const scalar cp1 = cp0; // TODO: new cp1

    // Add retained enthalpy to particle
    T1 += dhRet/(mass0*0.5*(cp0 + cp1));

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
    {
        this->U_ = U1;
        this->T_ = T1;
        this->cp_ = cp1;

        // Update particle density or diameter
        if (td.constProps().constantVolume())
        {
            this->rho_ = mass1/this->volume();
        }
        else
        {
            this->d_ = cbrt(mass1/this->rho_*6.0/mathematicalConstant::pi);
        }
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcPhaseChange
(
    TrackData& td,
    const scalar dt,
    const scalar T,
    scalarField& X,
    scalarList& dMassMT
)
{
    if (!td.cloud().phaseChange().active())
    {
        return;
    }

    // TODO: separate treatment for boiling

    scalar dMassTot = td.cloud().phaseChange().calculate
    (
        T,
        this->d_,
        X,
        dMassMT,
        this->U_ - this->Uc_,
        this->Tc_,
        pc_,
        this->muc_/this->rhoc_,
        dt
    );

    td.cloud().addToMassPhaseChange(dMassTot);
    // TODO: Re-calculate mass fractions


}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //

