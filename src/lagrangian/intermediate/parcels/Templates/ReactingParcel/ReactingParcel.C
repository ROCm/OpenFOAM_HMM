/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

    // Mass transfer due to surface reactions from particle to carrier phase
    // - components do not necessarily exist in particle already
    scalarList dMassSR(td.cloud().gases().size(), 0.0);


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


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate mass transfer
    // ~~~~~~~~~~~~~~~~~~~~~~~
    calcMassTransfer(td, dt, T0, T1, dMassMT);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialise enthalpy retention to zero
    scalar dhRet = 0.0;

    calcSurfaceReactions(td, dt, celli, T0, T1, dMassMT, dMassSR, dhRet);

    // New total mass
    const scalar mass1 = mass0 - sum(dMassMT) - sum(dMassSR);

    // Correct particle temperature to account for latent heat of
    // devolatilisation
    T1 -=
        td.constProps().Ldevol()
       *sum(dMassMT)
       /(0.5*(mass0 + mass1)*this->cp_);

    // Add retained enthalpy from surface reaction to particle and remove
    // from gas
    T1 += dhRet/(0.5*(mass0 + mass1)*this->cp_);
    dhTrans -= dhRet;

    // Correct dhTrans to account for enthalpy of evolved volatiles
    dhTrans +=
        sum(dMassMT)
       *td.cloud().composition().HGas(YGas_, 0.5*(T0 + T1));

    // Correct dhTrans to account for enthalpy of consumed solids
    dhTrans +=
        sum(dMassSR)
       *td.cloud().composition().HSolid(YSolid_, 0.5*(T0 + T1));


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Transfer mass lost from particle to carrier mass source
    forAll(dMassMT, i)
    {
        td.cloud().rhoTrans(i)[celli] += np0*(dMassMT[i] + dMassSR[i]);
    }

    // Update momentum transfer
    td.cloud().UTrans()[celli] += np0*dUTrans;

    // Accumulate coefficient to be applied in carrier phase momentum coupling
    td.cloud().UCoeff()[celli] += np0*mass0*Cud;

    // Update enthalpy transfer
    // - enthalpy of lost solids already accounted for
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
        td.cloud().hTrans()[celli] +=
            np0*mass1
           *(
                YMixture_[0]*td.cloud().composition().HGas(YGas_, T1)
              + YMixture_[2]*td.cloud().composition().HSolid(YSolid_, T1)
            );
        td.cloud().UTrans()[celli] += np0*mass1*U1;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
    {
        this->U_ = U1;
        this->T_ = T1;
        this->cp_ =
            YMixture_[0]*td.cloud().composition().cpGas(YGas_, T1)
          + YMixture_[1]*td.cloud().composition().cpLiquid(YLiquid_, pc_, T1)
          + YMixture_[2]*td.cloud().composition().cpSolid(YSolid_);

        // Update particle density or diameter
        if (td.cloud().massTransfer().changesVolume())
        {
            this->d_ = cbrt(mass1/this->rho_*6.0/mathematicalConstant::pi);
        }
        else
        {
            this->rho_ = mass1/this->volume();
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

    // Mass transfer due to surface reactions from particle to carrier phase
    // - components do not necessarily exist in particle already
    scalarList dMassSR(td.cloud().gases().size(), 0.0);


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


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate mass transfer
    // ~~~~~~~~~~~~~~~~~~~~~~~
    calcMassTransfer(td, dt, T0, T1, dMassMT);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialise enthalpy retention to zero
    scalar dhRet = 0.0;

    calcSurfaceReactions(td, dt, celli, T0, T1, dMassMT, dMassSR, dhRet);

    // New total mass
    const scalar mass1 = mass0 - sum(dMassMT) - sum(dMassSR);

    // New specific heat capacity
    const scalar cp1 =
        YMixture_[0]*td.cloud().composition().cpGas(YGas_, T1)
      + YMixture_[1]*td.cloud().composition().cpLiquid(YLiquid_, pc_, T1)
      + YMixture_[2]*td.cloud().composition().cpSolid(YSolid_);

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
        if (td.cloud().massTransfer().changesVolume())
        {
            this->d_ = cbrt(mass1/this->rho_*6.0/mathematicalConstant::pi);
        }
        else
        {
            this->rho_ = mass1/this->volume();
        }
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcMassTransfer
(
    TrackData& td,
    const scalar dt,
    const scalar T0,
    const scalar T1,
    scalarList& dMassMT
)
{
    if (td.cloud().composition().YMixture0()[1]>SMALL)
    {
        notImplemented
        (
            "void Foam::ReactingParcel<ParcelType>::"
            "calcMassTransfer(...): no treatment currently "
            "available for particles containing liquid species"
        )
    }

    // Check that model is active, and that the parcel temperature is
    // within necessary limits for mass transfer to occur
    if
    (
        !td.cloud().massTransfer().active()
     || this->T_<td.constProps().Tvap()
     || this->T_<td.constProps().Tbp()
    )
    {
        return;
    }

    // Determine mass to add to carrier phase
    const scalar mass = this->mass();
    const scalar dMassTot = td.cloud().massTransfer().calculate
    (
        dt,
        mass0_,
        mass,
        td.cloud().composition().YMixture0(),
        YMixture_,
        T0,
        canCombust_
    );

    // Update (total) mass fractions
    YMixture_[0] = (YMixture_[0]*mass - dMassTot)/(mass - dMassTot);
    YMixture_[1] = YMixture_[1]*mass/(mass - dMassTot);
    YMixture_[2] = 1.0 - YMixture_[0] - YMixture_[1];

    // Add to cummulative mass transfer
    forAll (YGas_, i)
    {
        label id = td.cloud().composition().gasGlobalIds()[i];

        // Mass transfer
        scalar volatileMass = YGas_[i]*dMassTot;
        dMassMT[id] += volatileMass;
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label celli,
    const scalar T0,
    const scalar T1,
    const scalarList& dMassMT,
    scalarList& dMassSR,
    scalar& dhRet
)
{
    // Check that model is active
    if (!td.cloud().surfaceReaction().active() || !canCombust_)
    {
        return;
    }

    // Update mass transfer(s)
    // - Also updates Y()'s
    td.cloud().surfaceReaction().calculate
    (
        dt,
        celli,
        this->d_,
        T0,
        T1,
        this->Tc_,
        pc_,
        this->rhoc_,
        this->mass(),
        dMassMT,
        YGas_,
        YLiquid_,
        YSolid_,
        YMixture_,
        dMassSR,
        dhRet
    );
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //

