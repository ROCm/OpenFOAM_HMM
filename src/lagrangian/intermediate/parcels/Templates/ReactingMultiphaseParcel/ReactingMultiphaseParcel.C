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

#include "ReactingMultiphaseParcel.H"

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::cpEff
(
    TrackData& td,
    const scalar p,
    const scalar T
) const
{
    return
        this->Y_[0]*td.cloud().composition().cp(0, YGas_, p, T)
      + this->Y_[1]*td.cloud().composition().cp(0, YLiquid_, p, T)
      + this->Y_[2]*td.cloud().composition().cp(0, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::HEff
(
    TrackData& td,
    const scalar p,
    const scalar T
) const
{
    return
        this->Y_[0]*td.cloud().composition().H(0, YGas_, p, T)
      + this->Y_[1]*td.cloud().composition().H(0, YLiquid_, p, T)
      + this->Y_[2]*td.cloud().composition().H(0, YSolid_, p, T);
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::updateCellQuantities
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ReactingParcel<ParcelType>::updateCellQuantities(td, dt, cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const vector U0 = this->U_;
    const scalar mass0 = this->mass();
    const scalar cp0 = this->cp_;
    const scalar np0 = this->nParticle_;
    scalar T0 = this->T_;
    const scalar pc = this->pc_;
    scalarField& Y = this->Y_;
    scalar mass1 = mass0;

    label idGas = td.cloud().composition().idGas();
    label idLiquid = td.cloud().composition().idLiquid();
    label idSolid = td.cloud().composition().idSolid();

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
    scalar T1 = calcHeatTransfer(td, dt, cellI, htc, dhTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate phase change - update mass, Y, cp, T, dhTrans
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar dMassPC = calcPhaseChange(td, dt, cellI, T1, Y[idLiquid], dMassMT);

    // Update particle mass
    mass1 -= dMassPC;

    // Update particle liquid component mass fractions
    this->updateMassFraction(mass0, dMassMT, YLiquid_);

    // New specific heat capacity of mixture
    scalar cp1 = cpEff(td, pc, T1);

    // Correct temperature due to evaporated components
    // TODO: use hl function in liquidMixture???
//    scalar dhPC = -dMassPCTot*td.cloud().composition().L(0, Y_, pc, T0);
    scalar Lvap = td.cloud().composition().L(idLiquid, this->YLiquid_, pc, T0);
    T1 -= Lvap*dMassPC/(0.5*(mass0 + mass1)*cp1);

    // Correct dhTrans to account for the change in enthalpy due to the
    // liquid phase change
    dhTrans +=
        dMassPC
       *td.cloud().composition().H(idLiquid, YLiquid_, pc, 0.5*(T0 + T1));

//????????????????????????????????????????????????????????????????????????????????
    // Store temperature for the start of the next process
    T0 = T1;
//????????????????????????????????????????????????????????????????????????????????


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate Devolatilisation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar dMassD = calcDevolatilisation(td, dt, T0, dMassMT);

    // Update particle mass
    mass1 -= dMassD;

    // New specific heat capacity of mixture
    cp1 = cpEff(td, pc, T1);

    // Update gas and solid component mass fractions
//    updateMassFraction(mass0, mass1, dMassMT, YGas_);
//    updateMassFraction(mass0, mass1, dMassMT, YSolid_);

    // Correct particle temperature to account for latent heat of
    // devolatilisation
    T1 -=
        td.constProps().Ldevol()
       *sum(dMassMT)
       /(0.5*(mass0 + mass1)*cp1);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialise enthalpy retention to zero
    scalar dhRet = 0.0;

    calcSurfaceReactions(td, dt, cellI, T0, T1, dMassMT, dMassSR, dhRet);

    // New total mass
    mass1 -= sum(dMassSR);


    // Enthalpy retention divided between particle and carrier phase by the
    // fraction fh and (1 - fh)
    T1 += td.constProps().fh()*dhRet/(0.5*(mass0 + mass1)*cp0);
    dhTrans -= (1.0 - td.constProps().fh())*dhRet;

    // Correct dhTrans to account for enthalpy of evolved volatiles
    dhTrans +=
        sum(dMassMT)
       *td.cloud().composition().H(idGas, YGas_, pc, 0.5*(T0 + T1));

    // Correct dhTrans to account for enthalpy of consumed solids
    dhTrans +=
        sum(dMassSR)
       *td.cloud().composition().H(idSolid, YSolid_, pc, 0.5*(T0 + T1));


    if (td.cloud().coupled())
    {
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Accumulate carrier phase source terms
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Transfer mass lost from particle to carrier mass source
        forAll(dMassMT, i)
        {
            td.cloud().rhoTrans(i)[cellI] += np0*(dMassMT[i] + dMassSR[i]);
        }

        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Coefficient to be applied in carrier phase momentum coupling
        td.cloud().UCoeff()[cellI] += np0*mass0*Cud;

        // Update enthalpy transfer
        // - enthalpy of lost solids already accounted for
        td.cloud().hTrans()[cellI] += np0*dhTrans;

        // Coefficient to be applied in carrier phase enthalpy coupling
        td.cloud().hCoeff()[cellI] += np0*htc*this->areaS();
    }


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;

        // Absorb particle(s) into carrier phase
        forAll(dMassMT, i)
        {
            td.cloud().rhoTrans(i)[cellI] += np0*dMassMT[i];
        }
        td.cloud().hTrans()[cellI] += np0*mass1*HEff(td, pc, T1);
        td.cloud().UTrans()[cellI] += np0*mass1*U1;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
    {
        this->U_ = U1;
        this->T_ = T1;
        this->cp_ = cpEff(td, pc, T1);

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
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::calcDevolatilisation
(
    TrackData& td,
    const scalar dt,
    const scalar T,
    scalarList& dMassMT
)
{
    // Check that model is active, and that the parcel temperature is
    // within necessary limits for devolatilisation to occur
    if
    (
        !td.cloud().devolatilisation().active()
     || this->T_<td.constProps().Tvap()
     || this->T_<td.constProps().Tbp()
    )
    {
        return 0.0;
    }

    // Determine mass to add to carrier phase
    const scalar mass = this->mass();
    scalarField& Y = this->Y_;
    const scalar dMassTot = td.cloud().devolatilisation().calculate
    (
        dt,
        this->mass0_,
        mass,
        td.cloud().composition().YMixture0(),
        Y,
        T,
        canCombust_
    );

    // Update (total) mass fractions
    Y[0] = (Y[0]*mass - dMassTot)/(mass - dMassTot);
    Y[1] = Y[1]*mass/(mass - dMassTot);
    Y[2] = 1.0 - Y[0] - Y[1];

    // Add to cummulative mass transfer
    label idGas = td.cloud().composition().idGas();
    forAll (YGas_, i)
    {
        label id = td.cloud().composition().globalIds(idGas)[i];

        // Volatiles mass transfer
        scalar volatileMass = YGas_[i]*dMassTot;
        dMassMT[id] += volatileMass;
    }

    td.cloud().addToMassDevolatilisation(dMassTot);

    return dMassTot;
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label cellI,
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
        cellI,
        this->d_,
        T0,
        T1,
        this->Tc_,
        this->pc_,
        this->rhoc_,
        this->mass(),
        dMassMT,
        YGas_,
        YLiquid_,
        YSolid_,
        this->Y_,
        dMassSR,
        dhRet
    );

    // TODO: td.cloud().addToMassSurfaceReaction(sum(dMassSR));
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelIO.C"

// ************************************************************************* //

