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
    const scalar mass0 = this->mass();
    const scalar np0 = this->nParticle_;
    const scalar T0 = this->T_;
    const scalar pc = this->pc_;
    scalarField& YMix = this->Y_;
    const label idL = td.cloud().composition().idLiquid();


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Mass transfer from particle to carrier phase
    scalarList dMassPC(td.cloud().gases().size(), 0.0);
    scalar shPC = 
        calcPhaseChange(td, dt, cellI, T0, idL, YMix[idL], YLiquid_, dMassPC);

    // Update particle component mass fractions
    updateMassFraction(mass0, dMassPC, YLiquid_);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate Devolatilisation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Mass transfer from particle to carrier phase
    // - components exist in particle already
    scalarList dMassDV(td.cloud().gases().size(), 0.0);
    scalar shDV = calcDevolatilisation(td, dt, T0, mass0, idGas, YMix, dMassDV);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Mass transfer of volatile components from particle to carrier phase
    const scalarList dMassMT = dMassPC + dMassDV;
    // Mass transfer due to surface reactions from particle to carrier phase
    // - components do not necessarily exist in particle already
    scalarList dMassSR(td.cloud().gases().size(), 0.0);
    // Initialise enthalpy retention to zero
    scalar dhRet = 0.0;
    calcSurfaceReactions(td, dt, cellI, T0, dMassMT, dMassSR, dhRet);

    // Enthalpy retention divided between particle and carrier phase by the
    // fraction fh and (1 - fh)
    scalar ShSR = td.constProps().fh()*dhRet;
    dhTrans -= (1.0 - td.constProps().fh())*dhRet;

    // Correct dhTrans to account for enthalpy of consumed solids
    dhTrans +=
        sum(dMassSR)*td.cloud().composition().H(idSolid, YSolid_, pc, T0);

    // Correct dhTrans to account for enthalpy of evolved volatiles
    dhTrans +=
        sum(dMassMT)*td.cloud().composition().H(idGas, YGas_, pc, T0);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer
    // ~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;

    // Total enthalpy source
    scalar Sh = ShPC + ShDV + ShSR;

    scalar T1 = calcHeatTransfer(td, dt, cellI, Sh, htc, shHT);


    // ~~~~~~~~~~~~~~~~~~
    // Calculate velocity
    // ~~~~~~~~~~~~~~~~~~
    // Update mass
    scalar mass1 = mass0 - massPC - massD - massSR;
    scalar Cud = 0.0;
    vector dUTrans = vector::zero;
    vector Fx = vector::zero;
    vector U1 = calcVelocity(td, dt, Fx, 0.5*(mass0 + mass1), Cud, dUTrans);



    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Collect contributions to determine new particle thermo properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Update specific heat capacity
    cp1 = cpEff(td, pc, T1);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().coupled())
    {
        // Transfer mass lost from particle to carrier mass source
        forAll(dMassMT, i)
        {
            td.cloud().rhoTrans(i)[cellI] +=
                np0*(dMassPC[i] + dMassDV[i] + dMassSR[i]);
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
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::calcDevolatilisation
(
    TrackData& td,
    const scalar dt,
    const scalar T,
    const scalar mass,
    const label idVolatile,
    scalarField_ YMixture,
    scalarList& dMassMT
)
{
    // Check that model is active, and that the parcel temperature is
    // within necessary limits for devolatilisation to occur
    if
    (
        !td.cloud().devolatilisation().active()
     || T < td.constProps().Tvap()
     || T < td.constProps().Tbp()
    )
    {
        return 0.0;
    }

    // Determine mass to add to carrier phase
    const scalar dMassTot = td.cloud().devolatilisation().calculate
    (
        dt,
        this->mass0_,
        mass,
        td.cloud().composition().YMixture0()[idVolatile],
        YMix[0],
        T,
        canCombust_
    );

    // Update (total) mass fractions
    YMix[0] = (YMix[0]*mass - dMassTot)/(mass - dMassTot);
    YMix[1] = YMix[1]*mass/(mass - dMassTot);
    YMix[2] = 1.0 - YMix[0] - YMix[1];

    // Add to cummulative mass transfer
    forAll (YGas_, i)
    {
        label id = td.cloud().composition().globalIds(idVolatile)[i];

        // Volatiles mass transfer
        scalar volatileMass = YGas_[i]*dMassTot;
        dMassMT[id] += volatileMass;
    }

    td.cloud().addToMassDevolatilisation(dMassTot);

    return = td.constProps().Ldevol()*dMassTot;
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar T,
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
        T,
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

