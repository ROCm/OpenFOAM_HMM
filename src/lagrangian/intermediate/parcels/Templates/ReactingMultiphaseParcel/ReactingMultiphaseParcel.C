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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::GAS(0);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::LIQUID(1);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::SOLID(2);


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::cpEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().cp(idG, YGas_, p, T)
      + this->Y_[LIQUID]*td.cloud().composition().cp(idL, YLiquid_, p, T)
      + this->Y_[SOLID]*td.cloud().composition().cp(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::HEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().H(idG, YGas_, p, T)
      + this->Y_[LIQUID]*td.cloud().composition().H(idL, YLiquid_, p, T)
      + this->Y_[SOLID]*td.cloud().composition().H(idS, YSolid_, p, T);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarList& YMix = this->Y_;

    scalar dMassGasTot = sum(dMassGas);
    scalar dMassLiquidTot = sum(dMassLiquid);
    scalar dMassSolidTot = sum(dMassSolid);

    this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    this->updateMassFraction(mass0*YMix[LIQUID], dMassLiquid, YLiquid_);
    this->updateMassFraction(mass0*YMix[SOLID], dMassSolid, YSolid_);

    scalar massNew = mass0 - (dMassGasTot + dMassLiquidTot + dMassSolidTot);

    YMix[GAS] = (mass0*YMix[GAS] - dMassGas)/massNew;
    YMix[LIQUID] = (mass0*YMix[LIQUID] - dMassLiquid)/massNew;
    YMix[SOLID] = 1.0 - YMix[GAS] - YMix[LIQUID];

    return massNew;
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
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar mass0 = this->mass();
    const scalar np0 = this->nParticle_;
    const scalar T0 = this->T_;
    const scalar pc = this->pc_;
    scalarField& YMix = this->Y_;
    const label idG = td.cloud().composition().idGas();
    const label idL = td.cloud().composition().idLiquid();
    const label idS = td.cloud().composition().idSolid();


    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum
    vector dUTrans = vector::zero;

    // Enthalpy
    scalar dhTrans = 0.0;

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);

    // Mass transfer due to devolatilisation
    scalarField dMassDV(YGas_.size(), 0.0);

    // Change in carrier phase composition due to surface reactions
    scalarField dMassSRGas(YGas_.size(), 0.0);
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField dMassSRCarrier(td.cloud().gases().size(), 0.0);


    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Return enthalpy source and calc mass transfer due to phase change
    scalar ShPC =
        calcPhaseChange(td, dt, cellI, T0, idL, YMix[idL], YLiquid_, dMassPC);


    // Devolatilisation
    // ~~~~~~~~~~~~~~~~

    // Return enthalpy source and calc mass transfer due to devolatilisation
    scalar ShDV = calcDevolatilisation(td, dt, T0, mass0, idG, YMix, dMassDV);


    // Surface reactions
    // ~~~~~~~~~~~~~~~~~

    // Return enthalpy source and calc mass transfer(s) due to surface reaction
    scalar ShSR =
        calcSurfaceReactions
        (
            td,
            dt,
            cellI,
            T0,
            dMassPC + dMassDV, // total volatiles evolved
            dMassSRGas,
            dMassSRLiquid,
            dMassSRSolid,
            dMassSRCarrier,
            dhTrans
        );


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Total enthalpy source
    scalar Sh = ShPC + ShDV + ShSR;

    // Calculate new particle temperature
    scalar htc = 0.0;
    scalar T1 = calcHeatTransfer(td, dt, cellI, Sh, htc, dhTrans);


    // Update component mass fractions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField dMassGas = dMassDV + dMassSRGas;
    scalarField dMassLiquid = dMassPC + dMassSRLiquid;
    scalarField dMassSolid = dMassSRSolid;

    scalar mass1 =
        updateMassFractions(mass0, dMassGas, dMassLiquid, dMassSolid);


    // Motion
    // ~~~~~~

    // No additional forces
    vector Fx = vector::zero;

    // Calculate new particle velocity
    scalar Cud = 0.0;
    vector U1 =
        calcVelocity(td, dt, cellI, Fx, 0.5*(mass0 + mass1), Cud, dUTrans);


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (td.cloud().coupled())
    {
        // Transfer mass lost from particle to carrier mass source
        forAll(YGas_, i)
        {
            label id = td.composition().localToGlobalGasId(GAS, i);
            td.cloud().rhoTrans(id)[cellI] += np0*dMassGas[i];
        }
        forAll(YLiquid_, i)
        {
            label id = td.composition().localToGlobalGasId(LIQUID, i);
            td.cloud().rhoTrans(id)[cellI] += np0*dMassLiquid[i];
        }
        // No mapping between solid components and carrier phase
//        forAll(YSolid_, i)
//        {
//            label id = td.composition().localToGlobalGasId(SOLID, i);
//            td.cloud().rhoTrans(id)[cellI] += np0*dMassSolid[i];
//        }
        forAll(dMassSRCarrier, i)
        {
            td.cloud().rhoTrans(i)[cellI] += np0*dMassSRCarrier[i];
        }

        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Coefficient to be applied in carrier phase momentum coupling
        td.cloud().UCoeff()[cellI] += np0*mass0*Cud;

        // Update enthalpy transfer
        // - enthalpy of lost solids already accounted for
        td.cloud().hTrans()[cellI] += np0*(dhTrans + Sh);

        // Coefficient to be applied in carrier phase enthalpy coupling
        td.cloud().hCoeff()[cellI] += np0*htc*this->areaS();
    }


    // Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;

        if (td.cloud().coupled())
        {
            // Absorb parcel into carrier phase
            forAll(YGas_, i)
            {
                label id = td.composition().localToGlobalGasId(GAS, i);
                td.cloud().rhoTrans(id)[cellI] += np0*mass1*YMix[GAS]*YGas_[i];
            }
            forAll(YLiquid_, i)
            {
                label id = td.composition().localToGlobalGasId(LIQUID, i);
                td.cloud().rhoTrans(id)[cellI] +=
                    np0
                   *mass1
                   *YMix[LIQUID]
                   *YLiquid_[i];
            }
            // No mapping between solid components and carrier phase
//            forAll(YSolid_, i)
//            {
//                label id = td.composition().localToGlobalGasId(SOLID, i);
//                td.cloud().rhoTrans(id)[cellI] +=
//                    np0
//                   *mass1
//                   *YMix[SOLID]
//                   *YSolid_[i];
//            }

            td.cloud().hTrans()[cellI] +=
                np0
               *mass1
               *HEff(td, pc, T1, idG, idL, idS);
            td.cloud().UTrans()[cellI] += np0*mass1*U1;
        }
    }


    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    else
    {
        this->U_ = U1;
        this->T_ = T1;
        this->cp_ = cpEff(td, pc, T1, idG, idL, idS);

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
    const scalarField& YMixture,
    scalarList& dMass
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
        T,
        td.cloud().composition().YMixture0()[idVolatile],
        YMixture[GAS],
        canCombust_
    );

    // Add to cummulative mass transfer
    forAll (YGas_, i)
    {
        label id = td.cloud().composition().globalIds(idVolatile)[i];

        // Volatiles mass transfer
        scalar volatileMass = YGas_[i]*dMassTot;
        dMass[i] += volatileMass;
    }

    td.cloud().addToMassDevolatilisation(dMassTot);

    return td.constProps().Ldevol()*dMassTot;
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar T,
    const scalarList& dMassVolatile,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarList& dMassSRCarrier,
    scalar& dhTrans
)
{
    // Check that model is active
    if (!td.cloud().surfaceReaction().active() || !canCombust_)
    {
        return;
    }

    // Update surface reactions
    scalar HReaction = td.cloud().surfaceReaction().calculate
    (
        dt,
        cellI,
        this->d_,
        T,
        this->Tc_,
        this->pc_,
        this->rhoc_,
        this->mass(),
        YGas_,
        YLiquid_,
        YSolid_,
        this->Y_,
        dMassVolatile,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier
    );

    // Heat of reaction divided between particle and carrier phase by the
    // fraction fh and (1 - fh)
    dhTrans -= (1.0 - td.constProps().fh())*HReaction;

/*
    // Correct dhTrans to account for enthalpy of consumed solids
    dhTrans +=
        sum(dMassSR)*td.cloud().composition().H(idSolid, YSolid_, pc, T0);

    // Correct dhTrans to account for enthalpy of evolved volatiles
    dhTrans +=
        sum(dMassMT)*td.cloud().composition().H(idGas, YGas_, pc, T0);

    // TODO: td.cloud().addToMassSurfaceReaction(sum(dMassSR));
*/

    return td.constProps().fh()*HReaction;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelIO.C"

// ************************************************************************* //

