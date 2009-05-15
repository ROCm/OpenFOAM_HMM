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
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::LIQ(1);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::SLD(2);


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
      + this->Y_[LIQ]*td.cloud().composition().cp(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().cp(idS, YSolid_, p, T);
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
      + this->Y_[LIQ]*td.cloud().composition().H(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().H(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::LEff
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
        this->Y_[GAS]*td.cloud().composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().L(idS, YSolid_, p, T);
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
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, ROOTVSMALL);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = 1.0 - YMix[GAS] - YMix[LIQ];

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
    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar rho0 = this->rho_;
    const scalar T0 = this->T_;
    const scalar cp0 = this->cp_;
    const scalar mass0 = this->mass();

    const scalar pc = this->pc_;

    const scalarField& YMix = this->Y_;
    const label idG = td.cloud().composition().idGas();
    const label idL = td.cloud().composition().idLiquid();
    const label idS = td.cloud().composition().idSolid();

    // Intial ethalpy state
    scalar H0H = HEff(td, pc, T0, idG, idL, idS);
    scalar H0L = LEff(td, pc, T0, idG, idL, idS);
    scalar H0 = H0H - H0L;


    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);

    // Return enthalpy source and calc mass transfer due to phase change
    scalar ShPC =
        calcPhaseChange
        (
            td,
            dt,
            cellI,
            d0,
            T0,
            U0,
            mass0,
            idL,
            YMix[LIQ],
            YLiquid_,
            dMassPC
        );


    // Devolatilisation
    // ~~~~~~~~~~~~~~~~

    // Mass transfer due to devolatilisation
    scalarField dMassDV(YGas_.size(), 0.0);

    // Return enthalpy source and calc mass transfer due to devolatilisation
    scalar ShDV =
        calcDevolatilisation
        (
            td,
            dt,
            T0,
            mass0,
            this->mass0_,
            idG,
            YMix[GAS],
            YGas_,
            canCombust_,
            dMassDV
        );


    // Surface reactions
    // ~~~~~~~~~~~~~~~~~

    // Change in carrier phase composition due to surface reactions
    scalarField dMassSRGas(YGas_.size(), 0.0);
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField dMassSRCarrier(td.cloud().gases().size(), 0.0);

    // Return enthalpy source and calc mass transfer(s) due to surface reaction
    scalar HReaction =
        calcSurfaceReactions
        (
            td,
            dt,
            cellI,
            d0,
            T0,
            mass0,
            canCombust_,
            dMassPC + dMassDV, // total mass of volatiles evolved
            YMix,
            YGas_,
            YLiquid_,
            YSolid_,
            dMassSRGas,
            dMassSRLiquid,
            dMassSRSolid,
            dMassSRCarrier
        );

    // Heat of reaction split between component retained by particle
    const scalar ShSR = td.constProps().hRetentionCoeff()*HReaction;

    // ...and component added to the carrier phase
    const scalar ShSRc = (1.0 - td.constProps().hRetentionCoeff())*HReaction;


    // Update component mass fractions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField dMassGas = dMassDV + dMassSRGas;
    scalarField dMassLiquid = dMassPC + dMassSRLiquid;
    scalarField dMassSolid = dMassSRSolid;

    scalar mass1 =
        updateMassFractions(mass0, dMassGas, dMassLiquid, dMassSolid);


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Total enthalpy source
    scalar Sh = ShPC + ShDV + ShSR;

    // Calculate new particle temperature
    scalar T1 = calcHeatTransfer(td, dt, cellI, d0, U0, rho0, T0, cp0, Sh);

    // Calculate new enthalpy state based on updated composition at new
    // temperature
    scalar H1H = HEff(td, pc, T1, idG, idL, idS);
    scalar H1L = LEff(td, pc, T1, idG, idL, idS);
    scalar H1 = H1H - H1L;


    // Motion
    // ~~~~~~

    // No additional forces
    vector Fx = vector::zero;

    // Calculate new particle velocity
    vector U1 = calcVelocity(td, dt, cellI, d0, U0, rho0, mass0, Fx);


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (td.cloud().coupled())
    {
        // Transfer mass lost from particle to carrier mass source
        forAll(YGas_, i)
        {
            label id = td.cloud().composition().localToGlobalGasId(GAS, i);
            td.cloud().rhoTrans(id)[cellI] += np0*dMassGas[i];
        }
        forAll(YLiquid_, i)
        {
            label id = td.cloud().composition().localToGlobalGasId(LIQ, i);
            td.cloud().rhoTrans(id)[cellI] += np0*dMassLiquid[i];
        }
//        // No mapping between solid components and carrier phase
//        forAll(YSolid_, i)
//        {
//            label id = td.cloud().composition().localToGlobalGasId(SLD, i);
//            td.cloud().rhoTrans(id)[cellI] += np0*dMassSolid[i];
//        }
        forAll(dMassSRCarrier, i)
        {
            td.cloud().rhoTrans(i)[cellI] += np0*dMassSRCarrier[i];
        }

        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*(mass0*U0 - mass1*U1);

        // Update enthalpy transfer
        td.cloud().hTrans()[cellI] += np0*(mass0*H0 - (mass1*H1 + ShSRc));
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
                label id = td.cloud().composition().localToGlobalGasId(GAS, i);
                td.cloud().rhoTrans(id)[cellI] += np0*mass1*YMix[GAS]*YGas_[i];
            }
            forAll(YLiquid_, i)
            {
                label id =
                    td.cloud().composition().localToGlobalGasId(LIQ, i);
                td.cloud().rhoTrans(id)[cellI] +=
                    np0*mass1*YMix[LIQ]*YLiquid_[i];
            }
//            // No mapping between solid components and carrier phase
//            forAll(YSolid_, i)
//            {
//                label id =
//                    td.cloud().composition().localToGlobalGasId(SLD, i);
//                td.cloud().rhoTrans(id)[cellI] +=
//                    np0*mass1*YMix[SLD]*YSolid_[i];
//            }

            td.cloud().hTrans()[cellI] += np0*mass1*H1;
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
    const scalar mass0,
    const label idVolatile,
    const scalar YVolatileTot,
    const scalarField& YVolatile,
    bool& canCombust,
    scalarField& dMassDV
) const
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

    // Total mass of volatiles evolved
    const scalar dMassTot = td.cloud().devolatilisation().calculate
    (
        dt,
        mass0,
        mass,
        T,
        td.cloud().composition().YMixture0()[idVolatile],
        YVolatileTot,
        canCombust
    );

    // Volatile mass transfer - equal components of each volatile specie
    dMassDV = YVolatile*dMassTot;

    td.cloud().addToMassDevolatilisation(this->nParticle_*dMassTot);

    return -dMassTot*td.constProps().LDevol();
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar d,
    const scalar T,
    const scalar mass,
    const bool canCombust,
    const scalarField& dMassVolatile,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier
) const
{
    // Check that model is active
    if (!td.cloud().surfaceReaction().active() || !canCombust)
    {
        return 0.0;
    }

    // Update surface reactions
    const scalar HReaction = td.cloud().surfaceReaction().calculate
    (
        dt,
        cellI,
        d,
        T,
        this->Tc_,
        this->pc_,
        this->rhoc_,
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        dMassVolatile,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier
    );

    td.cloud().addToMassSurfaceReaction
    (
        this->nParticle_
       *(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    );

    return HReaction;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelIO.C"

// ************************************************************************* //

