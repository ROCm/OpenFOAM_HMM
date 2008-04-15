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

#include "ReactingParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackingData>
void Foam::ReactingParcel<ParcelType>::calcCoupled
(
    TrackingData& td,
    const label celli,
    const scalar dt,
    const scalar rhoc,
    vector& Uc,
    const scalar muc,
    const scalar Tc,
    const scalar cpc,
    const scalar pc
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const vector U0 = this->U();
    const scalar T0 = this->T();
    const scalar mass0 = this->mass();
    const scalar cp0 = this->cp();
    const scalar np0 = this->nParticle();

    // Mass transfer from particle to carrier phase
    // - components exist in particle already
    scalarList dMassMT(td.cloud().gases().size(), 0.0);

    // Mass transfer due to surface reactions from particle to carrier phase
    // - components do not necessarily exist in particle already
    scalarList dMassSR(td.cloud().gases().size(), 0.0);

    // Total mass lost from particle due to surface reactions
    // - sub-model will adjust component mass fractions
    scalar dMassMTSR = 0.0;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    scalar T1 = calcHeatTransfer(td, celli, dt, rhoc, Uc, muc, Tc, cpc, htc);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, rhoc, Uc, muc, Cud);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate mass transfer
    // ~~~~~~~~~~~~~~~~~~~~~~~
    calcMassTransfer(td, dt, T0, T1, dMassMT);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    calcSurfaceReactions
    (
        td,
        dt,
        celli,
        rhoc,
        Tc,
        T0,
        T1,
        dMassMTSR,
        dMassSR
    );

    // New total mass
    const scalar mass1 = mass0 - sum(dMassMT) - dMassMTSR;

    // Ratio of mass devolatilised to the total volatile mass of the particle
    const scalar fVol = 1 -
        (YMixture_[0]*mass1)
       /(td.cloud().composition().YMixture0()[0]*mass0_);

    // Specific heat capacity of non-volatile components
    const scalar cpNonVolatile =
        (
            YMixture_[1]*td.cloud().composition().cpLiquid(YLiquid_, pc, Tc)
          + YMixture_[2]*td.cloud().composition().cpSolid(YSolid_)
        )/(YMixture_[1] + YMixture_[2]);

    // New specific heat capacity - linear variation until volatiles
    // have evolved
    const scalar cp1 = (cpNonVolatile - td.constProps().cp0())*fVol
       + td.constProps().cp0();


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Transfer mass lost from particle to carrier mass source
    forAll(dMassMT, i)
    {
        td.cloud().rhoTrans(i)[celli] +=
            np0*(dMassMT[i] + dMassSR[i]);
    }

    // Update momentum transfer
    td.cloud().UTrans()[celli] += np0*(mass0*U0 - mass1*U1);

    // Accumulate coefficient to be applied in carrier phase momentum coupling
    td.cloud().UCoeff()[celli] += np0*mass0*Cud;

    // Update enthalpy transfer
    td.cloud().hTrans()[celli] += np0*(mass0*cp0*T0 - mass1*cp1*T1);

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
        td.cloud().hTrans()[celli] += np0*mass1*cp1*T1;
        td.cloud().UTrans()[celli] += np0*mass1*U1;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
    {
        this->U() = U1;
        this->T() = T1;
        this->cp() = cp1;

        // Update particle density or diameter
        if (td.cloud().massTransfer().changesVolume())
        {
            this->d() = cbrt(mass1/this->rho()*6.0/mathematicalConstant::pi);
        }
        else
        {
            this->rho() = mass1/this->volume();
        }
    }
}


template<class ParcelType>
template<class TrackingData>
void Foam::ReactingParcel<ParcelType>::calcUncoupled
(
    TrackingData& td,
    const label celli,
    const scalar dt,
    const scalar rhoc,
    vector& Uc,
    const scalar muc,
    const scalar Tc,
    const scalar cpc,
    const scalar pc
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar T0 = this->T();
    const scalar mass0 = this->mass();
//    const scalar cp0 = this->cp();

    // Mass transfer from particle to carrier phase
    // - components exist in particle already
    scalarList dMassMT(td.cloud().gases().size(), 0.0);

    // Mass transfer due to surface reactions from particle to carrier phase
    // - components do not necessarily exist in particle already
    scalarList dMassSR(td.cloud().gases().size(), 0.0);

    // Total mass lost from particle due to surface reactions
    // - sub-model will adjust component mass fractions
    scalar dMassMTSR = 0.0;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    scalar T1 = calcHeatTransfer(td, celli, dt, rhoc, Uc, muc, Tc, cpc, htc);

    // Limit new temp max by vapourisarion temperature
    T1 = min(td.constProps().Tvap(), T1);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Fd = 0.0;
    const vector U1 = calcVelocity(td, dt, rhoc, Uc, muc, Fd);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate mass transfer
    // ~~~~~~~~~~~~~~~~~~~~~~~
    calcMassTransfer(td, dt, T0, T1, dMassMT);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    calcSurfaceReactions
    (
        td,
        dt,
        celli,
        rhoc,
        Tc,
        T0,
        T1,
        dMassMTSR,
        dMassSR
    );

    // New total mass
    const scalar mass1 = mass0 - sum(dMassMT) - dMassMTSR;

    // Ratio of mass devolatilised to the total volatile mass of the particle
    const scalar fVol = 1 -
        (YMixture_[0]*mass1)
       /(td.cloud().composition().YMixture0()[0]*mass0_);

    // Specific heat capacity of non-volatile components
    const scalar cpNonVolatile =
        (
            YMixture_[1]*td.cloud().composition().cpLiquid(YLiquid_, pc, Tc)
          + YMixture_[2]*td.cloud().composition().cpSolid(YSolid_)
        )/(YMixture_[1] + YMixture_[2]);

    // New specific heat capacity - linear variation until volatiles
    // have evolved
    const scalar cp1 = (cpNonVolatile - td.constProps().cp0())*fVol
       + td.constProps().cp0();


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
        this->U() = U1;
        this->T() = T1;
        this->cp() = cp1;

        // Update particle density or diameter
        if (td.cloud().massTransfer().changesVolume())
        {
            this->d() = cbrt(mass1/this->rho()*6.0/mathematicalConstant::pi);
        }
        else
        {
            this->rho() = mass1/this->volume();
        }
    }
}


template<class ParcelType>
template<class TrackingData>
void Foam::ReactingParcel<ParcelType>::calcMassTransfer
(
    TrackingData& td,
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
     || this->T()<td.constProps().Tvap()
     || this->T()<td.constProps().Tbp()
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
template<class TrackingData>
void Foam::ReactingParcel<ParcelType>::calcSurfaceReactions
(
    TrackingData& td,
    const scalar dt,
    const label celli,
    const scalar rhoc,
    const scalar Tc,
    const scalar T0,
    const scalar T1,
    scalar& dMassMTSR,
    scalarList& dMassMT
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
        this->d(),
        T0,
        T1,
        Tc,
        rhoc,
        this->mass(),
        YGas_,
        YLiquid_,
        YSolid_,
        YMixture_,
        dMassMTSR,
        dMassMT
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackingData>
bool Foam::ReactingParcel<ParcelType>::move
(
    TrackingData& td
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    const scalar deltaT = mesh.time().deltaT().value();
    scalar tEnd = (1.0 - this->stepFraction())*deltaT;
    const scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Remember which cell the parcel is in
        // since this will change if a face is hit
        label celli = this->cell();

        dt *= trackToFace(this->position() + dt*this->U(), td);

        tEnd -= dt;

        this->stepFraction() = 1.0 - tEnd/deltaT;

        // Avoid div0 in reacting sub-models
        if (dt < SMALL)
        {
            break;
        }

        cellPointWeight cpw
        (
            mesh,
            this->position(),
            celli,
            this->faceInterpolation()
        );
        scalar rhoc = td.rhoInterp().interpolate(cpw);
        vector Uc = td.UInterp().interpolate(cpw);
        scalar muc = td.muInterp().interpolate(cpw);
        scalar Tc = td.TInterp().interpolate(cpw);
        scalar cpc = td.cpInterp().interpolate(cpw);
        scalar pc = td.pInterp().interpolate(cpw);

        Uc = td.cloud().dispersion().update
        (
            dt,
            celli,
            this->U(),
            Uc,
            this->UTurb(),
            this->tTurb()
        );

        if (td.cloud().coupled())
        {
            calcCoupled(td, celli, dt, rhoc, Uc, muc, Tc, cpc, pc);
        }
        else
        {
            calcUncoupled(td, celli, dt, rhoc, Uc, muc, Tc, cpc, pc);
        }

        if (this->onBoundary() && td.keepParticle)
        {
            if (this->face() > -1)
            {
                if (isType<processorPolyPatch>(pbMesh[this->patch(this->face())]))
                {
                    td.switchProcessor = true;
                }
            }
        }
    }

    return td.keepParticle;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //

