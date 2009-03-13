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
    const label cellI
)
{
    ThermoParcel<ParcelType>::updateCellQuantities(td, dt, cellI);

    pc_ = td.pInterp().interpolate(this->position(), cellI);
}


template<class ParcelType>
void Foam::ReactingParcel<ParcelType>::updateMassFraction
(
    const scalar mass0,
    const scalarList& dMass,
    scalarField& Y
)
{
    scalar mass1 = mass0 + sum(dMass);

    forAll(Y, i)
    {
        Y[i] = (Y[i]*mass0 - dMass[i])/mass1;
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar mass0 = this->mass();
    const scalar np0 = this->nParticle_;
    const scalar T0 = this->T_;

    // ~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate velocity
    // ~~~~~~~~~~~~~~~~~~~~~

    scalar Cud = 0.0;
    vector dUTrans = vector::zero;
    const vector U1 = calcVelocity(td, dt, vector::zero, mass0, Cud, dUTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate phase change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Mass transfer from particle to carrier phase
    scalarList dMassPC(td.cloud().gases().size(), 0.0);
    const scalar dMassPCTot = calcPhaseChange(td, dt, cellI, T0, 1.0, dMassPC);

    // Enthalpy change due to change in particle composition (sink)
    scalar ShPC = -dMassPCTot*td.cloud().composition().L(0, Y_, pc, T0);

    // Enthalpy change due to species released into the carrier (source)
    scalar HEff = td.cloud().composition().H(0, Y_, pc, T0);
    ShPC += dMassPCTot*HEff;

    // Update particle component mass fractions
    updateMassFraction(mass0, dMassPC, Y_);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate heat transfer
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    scalar ShHT = 0.0;
    scalar T1 = calcHeatTransfer(td, dt, cellI, ShPC, htc, ShHT);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Collect contributions to determine new particle thermo properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New mass
    scalar mass1 = mass0 - dMassPCTot;

    // Total enthalpy transfer from the particle to the carrier phase
    scalar dhTrans = ShHT + ShPC;

    // New specific heat capacity of mixture - using new temperature
    scalar cp1 = td.cloud().composition().cp(0, Y_, pc_, T1);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 5. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().coupled())
    {
        // Transfer mass lost from particle to carrier mass source
        forAll(dMassPC, i)
        {
            td.cloud().rhoTrans(i)[cellI] += np0*dMassPC[i];
        }

        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Coefficient to be applied in carrier phase momentum coupling
        td.cloud().UCoeff()[cellI] += np0*mass0*Cud;

        // Update enthalpy transfer
        td.cloud().hTrans()[cellI] += np0*dhTrans;

        // Coefficient to be applied in carrier phase enthalpy coupling
        td.cloud().hCoeff()[cellI] += np0*htc*this->areaS();
    }


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 6a. Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;

        if (td.cloud().coupled())
        {
            // Absorb particle(s) into carrier phase
            forAll(Y_, i)
            {
                td.cloud().rhoTrans(i)[cellI] += np0*mass1*Y_[i];
            }
            td.cloud().UTrans()[cellI] += np0*mass1*U1;
            td.cloud().hTrans()[cellI] += np0*mass1*HEff;
        }
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 6b. Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
Foam::scalar Foam::ReactingParcel<ParcelType>::calcPhaseChange
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar T,
    const scalar YPhase, // TODO: NEEDED?????????????????????????????????????????
    scalarList& dMassMT
)
{
    if
    (
        !td.cloud().phaseChange().active()
     || T < td.constProps().Tvap()
     || YPhase > SMALL
    )
    {
        return 0.0;
    }

    scalar dMass = td.cloud().phaseChange().calculate
    (
        dt,
        cellI,
        T,
        pc_,
        this->d_,
        this->Tc_,
        this->muc_/this->rhoc_,
        this->U_ - this->Uc_,
        dMassMT
    );

    td.cloud().addToMassPhaseChange(dMass);

    return dMass;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //

