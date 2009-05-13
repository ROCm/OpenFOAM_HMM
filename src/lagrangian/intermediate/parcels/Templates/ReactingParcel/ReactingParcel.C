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
    if (pc_ < td.constProps().pMin())
    {
        WarningIn
        (
            "void Foam::ReactingParcel<ParcelType>::updateCellQuantities"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Limiting pressure in cell " << cellI << " to "
            << td.constProps().pMin() <<  nl << endl;

        pc_ = td.constProps().pMin();
    }

    // Apply correction to cell density to account for mass transfer
    scalar addedMass = 0.0;
    forAll(td.cloud().rhoTrans(), i)
    {
        addedMass += td.cloud().rhoTrans(i)[cellI];
    }
    this->rhoc_ += addedMass/td.cloud().pMesh().cellVolumes()[cellI];
}


template<class ParcelType>
void Foam::ReactingParcel<ParcelType>::updateMassFraction
(
    const scalar mass0,
    const scalarField& dMass,
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
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar rho0 = this->rho_;
    const scalar T0 = this->T_;
    const scalar cp0 = this->cp_;
    const scalar mass0 = this->mass();


    // Intialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum
    vector dUTrans = vector::zero;

    // Enthalpy
    scalar dhTrans = 0.0;

    // Mass transfer due to phase change
    scalarField dMassPC(Y_.size(), 0.0);


    // Phase change
    // ~~~~~~~~~~~~

    // Return enthalpy source and calc mass transfer due to phase change
    scalar ShPC =
        calcPhaseChange(td, dt, cellI, d0, T0, U0, 0, 1.0, Y_, dMassPC);

    // Update particle component mass fractions
    updateMassFraction(mass0, dMassPC, Y_);


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
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
            ShPC,
            htc,
            dhTrans
        );


    // Motion
    // ~~~~~~

    // Update mass
    scalar mass1 = mass0 - sum(dMassPC);

    // No additional forces
    vector Fx = vector::zero;

    // Calculate new particle velocity
    scalar Cud = 0.0;
    vector U1 =
        calcVelocity
        (
            td,
            dt,
            cellI,
            d0,
            U0,
            rho0,
            0.5*(mass0 + mass1),
            Fx,
            Cud,
            dUTrans
        );


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().coupled())
    {
        // Transfer mass lost from particle to carrier mass source
        forAll(dMassPC, i)
        {
            label id = td.cloud().composition().localToGlobalGasId(0, i);
            td.cloud().rhoTrans(id)[cellI] += np0*dMassPC[i];
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


    // Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;

        if (td.cloud().coupled())
        {
            // Absorb parcel into carrier phase
            forAll(Y_, i)
            {
                label id = td.cloud().composition().localToGlobalGasId(0, i);
                td.cloud().rhoTrans(id)[cellI] += np0*mass1*Y_[i];
            }
            td.cloud().UTrans()[cellI] += np0*mass1*U1;
            scalar HEff = td.cloud().composition().H(0, Y_, pc_, T1);
            td.cloud().hTrans()[cellI] += np0*mass1*HEff;
        }
    }


    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    else
    {
        this->U_ = U1;
        this->T_ = T1;
        this->cp_ = td.cloud().composition().cp(0, Y_, pc_, T1);

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
    const scalar d,
    const scalar T,
    const vector& U,
    const label idPhase,
    const scalar YPhase,
    const scalarField& YComponents,
    scalarField& dMassPC
)
{
    if
    (
        !td.cloud().phaseChange().active()
     || T < td.constProps().Tvap()
     || YPhase < SMALL
    )
    {
        return 0.0;
    }

    // Calculate mass transfer due to phase change
    td.cloud().phaseChange().calculate
    (
        dt,
        cellI,
        d,
        min(T, td.constProps().Tbp()), // Limit to boiling temperature
        pc_,
        this->Tc_,
        this->muc_/(this->rhoc_ + ROOTVSMALL),
        U - this->Uc_,
        dMassPC
    );

    scalar dMassTot = sum(dMassPC);

    // Add to cumulative phase change mass
    td.cloud().addToMassPhaseChange(this->nParticle_*dMassTot);

    // Effective latent heat of vaporisation
    scalar LEff = td.cloud().composition().L(idPhase, YComponents, pc_, T);

    return -dMassTot*LEff;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //

