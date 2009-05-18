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
void Foam::ReactingParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ThermoParcel<ParcelType>::setCellValues(td, dt, cellI);

    pc_ = td.pInterp().interpolate(this->position(), cellI);
    if (pc_ < td.constProps().pMin())
    {
        WarningIn
        (
            "void Foam::ReactingParcel<ParcelType>::setCellValues"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Limiting pressure in cell " << cellI << " to "
            << td.constProps().pMin() <<  nl << endl;

        pc_ = td.constProps().pMin();
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    scalar massCell = this->massCell(cellI);

    scalar addedMass = 0.0;
    forAll(td.cloud().rhoTrans(), i)
    {
        addedMass += td.cloud().rhoTrans(i)[cellI];
    }

    this->rhoc_ += addedMass/td.cloud().pMesh().cellVolumes()[cellI];

    scalar massCellNew = massCell + addedMass;
    this->Uc_ += td.cloud().UTrans()[cellI]/massCellNew;

    scalar cpEff = 0;
    if (addedMass > ROOTVSMALL)
    {
        forAll(td.cloud().rhoTrans(), i)
        {
            scalar Y = td.cloud().rhoTrans(i)[cellI]/addedMass;
            cpEff += Y*td.cloud().carrierSpecies()[i].Cp(this->Tc_);
        }
    }
    const scalar cpc = td.cpInterp().psi()[cellI];
    this->cpc_ = (massCell*cpc + addedMass*cpEff)/massCellNew;

    this->Tc_ += td.cloud().hsTrans()[cellI]/(this->cpc_*massCellNew);
}


template<class ParcelType>
Foam::scalar Foam::ReactingParcel<ParcelType>::updateMassFraction
(
    const scalar mass0,
    const scalarField& dMass,
    scalarField& Y
) const
{
    scalar mass1 = mass0 - sum(dMass);

    // only update the mass fractions if the new particle mass is finite
    if (mass1 > ROOTVSMALL)
    {
        forAll(Y, i)
        {
            Y[i] = (Y[i]*mass0 - dMass[i])/mass1;
        }
    }

    return mass1;
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

    // Intial ethalpy state
    scalar H0H = td.cloud().composition().H(0, Y_, pc_, T0);
    scalar H0L = td.cloud().composition().L(0, Y_, pc_, T0);
    scalar H0 = H0H - H0L;


    // Phase change
    // ~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(Y_.size(), 0.0);

    // Return enthalpy source and calc mass transfer due to phase change
    scalar ShPC =
        calcPhaseChange(td, dt, cellI, d0, T0, U0, mass0, 0, 1.0, Y_, dMassPC);

    // Update particle component mass and mass fractions
    scalar mass1 = updateMassFraction(mass0, dMassPC, Y_);


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    scalar T1 = calcHeatTransfer(td, dt, cellI, d0, U0, rho0, T0, cp0, ShPC);

    // Calculate new enthalpy state
    scalar H1H = td.cloud().composition().H(0, Y_, pc_, T1);
    scalar H1L = td.cloud().composition().L(0, Y_, pc_, T1);
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
        forAll(dMassPC, i)
        {
            label id = td.cloud().composition().localToGlobalGasId(0, i);
            td.cloud().rhoTrans(id)[cellI] += np0*dMassPC[i];
        }

        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*(mass0*U0 - mass1*U1);

        // Update sensible enthalpy transfer
        td.cloud().hsTrans()[cellI] += np0*(mass0*H0 - mass1*H1);
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
            td.cloud().hsTrans()[cellI] += np0*mass1*H1;
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
    const scalar mass,
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

    // Limit phase change mass by availability of each specie
    dMassPC = min(mass*YPhase*YComponents, dMassPC);

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

