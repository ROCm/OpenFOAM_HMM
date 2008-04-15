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

#include "ThermoParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackingData>
void Foam::ThermoParcel<ParcelType>::calcCoupled
(
    TrackingData& td,
    const label celli,
    const scalar dt,
    const scalar rhoc,
    vector& Uc,
    const scalar muc,
    const scalar Tc,
    const scalar cpc
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


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, rhoc, Uc, muc, Cud);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    const scalar T1 =
        calcHeatTransfer(td, celli, dt, rhoc, Uc, muc, Tc, cpc, htc);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Update momentum transfer
    td.cloud().UTrans()[celli] += np0*(mass0*(U0 - U1));

    // Accumulate coefficient to be applied in carrier phase momentum coupling
    td.cloud().UCoeff()[celli] += np0*mass0*Cud;

    // Update enthalpy transfer
    td.cloud().hTrans()[celli] += np0*mass0*cp0*(T0 - T1);

    // Accumulate coefficient to be applied in carrier phase enthalpy coupling
    td.cloud().hCoeff()[celli] += np0*htc*this->areaS();


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this->U() = U1;
    this->T() = T1;
}


template<class ParcelType>
template<class TrackingData>
void Foam::ThermoParcel<ParcelType>::calcUncoupled
(
    TrackingData& td,
    const label celli,
    const scalar dt,
    const scalar rhoc,
    vector& Uc,
    const scalar muc,
    const scalar Tc,
    const scalar cpc
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    this->U() = calcVelocity(td, dt, rhoc, Uc, muc, Cud);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    this->T() = calcHeatTransfer(td, celli, dt, rhoc, Uc, muc, Tc, cpc, htc);
}


template<class ParcelType>
template<class TrackingData>
Foam::scalar Foam::ThermoParcel<ParcelType>::calcHeatTransfer
(
    TrackingData& td,
    const label celli,
    const scalar dt,
    const scalar rhoc,
    const vector& Uc,
    const scalar muc,
    const scalar Tc,
    const scalar cpc,
    scalar& htc
)
{
    if (!td.cloud().heatTransfer().active())
    {
        return T_;
    }

    // Calc heat transfer coefficient
    htc = td.cloud().heatTransfer().h
    (
        this->d(),
        this->Ur(),
        rhoc,
        this->rho(),
        cpc,
        cp_,
        muc
    );

    // Determine ap and bp coefficients
    scalar ap = Tc;
    scalar bp = htc;
    if (td.cloud().radiation())
    {
        // Carrier phase incident radiation field
        // Currently the G field is not interpolated to the parcel position
        // - instead, the cell centre value is applied directly
        const scalarField& G = td.cloud().mesh().objectRegistry
            ::lookupObject<volScalarField>("G");

        // Helper variables
        const scalar sigma = radiation::sigmaSB.value();
        const scalar epsilon = td.constProps().epsilon0();
        const scalar epsilonSigmaT3 = epsilon*sigma*pow3(T_);
        ap = (htc*Tc + 0.25*epsilon*G[celli])/(htc + epsilonSigmaT3);
        bp += epsilonSigmaT3;
    }
    bp *= 6.0/(this->rho()*this->d()*cp_);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle temperature
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Analytical
    const scalar Tnew = ap + (T_ - ap)*exp(-bp*dt);

    // Euler-implicit
//    const scalar Tnew = (T_ + dt*ap*bp)/(1.0 + dt*bp);

    return Tnew;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackingData>
bool Foam::ThermoParcel<ParcelType>::move
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
            calcCoupled(td, celli, dt, rhoc, Uc, muc, Tc, cpc);
        }
        else
        {
            calcUncoupled(td, celli, dt, rhoc, Uc, muc, Tc, cpc);
        }

        if (this->onBoundary() && td.keepParticle)
        {
            if (this->face() > -1)
            {
                if
                (
                    isType<processorPolyPatch>
                        (pbMesh[this->patch(this->face())])
                )
                {
                    td.switchProcessor = true;
                }
            }
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackingData>
void Foam::ThermoParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackingData& td
)
{
    td.cloud().wallInteraction().correct(wpp, this->face(), this->U());
}



template<class ParcelType>
void Foam::ThermoParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    int&
)
{}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ThermoParcelIO.C"

// ************************************************************************* //

