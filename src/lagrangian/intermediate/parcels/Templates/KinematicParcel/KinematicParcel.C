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

#include "KinematicParcel.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::updateCellQuantities
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    rhoc_ = td.rhoInterp().interpolate(this->position(), cellI);
    if (rhoc_ < SMALL)
    {
        WarningIn
        (
            "void Foam::KinematicParcel<ParcelType>::updateCellQuantities"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Density < " << SMALL << " in cell " << cellI << nl << endl;
    }

    Uc_ = td.UInterp().interpolate(this->position(), cellI);
    muc_ = td.muInterp().interpolate(this->position(), cellI);

    // Apply dispersion components to carrier phase velocity
    Uc_ = td.cloud().dispersion().update
    (
        dt,
        cellI,
        U_,
        Uc_,
        UTurb_,
        tTurb_
    );
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar mass0 = mass();


    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum
    vector dUTrans = vector::zero;


    // Motion
    // ~~~~~~

    // No additional forces
    vector Fx = vector::zero;

    // Calculate new particle velocity
    scalar Cud = 0.0;
    vector U1 = calcVelocity(td, dt, cellI, Fx, mass0, Cud, dUTrans);


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += nParticle_*dUTrans;

        // Coefficient to be applied in carrier phase momentum coupling
        td.cloud().UCoeff()[cellI] += nParticle_*mass0*Cud;
    }


    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    U_ = U1;
}


template<class ParcelType>
template<class TrackData>
const Foam::vector Foam::KinematicParcel<ParcelType>::calcVelocity
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const vector& Fx,
    const scalar mass,
    scalar& Cud,
    vector& dUTrans
) const
{
    // Return linearised term from drag model
    Cud = td.cloud().drag().Cu(U_ - Uc_, d_, rhoc_, rho_, muc_);

    // Initialise total force (per unit mass)
    vector Ftot = vector::zero;

    // Gravity force
    if (td.cloud().forceGravity())
    {
        Ftot += td.g()*(1 - rhoc_/rho_);
    }

    // Virtual mass force
    if (td.cloud().forceVirtualMass())
    {
//        Ftot += td.constProps().Cvm()*rhoc_/rho_*d(Uc - U_)/dt;
    }

    // Pressure gradient force
    if (td.cloud().forcePressureGradient())
    {
        const vector& d = this->mesh().deltaCoeffs()[cellI];
        Ftot += rhoc_/rho_*(U_ & (d^Uc_));
    }


    // New particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~

    // Update velocity - treat as 3-D
    const vector ap = Uc_ + (Ftot + Fx)/(Cud + VSMALL);
    const scalar bp = Cud;

    vector Unew = td.cloud().UIntegrator().integrate(U_, dt, ap, bp).value();

    // Calculate the momentum transfer to the continuous phase
    // - do not include gravity impulse
    dUTrans = -mass*(Unew - U_ - dt*td.g());

    return Unew;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::KinematicParcel<ParcelType>::move
(
    TrackData& td
)
{
    ParcelType& p = static_cast<ParcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    const scalar deltaT = mesh.time().deltaT().value();
    scalar tEnd = (1.0 - p.stepFraction())*deltaT;
    const scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Remember which cell the Parcel is in
        // since this will change if a face is hit
        label cellI = p.cell();

        dt *= p.trackToFace(p.position() + dt*U_, td);

        tEnd -= dt;
        p.stepFraction() = 1.0 - tEnd/deltaT;

        // Update cell based properties
        p.updateCellQuantities(td, dt, cellI);

        // Avoid problems with extremely small timesteps
        if (dt > ROOTVSMALL)
        {
            p.calc(td, dt, cellI);
        }

        if (p.onBoundary() && td.keepParticle)
        {
            if (p.face() > -1)
            {
                if (isType<processorPolyPatch>(pbMesh[p.patch(p.face())]))
                {
                    td.switchProcessor = true;
                }
            }
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td
)
{
    td.cloud().wallInteraction().correct(wpp, this->face(), U_);
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch&,
    int&
)
{}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    int&
)
{}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const tensor& T
)
{
    Particle<ParcelType>::transformProperties(T);
    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    Particle<ParcelType>::transformProperties(separation);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "KinematicParcelIO.C"

// ************************************************************************* //

