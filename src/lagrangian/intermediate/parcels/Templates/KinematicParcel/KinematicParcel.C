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
    if (rhoc_ < td.constProps().rhoMin())
    {
        WarningIn
        (
            "void Foam::KinematicParcel<ParcelType>::updateCellQuantities"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Limiting density in cell " << cellI << " to "
            << td.constProps().rhoMin() <<  nl << endl;

        rhoc_ = td.constProps().rhoMin();
    }

    Uc_ = td.UInterp().interpolate(this->position(), cellI);

    // Apply correction to cell velocity to account for momentum transfer
    Uc_ += td.cloud().UTrans()[cellI]/(massCell(cellI));

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
    const scalar np0 = nParticle_;
    const scalar d0 = d_;
    const vector U0 = U_;
    const scalar rho0 = rho_;
    const scalar mass0 = mass();


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
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*mass0*(U0 - U1);
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
    const scalar d,
    const vector& U,
    const scalar rho,
    const scalar mass,
    const vector& Fx
) const
{
    const polyMesh& mesh = this->cloud().pMesh();

    // Return linearised term from drag model
    scalar Cud = td.cloud().drag().Cu(U - Uc_, d, rhoc_, rho, muc_);

    // Calculate particle forces
    vector Ftot = td.cloud().forces().calc(cellI, dt, rhoc_, rho, Uc_, U);


    // New particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~

    // Update velocity - treat as 3-D
    const vector ap = Uc_ + (Ftot + Fx)/(Cud + VSMALL);
    const scalar bp = Cud;

    vector Unew = td.cloud().UIntegrator().integrate(U, dt, ap, bp).value();

    // Apply correction to velocity for reduced-D cases
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);

    return Unew;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::KinematicParcel<ParcelType>::move(TrackData& td)
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
        // Apply correction to position for reduced-D cases
        meshTools::constrainToMeshCentre(mesh, p.position());

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
void Foam::KinematicParcel<ParcelType>::hitPatch(const polyPatch&, int&)
{}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties(const tensor& T)
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

