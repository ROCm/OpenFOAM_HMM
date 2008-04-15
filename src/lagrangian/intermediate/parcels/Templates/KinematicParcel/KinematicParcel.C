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

#include "KinematicParcel.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackingData>
void Foam::KinematicParcel<ParcelType>::calcCoupled
(
    TrackingData& td,
    const label celli,
    const scalar dt,
    const scalar rhoc,
    vector& Uc,
    const scalar muc
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar mass0 = mass();
    const vector U0 = U_;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, rhoc, Uc, muc, Cud);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Update momentum transfer
    td.cloud().UTrans()[celli] += nParticle_*mass0*(U0 - U1);

    // Accumulate coefficient to be applied in carrier phase momentum coupling
    td.cloud().UCoeff()[celli] += nParticle_*mass0*Cud;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this->U() = U1;
}


template<class ParcelType>
template<class TrackingData>
void Foam::KinematicParcel<ParcelType>::calcUncoupled
(
    TrackingData& td,
    const scalar dt,
    const scalar rhoc,
    vector& Uc,
    const scalar muc
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    this->U() = calcVelocity(td, dt, rhoc, Uc, muc, Cud);
}


template<class ParcelType>
template<class TrackingData>
Foam::vector Foam::KinematicParcel<ParcelType>::calcVelocity
(
    TrackingData& td,
    const scalar dt,
    const scalar rhoc,
    vector& Uc,
    const scalar muc,
    scalar& Cud
)
{
    // Correct carrier phase velocity for 2-D slab cases
    const polyMeshInfo& meshInfo = td.cloud().meshInfo();
    if (meshInfo.caseIs2dSlab())
    {
        Uc.component(meshInfo.emptyComponent()) = 0.0;
    }

    // Update relative velocity
    Ur_ = U_ - Uc;

    // Return linearised term from drag model
//    const scalar Cud = td.cloud().drag().Cu
    Cud = td.cloud().drag().Cu
    (
        Ur_,
        d_,
        rhoc,
        rho_,
        muc
    );

    // Update velocity - treat as 3-D
    const vector ap = (1.0 - rhoc/rho_)*td.g();
    const scalar bp = 1.0/Cud;


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Analytical
//    const scalar expTerm = exp(-dt/bp);
//    vector Unew = Uc + (U_ - Uc)*expTerm + ap*bp*(1.0 - expTerm);

    // Euler-implicit
    vector Unew = (U_ + dt*(ap + Uc/bp))/(1.0 + dt/bp);

    // Make corrections for 2-D cases
    if (meshInfo.caseIs2d())
    {
        if (meshInfo.caseIs2dSlab())
        {
            // Remove the slab normal parcel velocity component
            Unew.component(meshInfo.emptyComponent()) = 0.0;

            // Snap parcels to central plane
            this->position().component(meshInfo.emptyComponent()) =
                meshInfo.centrePoint().component(meshInfo.emptyComponent());
        }
        else if (meshInfo.caseIs2dWedge())
        {
            // Snap parcels to central plane
            this->position().component(meshInfo.emptyComponent()) = 0.0;
        }
        else
        {
            FatalErrorIn("void Foam::KinematicParcel::calcVelocity")
                << "Could not determine 2-D case geometry" << nl
                << abort(FatalError);
        }
    }

    return Unew;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ParcelType>
template <class TrackingData>
bool Foam::KinematicParcel<ParcelType>::move
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

        // Remember which cell the Parcel is in
        // since this will change if a face is hit
        label celli = this->cell();

        dt *= trackToFace(this->position() + dt*U_, td);

        tEnd -= dt;
        this->stepFraction() = 1.0 - tEnd/deltaT;

        cellPointWeight cpw
        (
            mesh,
            this->position(),
            celli,
            faceInterpolation()
        );
        scalar rhoc = td.rhoInterp().interpolate(cpw);
        vector Uc = td.UInterp().interpolate(cpw);
        scalar muc = td.muInterp().interpolate(cpw);

        Uc = td.cloud().dispersion().update
        (
            dt,
            celli,
            U_,
            Uc,
            UTurb_,
            tTurb_
        );

        if (td.cloud().coupled())
        {
            calcCoupled(td, celli, dt, rhoc, Uc, muc);
        }
        else
        {
            calcUncoupled(td, dt, rhoc, Uc, muc);
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


template <class ParcelType>
template <class TrackingData>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackingData& td
)
{
    td.switchProcessor = true;
}


template <class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


template <class ParcelType>
template <class TrackingData>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackingData& td
)
{
    td.cloud().wallInteraction().correct(wpp, this->face(), U_);
}


template <class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch&,
    int&
)
{}


template <class ParcelType>
template <class TrackingData>
void Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackingData& td
)
{
    td.keepParticle = false;
}


template <class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    int&
)
{}


template <class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const tensor& T
)
{
    Particle<ParcelType>::transformProperties(T);
    U_ = transform(T, U_);
}


template <class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    Particle<ParcelType>::transformProperties(separation);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "KinematicParcelIO.C"

// ************************************************************************* //

