/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "KinematicParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    tetIndices tetIs = this->currentTetIndices();

    rhoc_ = td.rhoInterp().interpolate(this->position(), tetIs);

    if (rhoc_ < td.cloud().constProps().rhoMin())
    {
        WarningIn
        (
            "void Foam::KinematicParcel<ParcelType>::setCellValues"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Limiting observed density in cell " << cellI << " to "
            << td.cloud().constProps().rhoMin() <<  nl << endl;

        rhoc_ = td.cloud().constProps().rhoMin();
    }

    Uc_ = td.UInterp().interpolate(this->position(), tetIs);

    muc_ = td.muInterp().interpolate(this->position(), tetIs);

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
void Foam::KinematicParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    Uc_ += td.cloud().UTrans()[cellI]/massCell(cellI);
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

    // Reynolds number
    const scalar Re = this->Re(U0, d0, rhoc_, muc_);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = vector::zero;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    scalar Cud = 0.0;
    vector U1 =
        calcVelocity
        (
            td,
            dt,
            cellI,
            Re,
            muc_,
            d0,
            U0,
            rho0,
            mass0,
            Su,
            dUTrans,
            Cud
        );


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().solution().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Update momentum transfer coefficient
        td.cloud().UCoeff()[cellI] += np0*mass0*Cud;
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
    const scalar Re,
    const scalar mu,
    const scalar d,
    const vector& U,
    const scalar rho,
    const scalar mass,
    const vector& Su,
    vector& dUTrans,
    scalar& Cud
) const
{
    const polyMesh& mesh = this->cloud().pMesh();

    // Momentum transfer coefficient
    const scalar utc = td.cloud().drag().utc(Re, d, mu) + ROOTVSMALL;

    tetIndices tetIs = this->currentTetIndices();

    // Momentum source due to particle forces
    const vector Fcp = mass*td.cloud().forces().calcCoupled
    (
        this->position(),
        tetIs,
        dt,
        rhoc_,
        rho,
        Uc_,
        U,
        d
    );

    const vector Fncp = mass*td.cloud().forces().calcNonCoupled
    (
        this->position(),
        tetIs,
        dt,
        rhoc_,
        rho,
        Uc_,
        U,
        d
    );


    // New particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~

    // Update velocity - treat as 3-D
    const scalar As = this->areaS(d);
    const vector ap = Uc_ + (Fcp + Fncp + Su)/(utc*As);
    const scalar bp = 6.0*utc/(rho*d);

    Cud = bp;

    IntegrationScheme<vector>::integrationResult Ures =
        td.cloud().UIntegrator().integrate(U, dt, ap, bp);

    vector Unew = Ures.value();

    dUTrans += dt*(utc*As*(Ures.average() - Uc_) - Fcp);

    // Apply correction to velocity and dUTrans for reduced-D cases
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);
    meshTools::constrainDirection(mesh, mesh.solutionD(), dUTrans);

    return Unew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const KinematicParcel<ParcelType>& p
)
:
    Particle<ParcelType>(p),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    f_(p.f_),
    angularMomentum_(p.angularMomentum_),
    torque_(p.torque_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    collisionRecords_(p.collisionRecords_),
    rhoc_(p.rhoc_),
    Uc_(p.Uc_),
    muc_(p.muc_)
{}


template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const KinematicParcel<ParcelType>& p,
    const KinematicCloud<ParcelType>& c
)
:
    Particle<ParcelType>(p, c),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    f_(p.f_),
    angularMomentum_(p.angularMomentum_),
    torque_(p.torque_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    collisionRecords_(p.collisionRecords_),
    rhoc_(p.rhoc_),
    Uc_(p.Uc_),
    muc_(p.muc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::KinematicParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    ParcelType& p = static_cast<ParcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    switch (td.part())
    {
        case TrackData::tpVelocityHalfStep:
        {
            // First and last leapfrog velocity adjust part, required
            // before and after tracking and force calculation

            p.U() += 0.5*trackTime*p.f()/p.mass();

            angularMomentum_ += 0.5*trackTime*torque_;

            break;
        }

        case TrackData::tpLinearTrack:
        {
            scalar tEnd = (1.0 - p.stepFraction())*trackTime;
            const scalar dtMax = tEnd;

            while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
            {
                // Apply correction to position for reduced-D cases
                meshTools::constrainToMeshCentre(mesh, p.position());

                // Set the Lagrangian time-step
                scalar dt = min(dtMax, tEnd);

                // Remember which cell the Parcel is in since this
                // will change if a face is hit
                label cellI = p.cell();

                if (p.active())
                {
                    dt *= p.trackToFace(p.position() + dt*U_, td);
                }

                tEnd -= dt;
                p.stepFraction() = 1.0 - tEnd/trackTime;

                // Avoid problems with extremely small timesteps
                if (dt > ROOTVSMALL)
                {
                    // Update cell based properties
                    p.setCellValues(td, dt, cellI);

                    if (td.cloud().solution().cellValueSourceCorrection())
                    {
                        p.cellValueSourceCorrection(td, dt, cellI);
                    }

                    p.calc(td, dt, cellI);
                }

                if (p.onBoundary() && td.keepParticle)
                {
                    if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
                    {
                        td.switchProcessor = true;
                    }
                }

                p.age() += dt;
            }

            break;
        }

        case TrackData::tpRotationalTrack:
        {
            notImplemented("TrackData::tpRotationalTrack");

            break;
        }

        default:
        {
            FatalErrorIn("KinematicParcel<ParcelType>::move(TrackData& td)")
                << td.part() << " is an invalid part of the tracking method."
                << abort(FatalError);
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackData>
bool Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch& pp,
    TrackData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    ParcelType& p = static_cast<ParcelType&>(*this);

    // Invoke post-processing model
    td.cloud().postProcessing().postPatch(p, patchI);

    // Invoke surface film model
    if (td.cloud().surfaceFilm().transferParcel(p, patchI))
    {
        // Parcel transferred to the surface film
        td.keepParticle = false;

        // All interactions done
        return true;
    }
    else
    {
        // Invoke patch interaction model
        return td.cloud().patchInteraction().correct
        (
            static_cast<ParcelType&>(*this),
            pp,
            td.keepParticle,
            trackFraction,
            tetIs
        );
    }
}


template<class ParcelType>
bool Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch& pp,
    int& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    return false;
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
    TrackData& td,
    const tetIndices&
)
{
    // Wall interactions handled by generic hitPatch function
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch&,
    int&,
    const tetIndices&
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
void Foam::KinematicParcel<ParcelType>::transformProperties(const tensor& T)
{
    Particle<ParcelType>::transformProperties(T);

    U_ = transform(T, U_);

    f_ = transform(T, f_);

    angularMomentum_ = transform(T, angularMomentum_);

    torque_ = transform(T, torque_);
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
