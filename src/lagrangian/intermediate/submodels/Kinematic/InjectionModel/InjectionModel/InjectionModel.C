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

#include "InjectionModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectionModel<CloudType>::prepareForNextTimeStep
(
    const scalar time0,
    const scalar time1,
    label& nParcels,
    scalar& volume
)
{
    // Initialise values
    nParcels = 0;
    volume = 0.0;

    // Return if not started injection event
    if (time1 < SOI_)
    {
        timeStep0_ = time1;
        return;
    }

    // Make times relative to SOI
    scalar t0 = timeStep0_ - SOI_;
    scalar t1 = time1 - SOI_;

    // Number of parcels to inject
    nParcels = nParcelsToInject(t0, t1);

    // Volume of parcels to inject
    volume = volumeToInject(t0, t1);

    // Hold previous time if no parcels, but non-zero volume fraction
    if ((nParcels == 0) && (volume > 0.0))
    {
        // hold value of timeStep0_
    }
    else
    {
        // advance value of timeStep0_
        timeStep0_ = time1;
    }
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::findInjectorCellAndPosition
(
    label& cellI,
    vector& position
)
{
    const vector p0 = position;

    bool foundCell = false;

    cellI = owner_.mesh().findCell(position);

    if (cellI >= 0)
    {
        const vector& C = owner_.mesh().C()[cellI];
        position += 1.0e-6*(C - position);

        foundCell = owner_.mesh().pointInCell(position, cellI);
    }
    reduce(foundCell, orOp<bool>());

    // Last chance - find nearest cell and try that one
    // - the point is probably on an edge
    if (!foundCell)
    {
        cellI = owner_.mesh().findNearestCell(position);

        if (cellI >= 0)
        {
            const vector& C = owner_.mesh().C()[cellI];
            position += 1.0e-6*(C - position);

            foundCell = owner_.mesh().pointInCell(position, cellI);
        }
        reduce(foundCell, orOp<bool>());
    }

    if (!foundCell)
    {
        FatalErrorIn
        (
            "InjectionModel<CloudType>::setInjectorCellAndPosition"
            "(label&, vector&)"
        )<< "Cannot find parcel injection cell. "
         << "Parcel position = " << p0 << nl
         << abort(FatalError);
    }
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::setNumberOfParticles
(
    const label nParcels,
    const scalar diameter,
    const scalar volumeFraction,
    const scalar rho,
    const scalar volume
)
{
    scalar nP = 0.0;
    switch (parcelBasis_)
    {
        case pbMass:
        {
            nP = volumeFraction*massTotal_/nParcels
               /(rho*mathematicalConstant::pi/6.0*pow3(diameter));
            break;
        }
        case pbNumber:
        {
            nP = volumeFraction*massTotal_/(rho*volume);
            break;
        }
        default:
        {
            nP = 0.0;
            FatalErrorIn
            (
                "void Foam::InjectionModel<CloudType>::setNumberOfParticles"
                "(const label, const scalar, const scalar, const scalar, "
                "const scalar)"
            )<< "Unknown parcelBasis type" << nl
             << exit(FatalError);
        }
    }

    return nP;
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::postInjectCheck()
{
    if (nParcelsAdded_ > 0)
    {
        Pout<< "\n--> Cloud: " << owner_.name() << nl
            << "    Added " << nParcelsAdded_
            <<  " new parcels" << nl << endl;
    }

    // Increment total number of parcels added
    nParcelsAddedTotal_ += nParcelsAdded_;

    // Reset parcel counters
    nParcelsAdded_ = 0;

    // Update time for start of next injection
    time0_ = owner_.db().time().value();

    // Increment number of injections
    nInjections_++;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:   dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    SOI_(readScalar(coeffDict_.lookup("SOI"))),
    volumeTotal_(0.0),
    massTotal_(dimensionedScalar(coeffDict_.lookup("massTotal")).value()),
    massInjected_(0.0),
    nInjections_(0),
    nParcelsAdded_(0),
    nParcelsAddedTotal_(0),
    parcelBasisType_(coeffDict_.lookup("parcelBasisType")),
    parcelBasis_(pbNumber),
    time0_(owner.db().time().value()),
    timeStep0_(0.0)
{
    if (parcelBasisType_ == "mass")
    {
        parcelBasis_ = pbMass;
    }
    else if (parcelBasisType_ == "number")
    {
        parcelBasis_ = pbNumber;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::InjectionModel<CloudType>::InjectionModel"
            "(const dictionary&, CloudType&, const word&)"
        )<< "parcelBasisType must be either 'number' or 'mass'" << nl
         << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::~InjectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class TrackData>
void Foam::InjectionModel<CloudType>::inject(TrackData& td)
{
    const scalar time = owner_.db().time().value();
    const scalar continuousDt = owner_.db().time().deltaT().value();

    // Prepare for next time step
    nParcelsAdded_ = 0;
    label nParcels = 0;
    scalar volume = 0.0;
    prepareForNextTimeStep(time0_, time, nParcels, volume);

    // Return if no parcels are required
    if (nParcels == 0)
    {
        postInjectCheck();
        return;
    }

    // Particle density given by constant properties
    const scalar rho = td.constProps().rho0();

    // Volume fraction to introduce during this timestep
    const scalar volFraction = volumeFraction(volume);

    // Duration of injection period during this timestep
    const scalar deltaT = min
    (
        continuousDt,
        min(time - SOI_, timeEnd() - time0_)
    );

    // Pad injection time if injection starts during this timestep
    const scalar padTime = max(0.0, SOI_ - time0_);

    // Introduce new parcels linearly with time
    for (label iParcel=0; iParcel<nParcels; iParcel++)
    {
        // Calculate the pseudo time of injection for parcel 'iParcel'
        scalar timeInj = time0_ + padTime + deltaT*iParcel/nParcels;

        // Determine injected parcel properties
        vector pos = position(iParcel, timeInj, owner_.meshInfo());

        // Diameter of parcels
        scalar d = d0(iParcel, timeInj);

        // Number of particles per parcel
        scalar nP = setNumberOfParticles
        (
            nParcels,
            d,
            volFraction,
            rho,
            volume
        );

        // Velocity of parcels
        vector U = velocity(iParcel, timeInj, owner_.meshInfo());

        // Determine the injection cell
        label cellI = -1;
        findInjectorCellAndPosition(cellI, pos);

        if (cellI >= 0)
        {
            scalar dt = time - timeInj;
            td.cloud().addNewParcel(pos, cellI, d, U, nP, dt);

            massInjected_ += nP*rho*mathematicalConstant::pi*pow3(d)/6.0;
            nParcelsAdded_++;
        }
    }

    postInjectCheck();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewInjectionModel.C"

// ************************************************************************* //
