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
#include "meshTools.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectionModel<CloudType>::readProps()
{
    IOobject propsDictHeader
    (
        "injectionProperties",
        owner_.db().time().timeName(),
        "uniform/Lagrangian"/owner_.name(),
        owner_.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (propsDictHeader.headerOk())
    {
        const IOdictionary propsDict(propsDictHeader);

        propsDict.readIfPresent("massInjected", massInjected_);
        propsDict.readIfPresent("nInjections", nInjections_);
        propsDict.readIfPresent("parcelsAddedTotal", parcelsAddedTotal_);
        propsDict.readIfPresent("timeStep0", timeStep0_);
    }
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::writeProps()
{
    if (owner_.db().time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                "injectionProperties",
                owner_.db().time().timeName(),
                "uniform/Lagrangian"/owner_.name(),
                owner_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        propsDict.add("massInjected", massInjected_);
        propsDict.add("nInjections", nInjections_);
        propsDict.add("parcelsAddedTotal", parcelsAddedTotal_);
        propsDict.add("timeStep0", timeStep0_);

        propsDict.regIOobject::write();
    }
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::prepareForNextTimeStep
(
    const scalar time0,
    const scalar time1,
    label& newParcels,
    scalar& newVolume
)
{
    // Initialise values
    newParcels = 0;
    newVolume = 0.0;

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
    newParcels = parcelsToInject(t0, t1);

    // Volume of parcels to inject
    newVolume = volumeToInject(t0, t1);

    // Hold previous time if no parcels, but non-zero volume fraction
    if ((newParcels == 0) && (newVolume > 0.0))
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
void Foam::InjectionModel<CloudType>::findCellAtPosition
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
        position += SMALL*(C - position);

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
            position += SMALL*(C - position);

            foundCell = owner_.mesh().pointInCell(position, cellI);
        }
        reduce(foundCell, orOp<bool>());
    }

    if (!foundCell)
    {
        FatalErrorIn
        (
            "Foam::InjectionModel<CloudType>::findCellAtPosition\n"
            "(\n"
            "    label&,\n"
            "    vector&\n"
            ")"
        )<< "Cannot find parcel injection cell. "
         << "Parcel position = " << p0 << nl
         << abort(FatalError);
    }
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::setNumberOfParticles
(
    const label parcels,
    const scalar volume,
    const scalar volumeFraction,
    const scalar diameter,
    const scalar rho
)
{
    scalar nP = 0.0;
    switch (parcelBasis_)
    {
        case pbMass:
        {
            nP = volumeFraction*massTotal_/parcels
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
                "Foam::scalar "
                "Foam::InjectionModel<CloudType>::setNumberOfParticles\n"
                "(\n"
                "    const label,\n"
                "    const scalar,\n"
                "    const scalar,\n"
                "    const scalar,\n"
                "    const scalar\n"
                ")"
            )<< "Unknown parcelBasis type" << nl
             << exit(FatalError);
        }
    }

    return nP;
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::postInjectCheck()
{
    if (parcelsAdded_ > 0)
    {
        Pout<< "\n--> Cloud: " << owner_.name() << nl
            << "    Added " << parcelsAdded_
            << " new parcels" << nl << endl;
    }

    // Increment total number of parcels added
    parcelsAddedTotal_ += parcelsAdded_;

    // Update time for start of next injection
    time0_ = owner_.db().time().value();

    // Increment number of injections
    nInjections_++;

    // Reset added parcels counter
    parcelsAdded_ = 0;

    // Write current state to properties file
    writeProps();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    SOI_(0.0),
    volumeTotal_(0.0),
    massTotal_(0.0),
    massInjected_(0.0),
    nInjections_(0),
    parcelsAdded_(0),
    parcelsAddedTotal_(0),
    parcelBasis_(pbNumber),
    time0_(0.0),
    timeStep0_(0.0)
{
    readProps();
}


template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    SOI_(readScalar(coeffDict_.lookup("SOI"))),
    volumeTotal_(0.0),
    massTotal_(dimensionedScalar(coeffDict_.lookup("massTotal")).value()),
    massInjected_(0.0),
    nInjections_(0),
    parcelsAdded_(0),
    parcelsAddedTotal_(0),
    parcelBasis_(pbNumber),
    time0_(owner.db().time().value()),
    timeStep0_(0.0)
{
    word parcelBasisType = coeffDict_.lookup("parcelBasisType");
    if (parcelBasisType == "mass")
    {
        parcelBasis_ = pbMass;
    }
    else if (parcelBasisType == "number")
    {
        parcelBasis_ = pbNumber;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::InjectionModel<CloudType>::InjectionModel\n"
            "(\n"
            "    const dictionary&,\n"
            "    CloudType&,\n"
            "    const word&\n"
            ")"
        )<< "parcelBasisType must be either 'number' or 'mass'" << nl
         << exit(FatalError);
    }

    readProps();
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
    parcelsAdded_ = 0;
    label newParcels = 0;
    scalar newVolume = 0.0;
    prepareForNextTimeStep(time0_, time, newParcels, newVolume);

    // Return if no parcels are required
    if (newParcels == 0)
    {
        postInjectCheck();
        return;
    }

    // Particle density given by constant properties
    const scalar rho = td.constProps().rho0();

    // Volume fraction to introduce during this timestep
    const scalar volFraction = volumeFraction(newVolume);

    // Duration of injection period during this timestep
    const scalar deltaT = min
    (
        continuousDt,
        min(time - SOI_, timeEnd() - time0_)
    );

    // Pad injection time if injection starts during this timestep
    const scalar padTime = max(0.0, SOI_ - time0_);

    // Introduce new parcels linearly with time
    for (label parcelI=0; parcelI<newParcels; parcelI++)
    {
        // Calculate the pseudo time of injection for parcel 'parcelI'
        scalar timeInj = time0_ + padTime + deltaT*parcelI/newParcels;

        // Determine the injection position and owner cell
        label cellI = -1;
        vector pos = vector::zero;
        setPositionAndCell(parcelI, timeInj, pos, cellI);

        if (cellI >= 0)
        {
            if (validInjection(parcelI))
            {
                // Diameter of parcels
                scalar d = d0(parcelI, timeInj);

                // Number of particles per parcel
                scalar nP = setNumberOfParticles
                (
                    newParcels,
                    newVolume,
                    volFraction,
                    d,
                    rho
                );

                // Velocity of parcels
                vector U = velocity(parcelI, timeInj);

                // Lagrangian timestep
                scalar dt = time - timeInj;

                // Apply corrections for 2-D cases
                meshTools::constrainToMeshCentre(owner_.mesh(), pos);
                meshTools::constrainDirection
                (
                    owner_.mesh(),
                    owner_.mesh().solutionD(),
                    U
                );

                // Add the new parcel
                td.cloud().addNewParcel(pos, cellI, d, U, nP, dt);
                massInjected_ += nP*rho*mathematicalConstant::pi*pow3(d)/6.0;
                parcelsAdded_++;
            }
        }
        else
        {
            WarningIn("Foam::InjectionModel<CloudType>::inject(TrackData& td)")
                << "Failed to inject new parcel:" << nl
                << "    id = " << parcelI << ", position = " << pos
                << nl << endl;
        }
    }

    postInjectCheck();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewInjectionModel.C"

// ************************************************************************* //
