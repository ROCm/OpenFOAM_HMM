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

#include "InteractingKinematicCloud.H"
#include "IntegrationScheme.H"
#include "interpolation.H"
#include "subCycleTime.H"

#include "DispersionModel.H"
#include "DragModel.H"
#include "InjectionModel.H"
#include "CollisionModel.H"
#include "PatchInteractionModel.H"
#include "PostProcessingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::InteractingKinematicCloud<ParcelType>::InteractingKinematicCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
    Cloud<ParcelType>(rho.mesh(), cloudName, false),
    kinematicCloud(),
    mesh_(rho.mesh()),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    constProps_(particleProperties_),
    parcelTypeId_(readLabel(particleProperties_.lookup("parcelTypeId"))),
    coupled_(particleProperties_.lookup("coupled")),
    cellValueSourceCorrection_
    (
        particleProperties_.lookup("cellValueSourceCorrection")
    ),
    rndGen_(label(0)),
    rho_(rho),
    U_(U),
    mu_(mu),
    g_(g),
    forces_(mesh_, particleProperties_, g_.value()),
    interpolationSchemes_(particleProperties_.subDict("interpolationSchemes")),
    dispersionModel_
    (
        DispersionModel<InteractingKinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    dragModel_
    (
        DragModel<InteractingKinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    injectionModel_
    (
        InjectionModel<InteractingKinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    collisionModel_
    (
        CollisionModel<InteractingKinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    patchInteractionModel_
    (
        PatchInteractionModel<InteractingKinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    postProcessingModel_
    (
        PostProcessingModel<InteractingKinematicCloud<ParcelType> >::New
        (
            this->particleProperties_,
            *this
        )
    ),
    UIntegrator_
    (
        vectorIntegrationScheme::New
        (
            "U",
            particleProperties_.subDict("integrationSchemes")
        )
    ),
    UTrans_
    (
        IOobject
        (
            this->name() + "UTrans",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimMass*dimVelocity, vector::zero)
    )
{
    if (readFields)
    {
        ParcelType::readFields(*this);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::InteractingKinematicCloud<ParcelType>::~InteractingKinematicCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::InteractingKinematicCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    if (!fullyDescribed)
    {
        parcel.rho() = constProps_.rho0();
    }

    scalar carrierDt = this->db().time().deltaT().value();
    parcel.stepFraction() = (carrierDt - lagrangianDt)/carrierDt;
}


template<class ParcelType>
void Foam::InteractingKinematicCloud<ParcelType>::resetSourceTerms()
{
    UTrans_.field() = vector::zero;
}


template<class ParcelType>
void Foam::InteractingKinematicCloud<ParcelType>::preEvolve()
{
    this->dispersion().cacheFields(true);
    forces_.cacheFields(true);
}


template<class ParcelType>
void Foam::InteractingKinematicCloud<ParcelType>::postEvolve()
{
    if (debug)
    {
        this->writePositions();
    }

    this->dispersion().cacheFields(false);
    forces_.cacheFields(false);

    this->postProcessing().post();
}


template<class ParcelType>
void Foam::InteractingKinematicCloud<ParcelType>::evolve()
{
    preEvolve();

    autoPtr<interpolation<scalar> > rhoInterpolator =
        interpolation<scalar>::New
        (
            interpolationSchemes_,
            rho_
        );

    autoPtr<interpolation<vector> > UInterpolator =
        interpolation<vector>::New
        (
            interpolationSchemes_,
            U_
        );

    autoPtr<interpolation<scalar> > muInterpolator =
        interpolation<scalar>::New
        (
            interpolationSchemes_,
            mu_
        );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterpolator(),
        UInterpolator(),
        muInterpolator(),
        g_.value()
    );

    this->injection().inject(td);

    if (coupled_)
    {
        resetSourceTerms();
    }

    // Sympletic leapfrog integration of particle forces:
    // + apply half deltaV with stored force
    // + move positions with new velocity
    // + calculate forces in new position
    // + apply half deltaV with new force

    label nSubCycles = collision().nSubCycles();

    if (nSubCycles > 1)
    {
        subCycleTime moveCollideSubCycle
        (
            const_cast<Time&>(this->db().time()),
            nSubCycles
        );

        while(!(++moveCollideSubCycle).end())
        {
            Info<< "subCycle time = " << this->db().time().timeName() << endl;

            moveCollide(td);
        }

        moveCollideSubCycle.endSubCycle();
    }
    else
    {
        moveCollide(td);
    }

    postEvolve();
}


template<class ParcelType>
void  Foam::InteractingKinematicCloud<ParcelType>::moveCollide
(
    typename ParcelType::trackData& td
)
{
    td.part() = ParcelType::trackData::tpVelocityHalfStep;
    Cloud<ParcelType>::move(td);

    td.part() = ParcelType::trackData::tpLinearTrack;
    Cloud<ParcelType>::move(td);

    // td.part() = ParcelType::trackData::tpRotationalTrack;
    // Cloud<ParcelType>::move(td);

    this->collision().collide();

    td.part() = ParcelType::trackData::tpVelocityHalfStep;
    Cloud<ParcelType>::move(td);
}


template<class ParcelType>
void Foam::InteractingKinematicCloud<ParcelType>::info() const
{
    vector linearMomentum = linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());

    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar rotationalKineticEnergy = rotationalKineticEnergyOfSystem();
    reduce(rotationalKineticEnergy, sumOp<scalar>());

    Info<< "Cloud: " << this->name() << nl
        << "    Total number of parcels added   = "
        << returnReduce(this->injection().parcelsAddedTotal(), sumOp<label>())
            << nl
        << "    Total mass introduced           = "
        << returnReduce(this->injection().massInjected(), sumOp<scalar>())
            << nl
        << "    Current number of parcels       = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl
        << "    Linear momentum                 = "
        << linearMomentum << nl
        << "   |Linear momentum|                = "
        << mag(linearMomentum) << nl
        << "    Linear kinetic energy           = "
        << linearKineticEnergy << nl
        << "    Rotational kinetic energy       = "
        << rotationalKineticEnergy << nl;
}


// ************************************************************************* //
