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

#include "KinematicCloud.H"
#include "IntegrationScheme.H"
#include "interpolation.H"

#include "DispersionModel.H"
#include "DragModel.H"
#include "InjectionModel.H"
#include "WallInteractionModel.H"


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::addNewParcel
(
    const vector& position,
    const label cellId,
    const scalar d,
    const vector& U,
    const scalar nParticles,
    const scalar lagrangianDt
)
{
    ParcelType* pPtr = new ParcelType
    (
        *this,
        position,
        cellId,
        parcelTypeId_,
        nParticles,
        d,
        U,
        constProps_
    );

    scalar continuousDt = this->db().time().deltaT().value();
    pPtr->stepFraction() = (continuousDt - lagrangianDt)/continuousDt;

    addParticle(pPtr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicCloud<ParcelType>::KinematicCloud
(
    const word& cloudType,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g
)
:
    Cloud<ParcelType>(rho.mesh(), cloudType, false),
    kinematicCloud(),
    cloudType_(cloudType),
    mesh_(rho.mesh()),
    particleProperties_
    (
        IOobject
        (
            cloudType + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    constProps_(particleProperties_),
    parcelTypeId_(readLabel(particleProperties_.lookup("parcelTypeId"))),
    coupled_(particleProperties_.lookup("coupled")),
    rndGen_(label(0)),
    rho_(rho),
    U_(U),
    mu_(mu),
    g_(g),
    forces_(mesh_, particleProperties_, g_.value()),
    interpolationSchemes_(particleProperties_.subDict("interpolationSchemes")),
    dispersionModel_
    (
        DispersionModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    dragModel_
    (
        DragModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    injectionModel_
    (
        InjectionModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    wallInteractionModel_
    (
        WallInteractionModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
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
        dimensionedVector("zero", dimensionSet(1, 1, -1, 0, 0), vector::zero)
    ),
    UCoeff_
    (
        IOobject
        (
            this->name() + "UCoeff",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -1, 0, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicCloud<ParcelType>::~KinematicCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::resetSourceTerms()
{
    UTrans_.field() = vector::zero;
    UCoeff_.field() = 0.0;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::evolve()
{
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

    if (debug)
    {
        this->writePositions();
    }

    if (coupled_)
    {
        resetSourceTerms();
    }

    Cloud<ParcelType>::move(td);
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::info() const
{
    Info<< "Cloud name: " << this->name() << nl
        << "    Parcels added during this run   = "
        << returnReduce(this->injection().parcelsAddedTotal(), sumOp<label>())
            << nl
        << "    Mass introduced during this run = "
        << returnReduce(this->injection().massInjected(), sumOp<scalar>())
            << nl
        << "    Current number of parcels       = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl;
}


// ************************************************************************* //
