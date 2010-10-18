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

#include "KinematicCloud.H"
#include "IntegrationScheme.H"
#include "interpolation.H"
#include "subCycleTime.H"

#include "DispersionModel.H"
#include "DragModel.H"
#include "InjectionModel.H"
#include "CollisionModel.H"
#include "PatchInteractionModel.H"
#include "PostProcessingModel.H"
#include "SurfaceFilmModel.H"

// * * * * * * * * * * * * * * cloudSolution * * * * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::cloudSolution::read()
{
    dict_.lookup("coupled") >> coupled_;
}


template<class ParcelType>
Foam::KinematicCloud<ParcelType>::cloudSolution::cloudSolution
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    active_(dict.lookup("active")),
    coupled_(false)
{
    if (active_)
    {
        read();
    }
}


template<class ParcelType>
Foam::KinematicCloud<ParcelType>::cloudSolution::cloudSolution
(
    const cloudSolution& cs
)
:
    mesh_(cs.mesh_),
    dict_(cs.dict_),
    active_(cs.active_),
    coupled_(cs.coupled_)
{}


template<class ParcelType>
Foam::KinematicCloud<ParcelType>::cloudSolution::cloudSolution
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dictionary::null),
    active_(false),
    coupled_(false)
{}


template<class ParcelType>
Foam::KinematicCloud<ParcelType>::cloudSolution::~cloudSolution()
{}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::preEvolve()
{
    this->dispersion().cacheFields(true);
    forces_.cacheFields(true, solution_.interpolationSchemes());
    updateCellOccupancy();
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::buildCellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        cellOccupancyPtr_.reset
        (
            new List<DynamicList<ParcelType*> >(mesh_.nCells())
        );
    }
    else if (cellOccupancyPtr_().size() != mesh_.nCells())
    {
        // If the size of the mesh has changed, reset the
        // cellOccupancy size

        cellOccupancyPtr_().setSize(mesh_.nCells());
    }

    List<DynamicList<ParcelType*> >& cellOccupancy = cellOccupancyPtr_();

    forAll(cellOccupancy, cO)
    {
        cellOccupancy[cO].clear();
    }

    forAllIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        cellOccupancy[iter().cell()].append(&iter());
    }
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::updateCellOccupancy()
{
    // Only build the cellOccupancy if the pointer is set, i.e. it has
    // been requested before.

    if (cellOccupancyPtr_.valid())
    {
        buildCellOccupancy();
    }
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::evolveCloud()
{
    autoPtr<interpolation<scalar> > rhoInterpolator =
        interpolation<scalar>::New
        (
            solution_.interpolationSchemes(),
            rho_
        );

    autoPtr<interpolation<vector> > UInterpolator =
        interpolation<vector>::New
        (
            solution_.interpolationSchemes(),
            U_
        );

    autoPtr<interpolation<scalar> > muInterpolator =
        interpolation<scalar>::New
        (
            solution_.interpolationSchemes(),
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

    label preInjectionSize = this->size();

    this->surfaceFilm().inject(td);

    // Update the cellOccupancy if the size of the cloud has changed
    // during the injection.
    if (preInjectionSize != this->size())
    {
        updateCellOccupancy();

        preInjectionSize = this->size();
    }

    this->injection().inject(td);

    if (solution_.coupled())
    {
        resetSourceTerms();
    }

    // Assume that motion will update the cellOccupancy as necessary
    // before it is required.
    motion(td);
}


template<class ParcelType>
void  Foam::KinematicCloud<ParcelType>::motion
(
    typename ParcelType::trackData& td
)
{
    // Sympletic leapfrog integration of particle forces:
    // + apply half deltaV with stored force
    // + move positions with new velocity
    // + calculate forces in new position
    // + apply half deltaV with new force

    label nSubCycles = collision().nSubCycles();

    if (nSubCycles > 1)
    {
        Info<< "    " << nSubCycles << " move-collide subCycles" << endl;

        subCycleTime moveCollideSubCycle
        (
            const_cast<Time&>(this->db().time()),
            nSubCycles
        );

        while(!(++moveCollideSubCycle).end())
        {
            moveCollide(td);
        }

        moveCollideSubCycle.endSubCycle();
    }
    else
    {
        moveCollide(td);
    }
}


template<class ParcelType>
void  Foam::KinematicCloud<ParcelType>::moveCollide
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

    updateCellOccupancy();

    this->collision().collide();

    td.part() = ParcelType::trackData::tpVelocityHalfStep;
    Cloud<ParcelType>::move(td);
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::postEvolve()
{
    if (debug)
    {
        this->writePositions();
    }

    this->dispersion().cacheFields(false);
    forces_.cacheFields(false, solution_.interpolationSchemes());

    this->postProcessing().post();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicCloud<ParcelType>::KinematicCloud
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
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    solution_(mesh_, particleProperties_.subDict("solution")),
    constProps_(particleProperties_),
    subModelProperties_(particleProperties_.subDict("subModels")),
    parcelTypeId_(readLabel(particleProperties_.lookup("parcelTypeId"))),
    cellValueSourceCorrection_
    (
        particleProperties_.lookup("cellValueSourceCorrection")
    ),
    rndGen_(label(0)),
    cellOccupancyPtr_(),
    rho_(rho),
    U_(U),
    mu_(mu),
    g_(g),
    forces_(mesh_, particleProperties_, g_.value()),
    dispersionModel_
    (
        DispersionModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        )
    ),
    dragModel_
    (
        DragModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        )
    ),
    injectionModel_
    (
        InjectionModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        )
    ),
    collisionModel_
    (
        CollisionModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        )
    ),
    patchInteractionModel_
    (
        PatchInteractionModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        )
    ),
    postProcessingModel_
    (
        PostProcessingModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        )
    ),
    surfaceFilmModel_
    (
        SurfaceFilmModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this,
            g
        )
    ),
    UIntegrator_
    (
        vectorIntegrationScheme::New
        (
            "U",
            solution_.integrationSchemes()
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
Foam::KinematicCloud<ParcelType>::~KinematicCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::checkParcelProperties
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

    scalar carrierDt = this->db().time().deltaTValue();
    parcel.stepFraction() = (carrierDt - lagrangianDt)/carrierDt;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::resetSourceTerms()
{
    UTrans_.field() = vector::zero;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::evolve()
{
    if (solution_.active())
    {
        preEvolve();

        evolveCloud();

        postEvolve();

        info();
        Info<< endl;
    }
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::info() const
{
    vector linearMomentum = linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());

    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar rotationalKineticEnergy = rotationalKineticEnergyOfSystem();
    reduce(rotationalKineticEnergy, sumOp<scalar>());

    Info<< "Cloud: " << this->name() << nl
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
    this->injection().info(Info);
    this->surfaceFilm().info(Info);
}


// ************************************************************************* //
