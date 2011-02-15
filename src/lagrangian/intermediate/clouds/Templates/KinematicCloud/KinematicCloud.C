/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

#include "CollisionModel.H"
#include "DispersionModel.H"
#include "InjectionModel.H"
#include "PatchInteractionModel.H"
#include "PostProcessingModel.H"
#include "SurfaceFilmModel.H"

// * * * * * * * * * * * * * * cloudSolution * * * * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::cloudSolution::read()
{
    dict_.lookup("transient") >> transient_;
    dict_.lookup("coupled") >> coupled_;
    dict_.lookup("cellValueSourceCorrection") >> cellValueSourceCorrection_;

    if (steadyState())
    {
        dict_.lookup("calcFrequency") >> calcFrequency_;
        dict_.lookup("maxCo") >> maxCo_;
        dict_.lookup("maxTrackTime") >> maxTrackTime_;
        dict_.subDict("sourceTerms").lookup("resetOnStartup")
            >> resetSourcesOnStartup_;
    }
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
    transient_(false),
    calcFrequency_(1),
    maxCo_(0.3),
    iter_(1),
    deltaT_(0.0),
    coupled_(false),
    cellValueSourceCorrection_(false),
    maxTrackTime_(0.0),
    resetSourcesOnStartup_(true)
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
    transient_(cs.transient_),
    calcFrequency_(cs.calcFrequency_),
    maxCo_(cs.maxCo_),
    iter_(cs.iter_),
    deltaT_(cs.deltaT_),
    coupled_(cs.coupled_),
    cellValueSourceCorrection_(cs.cellValueSourceCorrection_),
    maxTrackTime_(cs.maxTrackTime_),
    resetSourcesOnStartup_(cs.resetSourcesOnStartup_)
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
    transient_(false),
    calcFrequency_(0),
    maxCo_(GREAT),
    iter_(0),
    deltaT_(0.0),
    coupled_(false),
    cellValueSourceCorrection_(false),
    maxTrackTime_(0.0),
    resetSourcesOnStartup_(false)
{}


template<class ParcelType>
Foam::KinematicCloud<ParcelType>::cloudSolution::~cloudSolution()
{}


template<class ParcelType>
Foam::scalar Foam::KinematicCloud<ParcelType>::cloudSolution::relaxCoeff
(
    const word& fieldName
) const
{
    return readScalar(sourceTermDict().subDict(fieldName).lookup("alpha"));
}


template<class ParcelType>
bool Foam::KinematicCloud<ParcelType>::cloudSolution::semiImplicit
(
    const word& fieldName
) const
{
    return readBool(sourceTermDict().subDict(fieldName).lookup("semiImplicit"));
}


template<class ParcelType>
bool Foam::KinematicCloud<ParcelType>::cloudSolution::solveThisStep() const
{
    return
        active_
     && (
            mesh_.time().outputTime()
         || (mesh_.time().timeIndex() % calcFrequency_ == 0)
        );
}


template<class ParcelType>
bool Foam::KinematicCloud<ParcelType>::cloudSolution::canEvolve()
{
    // Set the calculation time step
    if (transient_)
    {
        deltaT_ = mesh_.time().deltaTValue();
    }
    else
    {
        deltaT_ = maxTrackTime_;
    }

    return solveThisStep();
}


template<class ParcelType>
bool Foam::KinematicCloud<ParcelType>::cloudSolution::output() const
{
    return active_ && mesh_.time().outputTime();
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::setModels()
{
    collisionModel_.reset
    (
        CollisionModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    dispersionModel_.reset
    (
        DispersionModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    injectionModel_.reset
    (
        InjectionModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    patchInteractionModel_.reset
    (
        PatchInteractionModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    postProcessingModel_.reset
    (
        PostProcessingModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    surfaceFilmModel_.reset
    (
        SurfaceFilmModel<KinematicCloud<ParcelType> >::New
        (
            subModelProperties_,
            *this,
            g_
        ).ptr()
    );

    UIntegrator_.reset
    (
        vectorIntegrationScheme::New
        (
            "U",
            solution_.integrationSchemes()
        ).ptr()
    );
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::solve
(
    typename ParcelType::trackData& td
)
{
    if (solution_.transient())
    {
        td.cloud().preEvolve();

        evolveCloud(td);
    }
    else
    {
        td.cloud().storeState();

        td.cloud().preEvolve();

        evolveCloud(td);

        td.cloud().relaxSources(td.cloud().cloudCopy());
    }

    td.cloud().info();

    td.cloud().postEvolve();

    if (solution_.steadyState())
    {
        td.cloud().restoreState();
    }
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::preEvolve()
{
    Info<< "\nSolving cloud " << this->name() << endl;

    this->dispersion().cacheFields(true);
    forces_.cacheFields(true);
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
void Foam::KinematicCloud<ParcelType>::evolveCloud
(
    typename ParcelType::trackData& td
)
{
    if (solution_.coupled())
    {
        td.cloud().resetSourceTerms();
    }

    if (solution_.transient())
    {
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


        // Assume that motion will update the cellOccupancy as necessary
        // before it is required.
        motion(td);
    }
    else
    {
//        this->surfaceFilm().injectSteadyState(td);

        this->injection().injectSteadyState(td, solution_.deltaT());

        td.part() = ParcelType::trackData::tpLinearTrack;
        Cloud<ParcelType>::move(td,  solution_.deltaT());
    }
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
    Cloud<ParcelType>::move(td,  solution_.deltaT());

    td.part() = ParcelType::trackData::tpLinearTrack;
    Cloud<ParcelType>::move(td,  solution_.deltaT());

    // td.part() = ParcelType::trackData::tpRotationalTrack;
    // Cloud<ParcelType>::move(td);

    updateCellOccupancy();

    this->collision().collide();

    td.part() = ParcelType::trackData::tpVelocityHalfStep;
    Cloud<ParcelType>::move(td,  solution_.deltaT());
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::postEvolve()
{
    Info<< endl;

    if (debug)
    {
        this->writePositions();
    }

    this->dispersion().cacheFields(false);
    forces_.cacheFields(false);

    this->postProcessing().post();

    solution_.nextIter();
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::cloudReset(KinematicCloud<ParcelType>& c)
{
    Cloud<ParcelType>::cloudReset(c);

    rndGen_ = c.rndGen_;

    forces_.transfer(c.forces_);

    collisionModel_.reset(c.collisionModel_.ptr());
    dispersionModel_.reset(c.dispersionModel_.ptr());
    injectionModel_.reset(c.injectionModel_.ptr());
    patchInteractionModel_.reset(c.patchInteractionModel_.ptr());
    postProcessingModel_.reset(c.postProcessingModel_.ptr());

    UIntegrator_.reset(c.UIntegrator_.ptr());
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
    cloudCopyPtr_(NULL),
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
    constProps_(particleProperties_, solution_.active()),
    subModelProperties_
    (
        particleProperties_.subOrEmptyDict("subModels", solution_.active())
    ),
    rndGen_
    (
        label(0),
        particleProperties_.lookupOrDefault<label>("randomSampleSize", 100000)
    ),
    cellOccupancyPtr_(),
    rho_(rho),
    U_(U),
    mu_(mu),
    g_(g),
    forces_
    (
        *this,
        mesh_,
        subModelProperties_.subOrEmptyDict
        (
            "particleForces",
            solution_.active()
        ),
        solution_.active()
    ),
    collisionModel_(NULL),
    dispersionModel_(NULL),
    injectionModel_(NULL),
    patchInteractionModel_(NULL),
    postProcessingModel_(NULL),
    surfaceFilmModel_(NULL),
    UIntegrator_(NULL),
    UTrans_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + "UTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimMass*dimVelocity, vector::zero)
        )
    ),
    UCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "UCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",  dimMass/dimTime, 0.0)
        )
    )
{
    if (solution_.active())
    {
        setModels();
    }

    if (readFields)
    {
        ParcelType::readFields(*this);
    }

    if (solution_.resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class ParcelType>
Foam::KinematicCloud<ParcelType>::KinematicCloud
(
    KinematicCloud<ParcelType>& c,
    const word& name
)
:
    Cloud<ParcelType>(c.mesh_, name, c),
    kinematicCloud(),
    cloudCopyPtr_(NULL),
    mesh_(c.mesh_),
    particleProperties_(c.particleProperties_),
    solution_(c.solution_),
    constProps_(c.constProps_),
    subModelProperties_(c.subModelProperties_),
    rndGen_(c.rndGen_, true),
    cellOccupancyPtr_(NULL),
    rho_(c.rho_),
    U_(c.U_),
    mu_(c.mu_),
    g_(c.g_),
    forces_(c.forces_),
    collisionModel_(c.collisionModel_->clone()),
    dispersionModel_(c.dispersionModel_->clone()),
    injectionModel_(c.injectionModel_->clone()),
    patchInteractionModel_(c.patchInteractionModel_->clone()),
    postProcessingModel_(c.postProcessingModel_->clone()),
    surfaceFilmModel_(c.surfaceFilmModel_->clone()),
    UIntegrator_(c.UIntegrator_->clone()),
    UTrans_
    (
        new DimensionedField<vector, volMesh>
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
            c.UTrans_()
        )
    ),
    UCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + "UCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.UCoeff_()
        )
    )
{}


template<class ParcelType>
Foam::KinematicCloud<ParcelType>::KinematicCloud
(
    const fvMesh& mesh,
    const word& name,
    const KinematicCloud<ParcelType>& c
)
:
    Cloud<ParcelType>(mesh, name, IDLList<ParcelType>()),
    kinematicCloud(),
    cloudCopyPtr_(NULL),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            name + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    solution_(mesh),
    constProps_(),
    subModelProperties_(dictionary::null),
    rndGen_(0, 0),
    cellOccupancyPtr_(NULL),
    rho_(c.rho_),
    U_(c.U_),
    mu_(c.mu_),
    g_(c.g_),
    forces_(*this, mesh),
    collisionModel_(NULL),
    dispersionModel_(NULL),
    injectionModel_(NULL),
    patchInteractionModel_(NULL),
    postProcessingModel_(NULL),
    surfaceFilmModel_(NULL),
    UIntegrator_(NULL),
    UTrans_(NULL),
    UCoeff_(NULL)
{}


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

    const scalar carrierDt = this->db().time().deltaTValue();
    parcel.stepFraction() = (carrierDt - lagrangianDt)/carrierDt;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<KinematicCloud<ParcelType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::resetSourceTerms()
{
    UTrans().field() = vector::zero;
    UCoeff().field() = 0.0;
}


template<class ParcelType>
template<class Type>
void Foam::KinematicCloud<ParcelType>::relax
(
    DimensionedField<Type, volMesh>& field,
    const DimensionedField<Type, volMesh>& field0,
    const word& name
) const
{
    const scalar coeff = solution_.relaxCoeff(name);

    field = field0 + coeff*(field - field0);
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::relaxSources
(
    const KinematicCloud<ParcelType>& cloudOldTime
)
{
    this->relax(UTrans_(), cloudOldTime.UTrans(), "U");
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::evolve()
{
    if (solution_.canEvolve())
    {
        typename ParcelType::trackData td(*this);

        solve(td);
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
