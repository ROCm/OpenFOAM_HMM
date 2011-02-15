/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

#include "ReactingCloud.H"

#include "CompositionModel.H"
#include "PhaseChangeModel.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::setModels()
{
    compositionModel_.reset
    (
        CompositionModel<ReactingCloud<ParcelType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    phaseChangeModel_.reset
    (
        PhaseChangeModel<ReactingCloud<ParcelType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::checkSuppliedComposition
(
    const scalarField& YSupplied,
    const scalarField& Y,
    const word& YName
)
{
    if (YSupplied.size() != Y.size())
    {
        FatalErrorIn
        (
            "ReactingCloud<ParcelType>::checkSuppliedComposition"
            "("
                "const scalarField&, "
                "const scalarField&, "
                "const word&"
            ")"
        )   << YName << " supplied, but size is not compatible with "
            << "parcel composition: " << nl << "    "
            << YName << "(" << YSupplied.size() << ") vs required composition "
            << YName << "(" << Y.size() << ")" << nl
            << abort(FatalError);
    }
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::cloudReset(ReactingCloud<ParcelType>& c)
{
    ThermoCloud<ParcelType>::cloudReset(c);

    compositionModel_.reset(c.compositionModel_.ptr());
    phaseChangeModel_.reset(c.phaseChangeModel_.ptr());

    dMassPhaseChange_ = c.dMassPhaseChange_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingCloud<ParcelType>::ReactingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    ThermoCloud<ParcelType>(cloudName, rho, U, g, thermo, false),
    reactingCloud(),
    cloudCopyPtr_(NULL),
    constProps_(this->particleProperties(), this->solution().active()),
    compositionModel_(NULL),
    phaseChangeModel_(NULL),
    rhoTrans_(thermo.carrier().species().size()),
    dMassPhaseChange_(0.0)
{
    if (this->solution().active())
    {
        setModels();
    }

    // Set storage for mass source fields and initialise to zero
    forAll(rhoTrans_, i)
    {
        const word& specieName = thermo.carrier().species()[i];
        rhoTrans_.set
        (
            i,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    this->name() + "rhoTrans_" + specieName,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
    }

    if (readFields)
    {
        ParcelType::readFields(*this);
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class ParcelType>
Foam::ReactingCloud<ParcelType>::ReactingCloud
(
    ReactingCloud<ParcelType>& c,
    const word& name
)
:
    ThermoCloud<ParcelType>(c, name),
    reactingCloud(),
    cloudCopyPtr_(NULL),
    constProps_(c.constProps_),
    compositionModel_(c.compositionModel_->clone()),
    phaseChangeModel_(c.phaseChangeModel_->clone()),
    rhoTrans_(c.rhoTrans_.size()),
    dMassPhaseChange_(c.dMassPhaseChange_)
{
    forAll(c.rhoTrans_, i)
    {
        const word& specieName = this->thermo().carrier().species()[i];
        rhoTrans_.set
        (
            i,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    this->name() + "rhoTrans_" + specieName,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.rhoTrans_[i]
            )
        );
    }
}


template<class ParcelType>
Foam::ReactingCloud<ParcelType>::ReactingCloud
(
    const fvMesh& mesh,
    const word& name,
    const ReactingCloud<ParcelType>& c
)
:
    ThermoCloud<ParcelType>(mesh, name, c),
    reactingCloud(),
    cloudCopyPtr_(NULL),
    constProps_(),
    compositionModel_(c.compositionModel_->clone()),
//    compositionModel_(NULL),
    phaseChangeModel_(NULL),
    rhoTrans_(0),
    dMassPhaseChange_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingCloud<ParcelType>::~ReactingCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    ThermoCloud<ParcelType>::checkParcelProperties
    (
        parcel,
        lagrangianDt,
        fullyDescribed
    );

    if (!fullyDescribed)
    {
        parcel.Y() = composition().YMixture0();
    }
    else
    {
        checkSuppliedComposition
        (
            parcel.Y(),
            composition().YMixture0(),
            "YMixture"
        );
    }

    // derived information - store initial mass
    parcel.mass0() = parcel.mass();
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ReactingCloud<ParcelType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::resetSourceTerms()
{
    ThermoCloud<ParcelType>::resetSourceTerms();
    forAll(rhoTrans_, i)
    {
        rhoTrans_[i].field() = 0.0;
    }
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::relaxSources
(
    const ReactingCloud<ParcelType>& cloudOldTime
)
{
    typedef DimensionedField<scalar, volMesh> dsfType;

    ThermoCloud<ParcelType>::relaxSources(cloudOldTime);

    forAll(rhoTrans_, fieldI)
    {
        dsfType& rhoT = rhoTrans_[fieldI];
        const dsfType& rhoT0 = cloudOldTime.rhoTrans()[fieldI];
        this->relax(rhoT, rhoT0, "rho");
    }
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename ParcelType::trackData td(*this);

        this->solve(td);
    }
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::addToMassPhaseChange(const scalar dMass)
{
    dMassPhaseChange_ += dMass;
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::info() const
{
    ThermoCloud<ParcelType>::info();

    Info<< "    Mass transfer phase change      = "
        << returnReduce(dMassPhaseChange_, sumOp<scalar>()) << nl;
}


// ************************************************************************* //
