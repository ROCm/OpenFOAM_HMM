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

#include "ThermoCloud.H"
#include "interpolationCellPoint.H"
#include "ThermoParcel.H"

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::cloudReset(ThermoCloud<ParcelType>& c)
{
    KinematicCloud<ParcelType>::cloudReset(c);

    heatTransferModel_ = c.heatTransferModel_->clone();
    TIntegrator_ = c.TIntegrator_->clone();

    radiation_ = c.radiation_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoCloud<ParcelType>::ThermoCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    KinematicCloud<ParcelType>
    (
        cloudName,
        rho,
        U,
        thermo.thermo().mu(),
        g,
        false
    ),
    thermoCloud(),
    cloudCopyPtr_(NULL),
    constProps_(this->particleProperties()),
    thermo_(thermo),
    T_(thermo.thermo().T()),
    p_(thermo.thermo().p()),
    heatTransferModel_
    (
        HeatTransferModel<ThermoCloud<ParcelType> >::New
        (
            this->subModelProperties(),
            *this
        )
    ),
    TIntegrator_
    (
        scalarIntegrationScheme::New
        (
            "T",
            this->solution().integrationSchemes()
        )
    ),
    radiation_(this->subModelProperties().lookup("radiation")),
    hsTrans_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy, 0.0)
        )
    ),
    hsCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimTemperature, 0.0)
        )
    )

{
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
Foam::ThermoCloud<ParcelType>::ThermoCloud
(
    ThermoCloud<ParcelType>& c,
    const word& name
)
:
    KinematicCloud<ParcelType>(c, name),
    thermoCloud(),
    cloudCopyPtr_(NULL),
    constProps_(c.particleProperties_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    heatTransferModel_(c.heatTransferModel_->clone()),
    TIntegrator_(c.TIntegrator_->clone()),
    radiation_(c.radiation_),
    hsTrans_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsTrans()
        )
    ),
    hsCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsCoeff()
        )
    )
{}


template<class ParcelType>
Foam::ThermoCloud<ParcelType>::ThermoCloud
(
    const fvMesh& mesh,
    const word& name,
    const ThermoCloud<ParcelType>& c
)
:
    KinematicCloud<ParcelType>(mesh, name, c),
    thermoCloud(),
    cloudCopyPtr_(NULL),
    constProps_(c.particleProperties_),
    thermo_(c.thermo()),
    T_(c.T()),
    p_(c.p()),
    heatTransferModel_(NULL),
    TIntegrator_(NULL),
    radiation_(false),
    hsTrans_(NULL),
    hsCoeff_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoCloud<ParcelType>::~ThermoCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    KinematicCloud<ParcelType>::checkParcelProperties
    (
        parcel,
        lagrangianDt,
        fullyDescribed
    );

    if (!fullyDescribed)
    {
        parcel.T() = constProps_.T0();
        parcel.Cp() = constProps_.Cp0();
    }
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ThermoCloud<ParcelType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::resetSourceTerms()
{
    KinematicCloud<ParcelType>::resetSourceTerms();
    hsTrans_->field() = 0.0;
    hsCoeff_->field() = 0.0;
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::relaxSources
(
    const ThermoCloud<ParcelType>& cloudOldTime
)
{
    KinematicCloud<ParcelType>::relaxSources(cloudOldTime);

    this->relax(hsTrans_(), cloudOldTime.hsTrans(), "hs");
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::evolve()
{
    if (this->solution().active())
    {
        typename ParcelType::trackData td(*this);

        this->solve(td);
    }
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::info() const
{
    KinematicCloud<ParcelType>::info();
}


// ************************************************************************* //
