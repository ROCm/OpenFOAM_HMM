/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "thermoSingleLayer.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "directMappedWallPolyPatch.H"
#include "specie.H"

// Sub-models
#include "heatTransferModel.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(thermoSingleLayer, 0);
        addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayer, mesh);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::surfaceFilmModels::thermoSingleLayer::read()
{
    // no additional properties to read
    return kinematicSingleLayer::read();
}


void Foam::surfaceFilmModels::thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar("zero", hsSp_.dimensions(), 0.0);
}


void Foam::surfaceFilmModels::thermoSingleLayer::correctThermoFields()
{
    switch (thermoModel_)
    {
        case tmConstant:
        {
            rho_ == dimensionedScalar(coeffs_.lookup("rho0"));
            mu_ == dimensionedScalar(coeffs_.lookup("mu0"));
            sigma_ == dimensionedScalar(coeffs_.lookup("sigma0"));
            Cp_ == dimensionedScalar(coeffs_.lookup("Cp0"));
            kappa_ == dimensionedScalar(coeffs_.lookup("kappa0"));

            break;
        }
        case tmSingleComponent:
        {
            const liquid& liq = thermo_.liquids().properties()[liquidId_];
            forAll(rho_, cellI)
            {
                const scalar T = T_[cellI];
                const scalar p = pPrimary_[cellI];
                rho_[cellI] = liq.rho(p, T);
                mu_[cellI] = liq.mu(p, T);
                sigma_[cellI] = liq.sigma(p, T);
                Cp_[cellI] = liq.Cp(p, T);
                kappa_[cellI] = liq.K(p, T);
            }

            rho_.correctBoundaryConditions();
            mu_.correctBoundaryConditions();
            sigma_.correctBoundaryConditions();
            Cp_.correctBoundaryConditions();
            kappa_.correctBoundaryConditions();

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void Foam::surfaceFilmModels::thermoSingleLayer::"
                "correctThermoFields()"
            )   << "Unknown thermoType enumeration" << abort(FatalError);
        }
    }
}


void Foam::surfaceFilmModels::thermoSingleLayer::updateSurfaceTemperatures()
{
    // Push boundary film temperature values into internal field
    for (label i=0; i<filmBottomPatchIDs_.size(); i++)
    {
        label patchI = filmBottomPatchIDs_[i];
        const polyPatch& pp = filmRegion_.boundaryMesh()[patchI];
        UIndirectList<scalar>(Tw_, pp.faceCells()) =
            T_.boundaryField()[patchI];
    }
    Tw_.correctBoundaryConditions();

    // Update film surface temperature
    Ts_ = T_;
    Ts_.correctBoundaryConditions();
}


void Foam::surfaceFilmModels::thermoSingleLayer::transferPrimaryRegionFields()
{
    kinematicSingleLayer::transferPrimaryRegionFields();

    // Update primary region fileds via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    hsSp_.correctBoundaryConditions();

    // Convert accummulated source terms into per unit area per unit time
    // Note: boundary values will still have original (neat) values
    const scalar deltaT = time_.deltaTValue();
    hsSp_.field() /= magSf_*deltaT;
}


void Foam::surfaceFilmModels::thermoSingleLayer::updateSubmodels()
{
    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update phase change
    phaseChange_->correct
    (
        time_.deltaTValue(),
        massPhaseChangeForPrimary_,
        energyPhaseChangeForPrimary_
    );
    massPhaseChangeForPrimary_.correctBoundaryConditions();

    // Update kinematic sub-models
    kinematicSingleLayer::updateSubmodels();

    // Update source fields
    hsSp_ -= energyPhaseChangeForPrimary_/magSf_/time_.deltaT();
}


Foam::tmp<Foam::fvScalarMatrix> Foam::surfaceFilmModels::thermoSingleLayer::q
(
    volScalarField& hs
) const
{
    const dimensionedScalar Tstd("Tstd", dimTemperature, specie::Tstd);

    return
    (
      - fvm::Sp(htcs_->h()/Cp_, hs) - htcs_->h()*(Tstd - TPrimary_)
      - fvm::Sp(htcw_->h()/Cp_, hs) - htcw_->h()*(Tstd - Tw_)
    );
}


void Foam::surfaceFilmModels::thermoSingleLayer::solveEnergy()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::solveEnergy()" << endl;
    }

    updateSurfaceTemperatures();

    volScalarField mLossCoeff
    (
        "mLossCoeff",
        massForPrimary_/magSf_/time_.deltaT()
    );

//    dimensionedScalar hs0("SMALL", hs_.dimensions(), SMALL);

    solve
    (
        fvm::ddt(deltaRho_, hs_)
      + fvm::div(phi_, hs_)
     ==
//        fvm::Sp(hsSp_/(hs_ + hs0), hs_)
        hsSp_
      + q(hs_)
      - fvm::Sp(mLossCoeff, hs_)
    );

    correctThermoFields();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::thermoSingleLayer::thermoSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    kinematicSingleLayer(modelType, mesh, g),
    thermo_(mesh.lookupObject<SLGThermo>("SLGThermo")),
    liquidId_(thermo_.liquidId(coeffs_.lookup("liquid"))),
    Cp_
    (
        IOobject
        (
            "Cp",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_,
        dimensionedScalar("Cp", dimEnergy/dimMass/dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_,
        dimensionedScalar
        (
            "kappa",
            dimEnergy/dimTime/dimLength/dimTemperature,
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    T_
    (
        IOobject
        (
            "Tf",
            time_.timeName(),
            filmRegion_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_
    ),
    Ts_
    (
        IOobject
        (
            "Tsf",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    Tw_
    (
        IOobject
        (
            "Twf",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    hs_
    (
        IOobject
        (
            "hsf",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar("zero", dimEnergy/dimMass, 0.0),
        T_.boundaryField().types()
    ),

    hsSp_
    (
        IOobject
        (
            "hsSp",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        pSp_.boundaryField().types()
    ),

    hsSpPrimary_
    (
        IOobject
        (
            hsSp_.name(), // must have same name as hSp_ to enable mapping
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", hsSp_.dimensions(), 0.0)
    ),

    TPrimary_
    (
        IOobject
        (
            "T", // must have same name as T on primary region to enable mapping
            time_.timeName(),
            filmRegion_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_
    ),

    YPrimary_(),

    htcs_
    (
        heatTransferModel::New(*this, coeffs_.subDict("upperSurfaceModels"))
    ),
    htcw_
    (
        heatTransferModel::New(*this, coeffs_.subDict("lowerSurfaceModels"))
    ),
    phaseChange_(phaseChangeModel::New(*this, coeffs_)),
    energyPhaseChangeForPrimary_
    (
        IOobject
        (
            "energyPhaseChangeForPrimary",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar("zero", dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        time_.timeName(),
                        filmRegion_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    filmRegion_,
                    dimensionedScalar("zero", dimless, 0.0),
                    pSp_.boundaryField().types()
                )
            );
        }
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::thermoSingleLayer::~thermoSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::thermoSingleLayer::addSources
(
    const label patchI,
    const label faceI,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchI,
        faceI,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    if (debug)
    {
        Info<< "    energy   = " << energySource << nl << endl;
    }

    hsSpPrimary_.boundaryField()[patchI][faceI] += energySource;
}


void Foam::surfaceFilmModels::thermoSingleLayer::preEvolveFilm()
{
    if (!initialisedThermo_)
    {
        // Retreive pressure from primary region
        pPrimary_.correctBoundaryConditions();

        // Correct (temperature dependent) thermo fields
        correctThermoFields();

        // Update derived fields
        hs_ == hs(T_);
        deltaRho_ == delta_*rho_;
        phi_ = fvc::interpolate(deltaRho_*U_) & filmRegion_.Sf();
        initialisedThermo_ = true;
    }
}


void Foam::surfaceFilmModels::thermoSingleLayer::evolveFilm()
{
    transferPrimaryRegionFields();

    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=0; oCorr<nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> tpu = this->pu();

        // Implicit pressure source coefficient
        tmp<volScalarField> tpp = this->pp();

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

        // Solve energy for hs_
        solveEnergy();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn());
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update temperature using latest hs_
    T_ == T(hs_);

    // Update film wall and surface velocities
    updateSurfaceVelocities();

    // Update film wall and surface temperatures
    updateSurfaceTemperatures();

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::Cp() const
{
    return Cp_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::kappa() const
{
    return kappa_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::T() const
{
    return T_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::Ts() const
{
    return Ts_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::Tw() const
{
    return Tw_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::hs() const
{
    return hs_;
}


void Foam::surfaceFilmModels::thermoSingleLayer::info() const
{
    kinematicSingleLayer::info();

    Info<< indent << "min/max(T)         = " << min(T_).value() << ", "
        << max(T_).value() << nl;

    phaseChange_->info();
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::surfaceFilmModels::thermoSingleLayer::Srho() const
{
    tmp<DimensionedField<scalar, volMesh> > tSrho
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Srho",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    scalarField& Srho = tSrho();
    const scalarField& V = mesh_.V();
    const scalar dt = time_.deltaTValue();

    forAll(filmBottomPatchIDs_, i)
    {
        const label primaryPatchI = primaryPatchIDs_[i];
        const directMappedWallPolyPatch& wpp =
            refCast<const directMappedWallPolyPatch>
            (
                 mesh_.boundaryMesh()[primaryPatchI]
            );

        const mapDistribute& distMap = wpp.map();

        const label filmPatchI = filmBottomPatchIDs_[i];

        scalarField patchMass
        (
            massPhaseChangeForPrimary_.boundaryField()[filmPatchI]
        );

        distMap.distribute(patchMass);

        const labelUList& cells = wpp.faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] = patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::surfaceFilmModels::thermoSingleLayer::Srho(const label i) const
{
    const label vapId =
        thermo_.carrierId(thermo_.liquids().components()[liquidId_]);

    tmp<DimensionedField<scalar, volMesh> > tSrho
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Srho(i)",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho();
        const scalarField& V = mesh_.V();
        const scalar dt = time_.deltaTValue();

        forAll(filmBottomPatchIDs_, i)
        {
            const label primaryPatchI = primaryPatchIDs_[i];
            const directMappedWallPolyPatch& wpp =
                refCast<const directMappedWallPolyPatch>
                (
                    mesh_.boundaryMesh()[primaryPatchI]
                );

            const mapDistribute& distMap = wpp.map();

            const label filmPatchI = filmBottomPatchIDs_[i];

            scalarField patchMass
            (
                massPhaseChangeForPrimary_.boundaryField()[filmPatchI]
            );

            distMap.distribute(patchMass);

            const labelUList& cells = wpp.faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] = patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::surfaceFilmModels::thermoSingleLayer::Sh() const
{
    tmp<DimensionedField<scalar, volMesh> > tSh
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Sh",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    // All of enthalpy change due to phase change added to film energy equation

/*
    scalarField& Sh = tSh();
    const scalarField& V = mesh_.V();
    const scalar dt = time_.deltaTValue();

    forAll(filmBottomPatchIDs_, i)
    {
        const label primaryPatchI = primaryPatchIDs_[i];
        const directMappedWallPolyPatch& wpp =
            refCast<const directMappedWallPolyPatch>
            (
                mesh_.boundaryMesh()[primaryPatchI]
            );

        const mapDistribute& distMap = wpp.map();

        const label filmPatchI = filmBottomPatchIDs_[i];

        scalarField patchEnergy
        (
            energyPhaseChangeForPrimary_.boundaryField()[filmPatchI]
        );
        distMap.distribute(patchEnergy);

        const labelUList& cells = wpp.faceCells();

        forAll(patchMass, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }
*/
    return tSh;
}


// ************************************************************************* //
