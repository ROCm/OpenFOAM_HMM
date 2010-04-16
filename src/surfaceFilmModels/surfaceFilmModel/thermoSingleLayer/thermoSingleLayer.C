/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
#include "phaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"

// Sub-models
#include "injectionModel.H"

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
    if (kinematicSingleLayer::read())
    {
        coeffs_.lookup("htcw") >> htcw_;
        coeffs_.lookup("htcs") >> htcs_;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::surfaceFilmModels::thermoSingleLayer::initialise()
{
    if (debug)
    {
        Pout<< "thermoSingleLayer::initialise()" << endl;
    }

    kinematicSingleLayer::initialise();

    hs_ == hs(T_);
}


void Foam::surfaceFilmModels::thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar("zero", hsSp_.dimensions(), 0.0);
}


void Foam::surfaceFilmModels::thermoSingleLayer::transferPrimaryRegionFields()
{
    kinematicSingleLayer::transferPrimaryRegionFields();

    // Update temperature from primary region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    hsSp_.correctBoundaryConditions();

    // Convert accummulated source terms into per unit area per unit time
    // Note: boundary values will still have original (neat) values
    const scalar deltaT = filmRegion_.time().deltaTValue();
    hsSp_.field() /= magSf_*deltaT;
}


void Foam::surfaceFilmModels::thermoSingleLayer::updateSubmodels()
{
    kinematicSingleLayer::updateSubmodels();

    const dimensionedScalar deltaT = filmRegion_.time().deltaT();
    hsSp_ -= massForPrimary_*hs_/magSf_/deltaT;

    // Update the sub-models
    phaseChange_->correct();
}


Foam::tmp<Foam::fvScalarMatrix> Foam::surfaceFilmModels::thermoSingleLayer::q
(
    volScalarField& hs
) const
{
    DimensionedField<scalar, volMesh> Tw
    (
        IOobject
        (
            "Tw",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        filmRegion_,
        dimensionedScalar("zero", dimTemperature, 0.0)
    );

    for (label i=0; i<filmBottomPatchIDs_.size(); i++)
    {
        label patchI = filmBottomPatchIDs_[i];
        const polyPatch& pp = filmRegion_.boundaryMesh()[patchI];
        UIndirectList<scalar>(Tw, pp.faceCells()) =
            TPrimary_.boundaryField()[patchI];
    }

    // TODO: Use T at film thickness instead of cell value
    DimensionedField<scalar, volMesh> Ts =
        TPrimary_.dimensionedInternalField();

    return
    (
      - fvm::Sp(htcs_/cp_, hs)
      + htcs_*(dimensionedScalar("Tstd", dimTemperature, 298.15) - Ts)
      - fvm::Sp(htcw_/cp_, hs)
      + htcw_*(dimensionedScalar("Tstd", dimTemperature, 298.15) - Tw)
    );
}


void Foam::surfaceFilmModels::thermoSingleLayer::solveEnergy()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::solveEnergy()" << endl;
    }

    solve
    (
        fvm::ddt(deltaRho_, hs_)
      + fvm::div(phi_, hs_)
     ==
        hsSp_
      + q(hs_)
    );
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
    cp_
    (
        IOobject
        (
            "cp",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar(coeffs_.lookup("cp"))
    ),

    htcw_
    (
        dimensionedScalar("htc", dimEnergy/dimTime/dimArea/dimTemperature, 0.0)
    ),
    htcs_
    (
        dimensionedScalar("htc", dimEnergy/dimTime/dimArea/dimTemperature, 0.0)
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

    phaseChange_(phaseChangeModel::New(*this, coeffs_)),

    hsSpDetach_(filmRegion_.nCells(), 0.0)
{
    read();

    initialise();
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


void Foam::surfaceFilmModels::thermoSingleLayer::evolveFilm()
{
    transferPrimaryRegionFields();

    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=0; oCorr<nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> pu = this->pu();

        // Implicit pressure source coefficient
        tmp<volScalarField> pp = this->pp();

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(pu(), pp());

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve energy for hs_
            solveEnergy();

            // Solve thickness for delta_
            solveThickness(pu(), pp(), UEqn());
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update temperature using latest hs_
    T_ == T(hs_);

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::T() const
{
    return T_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::cp() const
{
    return cp_;
}


void Foam::surfaceFilmModels::thermoSingleLayer::info() const
{
    kinematicSingleLayer::info();

    Info<< indent<< "min/max(T)      = " << min(T_).value() << ", "
        << max(T_).value() << nl;
}


// ************************************************************************* //
