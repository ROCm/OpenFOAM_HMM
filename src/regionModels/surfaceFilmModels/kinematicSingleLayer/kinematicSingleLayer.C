/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "kinematicSingleLayer.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"
#include "directMappedWallPolyPatch.H"
#include "mapDistribute.H"

// Sub-models
#include "injectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kinematicSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmModel, kinematicSingleLayer, mesh);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool kinematicSingleLayer::read()
{
    if (surfaceFilmModel::read())
    {
        const dictionary& solution = this->solution().subDict("PISO");
        solution.lookup("momentumPredictor") >> momentumPredictor_;
        solution.lookup("nOuterCorr") >> nOuterCorr_;
        solution.lookup("nCorr") >> nCorr_;
        solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;

        coeffs_.lookup("Cf") >> Cf_;
        coeffs_.lookup("deltaStable") >> deltaStable_;

        return true;
    }
    else
    {
        return false;
    }
}


void kinematicSingleLayer::correctThermoFields()
{
    if (thermoModel_ == tmConstant)
    {
        rho_ == dimensionedScalar(coeffs_.lookup("rho0"));
        mu_ == dimensionedScalar(coeffs_.lookup("mu0"));
        sigma_ == dimensionedScalar(coeffs_.lookup("sigma0"));
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::surfaceFilmModels::kinematicSingleLayer::"
            "correctThermo()"
        )   << "Kinematic surface film must use "
            << thermoModelTypeNames_[thermoModel_] << "thermodynamics" << endl;
    }
}


void kinematicSingleLayer::resetPrimaryRegionSourceTerms()
{
    rhoSpPrimary_ == dimensionedScalar("zero", rhoSp_.dimensions(), 0.0);
    USpPrimary_ == dimensionedVector("zero", USp_.dimensions(), vector::zero);
    pSpPrimary_ == dimensionedScalar("zero", pSp_.dimensions(), 0.0);
}


void kinematicSingleLayer::transferPrimaryRegionThermoFields()
{
    // Update fields from primary region via direct mapped
    // (coupled) boundary conditions
    UPrimary_.correctBoundaryConditions();
    pPrimary_.correctBoundaryConditions();
    rhoPrimary_.correctBoundaryConditions();
    muPrimary_.correctBoundaryConditions();
}


void kinematicSingleLayer::transferPrimaryRegionSourceFields()
{
    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    rhoSp_.correctBoundaryConditions();
    USp_.correctBoundaryConditions();
    pSp_.correctBoundaryConditions();

    // Convert accummulated source terms into per unit area per unit time
    // Note: boundary values will still have original (neat) values
    const scalar deltaT = time_.deltaTValue();
    rhoSp_.field() /= magSf()*deltaT;
    USp_.field() /= magSf()*deltaT;
    pSp_.field() /= magSf()*deltaT;

    // reset transfer to primary fields
    massForPrimary_ == dimensionedScalar("zero", dimMass, 0.0);
    diametersForPrimary_ == dimensionedScalar("zero", dimLength, -1.0);
}


tmp<volScalarField> kinematicSingleLayer::pu()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pu",
                time_.timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pPrimary_                  // pressure (mapped from primary region)
          + pSp_                           // accumulated particle impingement
          - fvc::laplacian(sigma_, delta_) // surface tension
        )
    );
}


tmp<volScalarField> kinematicSingleLayer::pp()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pp",
                time_.timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -rho_*gNormClipped() // hydrostatic effect only
        )
    );
}


void kinematicSingleLayer::correctDetachedFilm()
{
    tmp<volScalarField> tgNorm(this->gNorm());
    const scalarField& gNorm = tgNorm();
    const scalarField& magSf = this->magSf();

    forAll(gNorm, i)
    {
        if (gNorm[i] > SMALL)
        {
            const scalar ddelta = max(0.0, delta_[i] - deltaStable_.value());
            massForPrimary_[i] =
                max
                (
                    0.0,
                    ddelta*rho_[i]*magSf[i] - massPhaseChangeForPrimary_[i]
                );
        }
    }
}


void kinematicSingleLayer::updateSubmodels()
{
    correctDetachedFilm();

    // Update injection model - mass returned is actual mass injected
    injection_->correct(massForPrimary_, diametersForPrimary_);

    // Update cumulative detached mass counter
    injectedMassTotal_ += sum(massForPrimary_.field());

    // Push values to boundaries ready for transfer to the primary region
    massForPrimary_.correctBoundaryConditions();
    diametersForPrimary_.correctBoundaryConditions();

    // Update source fields
    const dimensionedScalar deltaT = time().deltaT();
    rhoSp_ -= (massForPrimary_ + massPhaseChangeForPrimary_)/magSf()/deltaT;
}


void kinematicSingleLayer::continuityCheck()
{
    const volScalarField deltaRho0 = deltaRho_;

    solveContinuity();

    if (debug)
    {
        const volScalarField mass(deltaRho_*magSf());
        const dimensionedScalar totalMass =
            fvc::domainIntegrate(mass)
          + dimensionedScalar("SMALL", dimMass*dimVolume, ROOTVSMALL);

        const scalar sumLocalContErr =
            (
                fvc::domainIntegrate(mag(mass - magSf()*deltaRho0))/totalMass
            ).value();

        const scalar globalContErr =
            (
                fvc::domainIntegrate(mass - magSf()*deltaRho0)/totalMass
            ).value();

        cumulativeContErr_ += globalContErr;

        Info<< "Surface film: " << type() << nl
            << "    time step continuity errors: sum local = "
            << sumLocalContErr << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr_ << endl;
    }
}


void kinematicSingleLayer::solveContinuity()
{
    if (debug)
    {
        Info<< "kinematicSingleLayer::solveContinuity()" << endl;
    }

    solve
    (
        fvm::ddt(deltaRho_)
      + fvc::div(phi_)
     ==
        rhoSp_
    );
}


void kinematicSingleLayer::updateSurfaceVelocities()
{
    // Push boundary film velocity values into internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchI = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchI];
        UIndirectList<vector>(Uw_, pp.faceCells()) =
            U_.boundaryField()[patchI];
    }
    Uw_ -= nHat()*(Uw_ & nHat());
    Uw_.correctBoundaryConditions();

    // TODO: apply quadratic profile to determine surface velocity
    Us_ = U_;
    Us_.correctBoundaryConditions();
}


tmp<fvVectorMatrix> kinematicSingleLayer::tau(volVectorField& U) const
{
    // Calculate shear stress
    volScalarField Cs("Cs", rho_*Cf_*mag(Us_ - U));
    volScalarField Cw
    (
        "Cw",
        mu_/(0.3333*(delta_ + dimensionedScalar("SMALL", dimLength, SMALL)))
    );
    Cw.min(1.0e+06);

    return
    (
       - fvm::Sp(Cs, U) + Cs*Us_
       - fvm::Sp(Cw, U) + Cw*Uw_
    );
}


tmp<Foam::fvVectorMatrix> kinematicSingleLayer::solveMomentum
(
    const volScalarField& pu,
    const volScalarField& pp
)
{
    if (debug)
    {
        Info<< "kinematicSingleLayer::solveMomentum()" << endl;
    }

    updateSurfaceVelocities();

    // Momentum
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(deltaRho_, U_)
      + fvm::div(phi_, U_)
     ==
        USp_
      + tau(U_)
      + fvc::grad(sigma_)
      - fvm::Sp
        (
            (massForPrimary_ + massPhaseChangeForPrimary_)
           /magSf()/time().deltaT(),
            U_
        )
    );

    fvVectorMatrix& UEqn = tUEqn();

    UEqn.relax();

    if (momentumPredictor_)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
              - fvc::interpolate(delta_)
              * (
                    regionMesh().magSf()
                  * (
                        fvc::snGrad(pu, "snGrad(p)")
                      + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
                      + fvc::snGrad(delta_)*fvc::interpolate(pp)
                    )
                  - (fvc::interpolate(rho_*gTan()) & regionMesh().Sf())
                )
            )
        );

        // Remove any patch-normal components of velocity
        U_ -= nHat()*(nHat() & U_);
        U_.correctBoundaryConditions();
    }

    return tUEqn;
}


void kinematicSingleLayer::solveThickness
(
    const volScalarField& pu,
    const volScalarField& pp,
    const fvVectorMatrix& UEqn
)
{
    if (debug)
    {
        Info<< "kinematicSingleLayer::solveThickness()" << endl;
    }

    volScalarField rUA(1.0/UEqn.A());
    U_ = rUA*UEqn.H();

    surfaceScalarField deltarUAf(fvc::interpolate(delta_*rUA));
    surfaceScalarField rhof(fvc::interpolate(rho_));

    surfaceScalarField phiAdd
    (
        "phiAdd",
        regionMesh().magSf()
      * (
            fvc::snGrad(pu, "snGrad(p)")
          + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
        )
      - (fvc::interpolate(rho_*gTan()) & regionMesh().Sf())
    );
    constrainFilmField(phiAdd, 0.0);

    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(U_*rho_) & regionMesh().Sf())
      - deltarUAf*phiAdd*rhof
    );
    constrainFilmField(phid, 0.0);

    surfaceScalarField ddrhorUAppf
    (
        fvc::interpolate(delta_)*deltarUAf*rhof*fvc::interpolate(pp)
    );
//    constrainFilmField(ddrhorUAppf, 0.0);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        // Film thickness equation
        fvScalarMatrix deltaEqn
        (
            fvm::ddt(rho_, delta_)
          + fvm::div(phid, delta_)
          - fvm::laplacian(ddrhorUAppf, delta_)
         ==
            rhoSp_
        );

        deltaEqn.solve();

        if (nonOrth == nNonOrthCorr_)
        {
            phiAdd +=
                fvc::interpolate(pp)
              * fvc::snGrad(delta_)
              * regionMesh().magSf();

            phi_ == deltaEqn.flux();
        }
    }

    // Bound film thickness by a minimum of zero
    delta_.max(0.0);

    // Update U field
    U_ -= fvc::reconstruct(deltarUAf*phiAdd);

    // Remove any patch-normal components of velocity
    U_ -= nHat()*(nHat() & U_);

    U_.correctBoundaryConditions();

    // Continuity check
    continuityCheck();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kinematicSingleLayer::kinematicSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const bool readFields
)
:
    surfaceFilmModel(modelType, mesh, g),

    momentumPredictor_(solution().subDict("PISO").lookup("momentumPredictor")),
    nOuterCorr_(readLabel(solution().subDict("PISO").lookup("nOuterCorr"))),
    nCorr_(readLabel(solution().subDict("PISO").lookup("nCorr"))),
    nNonOrthCorr_(readLabel(solution().subDict("PISO").lookup("nNonOrthCorr"))),

    cumulativeContErr_(0.0),

    Cf_(readScalar(coeffs().lookup("Cf"))),
    deltaStable_(coeffs().lookup("deltaStable")),

    rho_
    (
        IOobject
        (
            "rhof",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimDensity, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    mu_
    (
        IOobject
        (
            "muf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimPressure*dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigma_
    (
        IOobject
        (
            "sigmaf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/sqr(dimTime), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    delta_
    (
        IOobject
        (
            "deltaf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    U_
    (
        IOobject
        (
            "Uf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Us_
    (
        IOobject
        (
            "Usf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_,
        zeroGradientFvPatchScalarField::typeName
    ),
    Uw_
    (
        IOobject
        (
            "Uwf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_,
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaRho_
    (
        IOobject
        (
            delta_.name() + "*" + rho_.name(),
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", delta_.dimensions()*rho_.dimensions(), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    phi_
    (
        IOobject
        (
            "phi",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimLength*dimMass/dimTime
    ),

    massForPrimary_
    (
        IOobject
        (
            "massForPrimary",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    diametersForPrimary_
    (
        IOobject
        (
            "diametersForPrimary",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimLength, -1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    massPhaseChangeForPrimary_
    (
        IOobject
        (
            "massPhaseChangeForPrimary",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    USp_
    (
        IOobject
        (
            "USpf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector
        (
            "zero", dimMass*dimVelocity/dimArea/dimTime, vector::zero
        ),
        this->mappedPushedFieldPatchTypes<vector>()
    ),
    pSp_
    (
        IOobject
        (
            "pSpf",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),
    rhoSp_
    (
        IOobject
        (
            "rhoSpf",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),

    USpPrimary_
    (
        IOobject
        (
            USp_.name(), // must have same name as USp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedVector("zero", USp_.dimensions(), vector::zero)
    ),
    pSpPrimary_
    (
        IOobject
        (
            pSp_.name(), // must have same name as pSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar("zero", pSp_.dimensions(), 0.0)
    ),
    rhoSpPrimary_
    (
        IOobject
        (
            rhoSp_.name(), // must have same name as rhoSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar("zero", rhoSp_.dimensions(), 0.0)
    ),

    UPrimary_
    (
        IOobject
        (
            "U", // must have same name as U to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector("zero", dimVelocity, vector::zero),
        this->mappedFieldAndInternalPatchTypes<vector>()
    ),
    pPrimary_
    (
        IOobject
        (
            "p", // must have same name as p to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    rhoPrimary_
    (
        IOobject
        (
            "rho", // must have same name as rho to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimDensity, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    muPrimary_
    (
        IOobject
        (
            "mu", // must have same name as mu to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimPressure*dimTime, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    injection_(injectionModel::New(*this, coeffs_)),

    addedMassTotal_(0.0),
    injectedMassTotal_(0.0)
{
    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctThermoFields();

        deltaRho_ == delta_*rho_;
        phi_ = fvc::interpolate(deltaRho_*U_) & regionMesh().Sf();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kinematicSingleLayer::~kinematicSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void kinematicSingleLayer::addSources
(
    const label patchI,
    const label faceI,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    if (debug)
    {
        Info<< "\nSurface film: " << type() << ": adding to film source:" << nl
            << "    mass     = " << massSource << nl
            << "    momentum = " << momentumSource << nl
            << "    pressure = " << pressureSource << endl;
    }

    rhoSpPrimary_.boundaryField()[patchI][faceI] += massSource;
    USpPrimary_.boundaryField()[patchI][faceI] += momentumSource;
    pSpPrimary_.boundaryField()[patchI][faceI] += pressureSource;

    addedMassTotal_ += massSource;
}


void kinematicSingleLayer::preEvolveRegion()
{
    transferPrimaryRegionThermoFields();

    correctThermoFields();

    transferPrimaryRegionSourceFields();
}


void kinematicSingleLayer::evolveRegion()
{
    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    // Implicit pressure source coefficient
    tmp<volScalarField> tpp(this->pp());

    for (int oCorr=0; oCorr<nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution - varies with delta_
        tmp<volScalarField> tpu(this->pu());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn());
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update film wall and surface velocities
    updateSurfaceVelocities();

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


scalar kinematicSingleLayer::CourantNumber() const
{
    scalar CoNum = 0.0;

    if (regionMesh().nInternalFaces() > 0)
    {
        const scalar deltaT = time_.deltaTValue();

        const surfaceScalarField SfUfbyDelta =
            regionMesh().surfaceInterpolation::deltaCoeffs()*mag(phi_);
        const surfaceScalarField rhoDelta = fvc::interpolate(rho_*delta_);
        const surfaceScalarField& magSf = regionMesh().magSf();

        forAll(rhoDelta, i)
        {
            if (rhoDelta[i] > ROOTVSMALL)
            {
                CoNum = max(CoNum, SfUfbyDelta[i]/rhoDelta[i]/magSf[i]*deltaT);
            }
        }
    }

    reduce(CoNum, maxOp<scalar>());

    Info<< "Film max Courant number: " << CoNum << endl;

    return CoNum;
}


const volVectorField& kinematicSingleLayer::U() const
{
    return U_;
}


const volVectorField& kinematicSingleLayer::Us() const
{
    return Us_;
}


const volVectorField& kinematicSingleLayer::Uw() const
{
    return Uw_;
}


const volScalarField& kinematicSingleLayer::rho() const
{
    return rho_;
}


const volScalarField& kinematicSingleLayer::T() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::T() const"
    )   << "T field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Ts() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::Ts() const"
    )   << "Ts field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Tw() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::Tw() const"
    )   << "Tw field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Cp() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::Cp() const"
    )   << "Cp field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::kappa() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::kappa() const"
    )   << "kappa field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::massForPrimary() const
{
    return massForPrimary_;
}


const volScalarField& kinematicSingleLayer::diametersForPrimary() const
{
    return diametersForPrimary_;
}


const volScalarField& kinematicSingleLayer::massPhaseChangeForPrimary() const
{
    return massPhaseChangeForPrimary_;
}


void kinematicSingleLayer::info() const
{
    Info<< "\nSurface film: " << type() << endl;

    Info<< indent << "added mass         = "
        << returnReduce<scalar>(addedMassTotal_, sumOp<scalar>()) << nl
        << indent << "current mass       = "
        << gSum((deltaRho_*magSf())()) << nl
        << indent << "injected mass      = "
        << returnReduce<scalar>(injectedMassTotal_, sumOp<scalar>()) << nl
        << indent << "min/max(mag(U))    = " << min(mag(U_)).value() << ", "
        << max(mag(U_)).value() << nl
        << indent << "min/max(delta)     = " << min(delta_).value() << ", "
        << max(delta_).value() << nl;
}


tmp<DimensionedField<scalar, volMesh> >  kinematicSingleLayer::Srho() const
{
    tmp<DimensionedField<scalar, volMesh> > tSrho
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "kinematicSingleLayer::Srho",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    scalarField& Srho = tSrho();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchI = intCoupledPatchIDs()[i];
        const mapDistribute& distMap = mappedPatches_[filmPatchI].map();

        scalarField patchMass =
            massPhaseChangeForPrimary_.boundaryField()[filmPatchI];
        distMap.distribute(patchMass);

        const label primaryPatchI = primaryPatchIDs()[i];
        const unallocLabelList& cells =
            primaryMesh().boundaryMesh()[primaryPatchI].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] = patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<DimensionedField<scalar, volMesh> > kinematicSingleLayer::Srho
(
    const label
) const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "kinematicSingleLayer::Srho(i)",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );
}


tmp<DimensionedField<scalar, volMesh> > kinematicSingleLayer::Sh() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "kinematicSingleLayer::Sh",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
