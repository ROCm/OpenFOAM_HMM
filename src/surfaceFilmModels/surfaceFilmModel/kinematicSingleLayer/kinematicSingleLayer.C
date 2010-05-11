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

#include "kinematicSingleLayer.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"
#include "directMappedWallPolyPatch.H"

// Sub-models
#include "injectionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(kinematicSingleLayer, 0);
        addToRunTimeSelectionTable
        (
            surfaceFilmModel,
            kinematicSingleLayer,
            mesh
        );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::surfaceFilmModels::kinematicSingleLayer::read()
{
    if (surfaceFilmModel::read())
    {
        const dictionary& solution = filmRegion_.solutionDict().subDict("PISO");
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


void Foam::surfaceFilmModels::kinematicSingleLayer::initialise()
{
    if (debug)
    {
        Pout<< "kinematicSingleLayer::initialise()" << endl;
    }

    label nBoundaryFaces = 0;
    DynamicList<label> primaryPatchIDs;
    DynamicList<label> filmBottomPatchIDs;
    const polyBoundaryMesh& bm = filmRegion_.boundaryMesh();
    forAll(bm, patchI)
    {
        const polyPatch& pp = bm[patchI];
        if (isA<directMappedWallPolyPatch>(pp))
        {
            if (debug)
            {
                Pout<< "found " << directMappedWallPolyPatch::typeName
                    <<  " " << pp.name() << endl;
            }

            filmBottomPatchIDs.append(patchI);
            const directMappedWallPolyPatch& dwpp =
                refCast<const directMappedWallPolyPatch>(pp);

            primaryPatchIDs.append
            (
                mesh_.boundaryMesh().findPatchID(dwpp.samplePatch())
            );

            const labelList& fCells = pp.faceCells();
            nBoundaryFaces += fCells.size();

            // Cache patch normals
            UIndirectList<vector>(nHat_, fCells) = pp.faceNormals();

            // Cache mesh face areas
            UIndirectList<scalar>(magSf_, fCells) = mag(pp.faceAreas());
        }
    }
    nHat_.correctBoundaryConditions();
    magSf_.correctBoundaryConditions();

    primaryPatchIDs_.transfer(primaryPatchIDs);
    filmBottomPatchIDs_.transfer(filmBottomPatchIDs);

    if (nBoundaryFaces == 0)
    {
        WarningIn("kinematicSingleLayer::initialise()")
            << "Film model being applied without direct mapped boundary "
            << "conditions" << endl;
    }

    if (nBoundaryFaces != filmRegion_.nCells())
    {
        FatalErrorIn("kinematicSingleLayer::initialise()")
            << "Number of primary region coupled boundary faces not equal to "
            << "the number of cells in the film region" << nl
            << abort(FatalError);
    }

    scalarField topMagSf(magSf_.size(), 0.0);
    filmTopPatchIDs_.setSize(filmBottomPatchIDs_.size(), -1);
    forAll(filmBottomPatchIDs_, i)
    {
        const label patchI = filmBottomPatchIDs_[i];
        const polyPatch& ppBottom = bm[patchI];
        if (ppBottom.size() > 0)
        {
            label cellId = bm[patchI].faceCells()[0];
            const cell& cFaces = filmRegion_.cells()[cellId];

            label faceBottom = ppBottom.start();
            label faceTop =
                 cFaces.opposingFaceLabel(faceBottom, filmRegion_.faces());

            label topPatchI = bm.whichPatch(faceTop);
            filmTopPatchIDs_[i] = topPatchI;
            const polyPatch& ppTop = bm[topPatchI];
            UIndirectList<scalar>(topMagSf, ppTop.faceCells()) =
                mag(ppTop.faceAreas());
        }
    }

    Pstream::listCombineGather(filmTopPatchIDs_, maxEqOp<label>());
    Pstream::listCombineScatter(filmTopPatchIDs_);

    magSf_.field() = 0.5*(magSf_ + topMagSf);
    magSf_.correctBoundaryConditions();
}


void Foam::surfaceFilmModels::kinematicSingleLayer::
resetPrimaryRegionSourceTerms()
{
    rhoSpPrimary_ == dimensionedScalar("zero", rhoSp_.dimensions(), 0.0);
    USpPrimary_ == dimensionedVector("zero", USp_.dimensions(), vector::zero);
    pSpPrimary_ == dimensionedScalar("zero", pSp_.dimensions(), 0.0);
}


void Foam::surfaceFilmModels::kinematicSingleLayer::
transferPrimaryRegionFields()
{
    // Update pressure and velocity from primary region via direct mapped
    // (coupled) boundary conditions
    UPrimary_.correctBoundaryConditions();
    pPrimary_.correctBoundaryConditions();

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    rhoSp_.correctBoundaryConditions();
    USp_.correctBoundaryConditions();
    pSp_.correctBoundaryConditions();

    // Convert accummulated source terms into per unit area per unit time
    // Note: boundary values will still have original (neat) values
    const scalar deltaT = filmRegion_.time().deltaTValue();
    rhoSp_.field() /= magSf_*deltaT;
    USp_.field() /= magSf_*deltaT;
    pSp_.field() /= magSf_*deltaT;
}


Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::kinematicSingleLayer::pu()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pu",
                filmRegion_.time().timeName(),
                filmRegion_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pPrimary_                  // pressure (mapped from primary region)
          + pSp_                           // accumulated particle impingement
          - fvc::laplacian(sigma_, delta_) // surface tension
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::kinematicSingleLayer::pp()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pp",
                filmRegion_.time().timeName(),
                filmRegion_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -rho_*gNormClipped() // hydrostatic effect only
        )
    );
}


void Foam::surfaceFilmModels::kinematicSingleLayer::correctDetachedFilm()
{
    massForPrimary_ == dimensionedScalar("zero", dimMass, 0.0);
    diametersForPrimary_ == dimensionedScalar("zero", dimLength, -1.0);

    const scalarField gNorm = this->gNorm();

    forAll(gNorm, i)
    {
        if (gNorm[i] > SMALL)
        {
            scalar ddelta = max(0.0, delta_[i] - deltaStable_.value());
            massForPrimary_[i] = ddelta*rho_[i]*magSf_[i];
        }
    }
}


void Foam::surfaceFilmModels::kinematicSingleLayer::updateSubmodels()
{
    correctDetachedFilm();

    // Update injection model - mass returned is actual mass injected
    injection_->correct(massForPrimary_, diametersForPrimary_);

    // Update cumulative detached mass counter
    detachedMass_ += sum(massForPrimary_.field());

    // Push values to boundaries ready for transfer to the primary region
    massForPrimary_.correctBoundaryConditions();
    diametersForPrimary_.correctBoundaryConditions();

    // Update source fields
    const dimensionedScalar deltaT = filmRegion_.time().deltaT();
    rhoSp_ -= massForPrimary_/magSf_/deltaT;
    USp_ -= massForPrimary_*U_/magSf_/deltaT;
}


Foam::scalar
Foam::surfaceFilmModels::kinematicSingleLayer::CourantNumber() const
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    if (filmRegion_.nInternalFaces())
    {
        const scalar deltaT = time_.deltaTValue();

        surfaceScalarField SfUfbyDelta =
            filmRegion_.surfaceInterpolation::deltaCoeffs()*mag(phi_)
           /fvc::interpolate
            (
                rho_*(delta_ + dimensionedScalar("SMALL", dimLength, SMALL))
            );

        CoNum = max(SfUfbyDelta/filmRegion_.magSf()).value()*deltaT;

        meanCoNum = (sum(SfUfbyDelta)/sum(filmRegion_.magSf())).value()*deltaT;
    }

    Info<< "    Courant number mean: " << meanCoNum << " max: " << CoNum
        << endl;

    return CoNum;
}


void Foam::surfaceFilmModels::kinematicSingleLayer::continuityCheck()
{
    const volScalarField deltaRho0 = deltaRho_;

    solveContinuity();

    if (debug)
    {
        volScalarField mass = deltaRho_*magSf_;
        dimensionedScalar totalMass =
            fvc::domainIntegrate(mass)
          + dimensionedScalar("SMALL", dimMass*dimVolume, ROOTVSMALL);

        scalar sumLocalContErr =
            (
                fvc::domainIntegrate(mag(mass - magSf_*deltaRho0))/totalMass
            ).value();

        scalar globalContErr =
            (
                fvc::domainIntegrate(mass - magSf_*deltaRho0)/totalMass
            ).value();

        cumulativeContErr_ += globalContErr;

        Info<< "Surface film: " << type() << nl
            << "    time step continuity errors: sum local = "
            << sumLocalContErr << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr_ << endl;
    }
}


void Foam::surfaceFilmModels::kinematicSingleLayer::solveContinuity()
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


Foam::tmp<Foam::fvVectorMatrix>
Foam::surfaceFilmModels::kinematicSingleLayer::tau
(
    volVectorField& U
) const
{
    DimensionedField<vector, volMesh> Uw
    (
        IOobject
        (
            "Uw",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        filmRegion_,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );

    forAll(filmBottomPatchIDs_, i)
    {
        label patchI = filmBottomPatchIDs_[i];
        const polyPatch& pp = filmRegion_.boundaryMesh()[patchI];
        UIndirectList<vector>(Uw, pp.faceCells()) =
            UPrimary_.boundaryField()[patchI];
    }
    Uw -= nHat_*(Uw & nHat_);

    volVectorField Us
    (
        IOobject
        (
            "Us",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        filmRegion_,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );

    // TODO: use wall functions
    Us.dimensionedInternalField() = UPrimary_.dimensionedInternalField();
    Us -= nHat_*(Us & nHat_);

    volScalarField Cs("Cs", rho_*Cf_*mag(Us - U_));
    volScalarField Cw
    (
        "Cw",
        mu_/(0.5*(delta_ + dimensionedScalar("SMALL", dimLength, SMALL)))
    );
    Cw.min(1.0e+06);

    return
    (
       - fvm::Sp(Cs, U) + Cs*Us
       - fvm::Sp(Cw, U) + Cw*Uw
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::surfaceFilmModels::kinematicSingleLayer::solveMomentum
(
    const volScalarField& pu,
    const volScalarField& pp
)
{
    if (debug)
    {
        Info<< "kinematicSingleLayer::solveMomentum()" << endl;
    }

    // Momentum
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(deltaRho_, U_)
      + fvm::div(phi_, U_)
     ==
        USp_
      + tau(U_)
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
                    filmRegion_.magSf()
                  * (
                        fvc::snGrad(pu, "snGrad(p)")
                      + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
                      + fvc::snGrad(delta_)*fvc::interpolate(pp)
                    )
                  - (fvc::interpolate(rho_*gTan()) & filmRegion_.Sf())
                )
            )
        );

        // Remove any patch-normal components of velocity
        U_ -= nHat_*(nHat_ & U_);
        U_.correctBoundaryConditions();
    }

    return tUEqn;
}


void Foam::surfaceFilmModels::kinematicSingleLayer::solveThickness
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

    volScalarField rUA = 1.0/UEqn.A();
    U_ = rUA*UEqn.H();

    surfaceScalarField deltarUAf = fvc::interpolate(delta_*rUA);
    surfaceScalarField rhof = fvc::interpolate(rho_);

    surfaceScalarField phiAdd
    (
        "phiAdd",
        filmRegion_.magSf()
      * (
            fvc::snGrad(pu, "snGrad(p)")
          + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
        )
      - (fvc::interpolate(rho_*gTan()) & filmRegion_.Sf())
    );
    constrainFilmField(phiAdd, 0.0);

    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(U_*rho_) & filmRegion_.Sf())
      - deltarUAf*phiAdd*rhof
    );
    constrainFilmField(phid, 0.0);

    surfaceScalarField ddrhorUAppf =
        fvc::interpolate(delta_)*deltarUAf*rhof*fvc::interpolate(pp);
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
              * filmRegion_.magSf();

            phi_ == deltaEqn.flux();
        }
    }

    // Bound film thickness by a minimum of zero
    delta_.max(0.0);

    // Update U field
    U_ -= fvc::reconstruct(deltarUAf*phiAdd);

    // Remove any patch-normal components of velocity
    U_ -= nHat_*(nHat_ & U_);

    U_.correctBoundaryConditions();

    // Continuity check
    continuityCheck();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::kinematicSingleLayer::kinematicSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    surfaceFilmModel(modelType, mesh, g),
    filmRegion_
    (
        IOobject
        (
            filmRegionName_,
            time_.timeName(),
            time_,
            IOobject::MUST_READ
        )
    ),
    nHat_
    (
        IOobject
        (
            "nHat",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedVector("zero", dimless, vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    magSf_
    (
        IOobject
        (
            "magSf",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar("zero", dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    primaryPatchIDs_(0),
    filmTopPatchIDs_(0),
    filmBottomPatchIDs_(0),

    momentumPredictor_(true),
    nOuterCorr_(-1),
    nCorr_(-1),
    nNonOrthCorr_(-1),
    cumulativeContErr_(0.0),

    Cf_(0.0),
    deltaStable_("deltaStable", dimLength, 0.0),

    rho_
    (
        IOobject
        (
            "rhof",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar(coeffs_.lookup("rho"))
    ),
    mu_
    (
        IOobject
        (
            "muf",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar(coeffs_.lookup("mu"))
    ),
    sigma_
    (
        IOobject
        (
            "sigmaf",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar(coeffs_.lookup("sigma"))
    ),

    delta_
    (
        IOobject
        (
            "deltaf",
            time_.timeName(),
            filmRegion_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_
    ),
    U_
    (
        IOobject
        (
            "Uf",
            time_.timeName(),
            filmRegion_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_
    ),
    deltaRho_
    (
        IOobject
        (
            delta_.name() + "*" + rho_.name(),
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        delta_*rho_,
        zeroGradientFvPatchScalarField::typeName
    ),

    phi_
    (
        IOobject
        (
            "phi",
            time_.timeName(),
            filmRegion_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(delta_*rho_*U_) & filmRegion_.Sf()
    ),

    massForPrimary_
    (
        IOobject
        (
            "massForPrimary",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar("zero", dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    diametersForPrimary_
    (
        IOobject
        (
            "diametersForPrimary",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar("zero", dimLength, -1.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    USp_
    (
        IOobject
        (
            "USpf",
            time_.timeName(),
            filmRegion_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_
    ),
    pSp_
    (
        IOobject
        (
            "pSpf",
            time_.timeName(),
            filmRegion_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_
    ),
    rhoSp_
    (
        IOobject
        (
            "rhoSpf",
            time_.timeName(),
            filmRegion_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmRegion_,
        dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0),
        pSp_.boundaryField().types()
    ),

    USpPrimary_
    (
        IOobject
        (
            USp_.name(), // must have same name as USp_ to enable mapping
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", USp_.dimensions(), vector::zero)
    ),
    pSpPrimary_
    (
        IOobject
        (
            pSp_.name(), // must have same name as pSp_ to enable mapping
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", pSp_.dimensions(), 0.0)
    ),
    rhoSpPrimary_
    (
        IOobject
        (
            rhoSp_.name(), // must have same name as rhoSp_ to enable mapping
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", rhoSp_.dimensions(), 0.0)
    ),

    UPrimary_
    (
        IOobject
        (
            "U", // must have same name as U on primary region to enable mapping
            time_.timeName(),
            filmRegion_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_
    ),
    pPrimary_
    (
        IOobject
        (
            "p", // must have same name as p on primary region to enable mapping
            time_.timeName(),
            filmRegion_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        filmRegion_
    ),

    injection_(injectionModel::New(*this, coeffs_)),

    addedMass_(0.0),
    detachedMass_(0.0)
{
    read();

    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::kinematicSingleLayer::~kinematicSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::fvMesh&
Foam::surfaceFilmModels::kinematicSingleLayer::film() const
{
    return filmRegion_;
}


const Foam::labelList&
Foam::surfaceFilmModels::kinematicSingleLayer::filmBottomPatchIDs() const
{
    return filmBottomPatchIDs_;
}


const Foam::labelList&
Foam::surfaceFilmModels::kinematicSingleLayer::filmTopPatchIDs() const
{
    return filmTopPatchIDs_;
}


const Foam::labelList&
Foam::surfaceFilmModels::kinematicSingleLayer::primaryPatchIDs() const
{
    return primaryPatchIDs_;
}


bool Foam::surfaceFilmModels::kinematicSingleLayer::isFilmPatch
(
    const label patchI
) const
{
    if (!active_)
    {
        return false;
    }

    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchI)
        {
            return true;
        }
    }

    return false;
}


void Foam::surfaceFilmModels::kinematicSingleLayer::addSources
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

    addedMass_ += massSource;
}


void Foam::surfaceFilmModels::kinematicSingleLayer::evolveFilm()
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
            // Solve thickness for delta_
            solveThickness(pu(), pp(), UEqn());
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


const Foam::volVectorField&
Foam::surfaceFilmModels::kinematicSingleLayer::U() const
{
    return U_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::rho() const
{
    return rho_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::T() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::T() const"
    )   << "T field not available for " << type() << abort(FatalError);

    return reinterpret_cast<const volScalarField&>(null);
}


const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::cp() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::cp() const"
    )   << "cp field not available for " << type() << abort(FatalError);

    return reinterpret_cast<const volScalarField&>(null);
}



const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::massForPrimary() const
{
    return massForPrimary_;
}


const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::diametersForPrimary() const
{
    return diametersForPrimary_;
}


void Foam::surfaceFilmModels::kinematicSingleLayer::info() const
{
    Info<< "\nSurface film: " << type() << endl;

    // Output Courant number for info only - does not change time step
    CourantNumber();

    Info<< indent << "added mass      = "
        << returnReduce<scalar>(addedMass_, sumOp<scalar>()) << nl
        << indent << "current mass    = "
        << gSum((deltaRho_*magSf_)()) << nl
        << indent << "detached mass   = "
        << returnReduce<scalar>(detachedMass_, sumOp<scalar>()) << nl
        << indent << "min/max(mag(U)) = " << min(mag(U_)).value() << ", "
        << max(mag(U_)).value() << nl
        << indent << "min/max(delta)  = " << min(delta_).value() << ", "
        << max(delta_).value() << nl;
}


// ************************************************************************* //
