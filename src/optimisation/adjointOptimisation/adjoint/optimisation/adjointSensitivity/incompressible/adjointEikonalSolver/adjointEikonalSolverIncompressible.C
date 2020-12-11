/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "adjointEikonalSolverIncompressible.H"
#include "wallFvPatch.H"
#include "patchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointEikonalSolver, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

wordList adjointEikonalSolver::patchTypes() const
{
    wordList daTypes
    (
        mesh_.boundary().size(),
        fixedValueFvPatchScalarField::typeName
    );

    for (const label patchi : wallPatchIDs_)
    {
        daTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }

    return daTypes;
}


void adjointEikonalSolver::read()
{
    nEikonalIters_ = dict_.getOrDefault<label>("iters", 1000);
    tolerance_ = dict_.getOrDefault<scalar>("tolerance", 1e-6);
    epsilon_ = dict_.getOrDefault<scalar>("epsilon", 0.1);
}


tmp<surfaceScalarField> adjointEikonalSolver::computeYPhi()
{
    // Primal distance field
    const volScalarField& d = RASModelVars_().d();

    volVectorField ny
    (
        IOobject
        (
            "ny",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        patchDistMethod::patchTypes<vector>(mesh_, wallPatchIDs_)
    );

    const fvPatchList& patches = mesh_.boundary();
    volVectorField::Boundary& nybf = ny.boundaryFieldRef();

    for (const label patchi : wallPatchIDs_)
    {
        nybf[patchi] == -patches[patchi].nf();
    }

    ny = fvc::grad(d);

    surfaceVectorField nf(fvc::interpolate(ny));

    return tmp<surfaceScalarField>::New("yPhi", mesh_.Sf() & nf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointEikonalSolver::adjointEikonalSolver
(
    const fvMesh& mesh,
    const dictionary& dict,
    const autoPtr<incompressible::RASModelVariables>& RASModelVars,
    incompressibleAdjointVars& adjointVars,
    const labelHashSet& sensitivityPatchIDs
)
:
    mesh_(mesh),
    dict_(dict.subOrEmptyDict("adjointEikonalSolver")),
    RASModelVars_(RASModelVars),
    adjointTurbulence_(adjointVars.adjointTurbulence()),
    sensitivityPatchIDs_(sensitivityPatchIDs),
    nEikonalIters_(-1),
    tolerance_(-1),
    epsilon_(Zero),
    wallPatchIDs_(mesh_.boundaryMesh().findPatchIDs<wallPolyPatch>()),
    da_
    (
        IOobject
        (
            word
            (
                adjointVars.useSolverNameForFields() ?
                "da" + adjointTurbulence_().adjointSolverName() :
                "da"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(sqr(dimLength)/pow3(dimTime), Zero),
        patchTypes()
    ),
    source_
    (
        IOobject
        (
            "sourceEikonal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength/pow3(dimTime), Zero)
    ),
    distanceSensPtr_(createZeroBoundaryPtr<vector>(mesh_))
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool adjointEikonalSolver::readDict(const dictionary& dict)
{
    dict_ = dict.subOrEmptyDict("adjointEikonalSolver");

    return true;
}


void adjointEikonalSolver::accumulateIntegrand(const scalar dt)
{
    // Accumulate integrand from the current time step
    source_ += adjointTurbulence_->distanceSensitivities()*dt;
}


void adjointEikonalSolver::solve()
{
    read();

    // Primal distance field
    const volScalarField& d = RASModelVars_().d();

    // Convecting flux
    tmp<surfaceScalarField> tyPhi = computeYPhi();
    const surfaceScalarField& yPhi = tyPhi();

    // Iterate the adjoint to the eikonal equation
    for (label iter = 0; iter < nEikonalIters_; ++iter)
    {
        read();

        Info<< "Adjoint Eikonal Iteration : " << iter << endl;

        fvScalarMatrix daEqn
        (
            2*fvm::div(-yPhi, da_)
          + fvm::SuSp(-epsilon_*fvc::laplacian(d), da_)
          - epsilon_*fvm::laplacian(d, da_)
          + source_
        );

        daEqn.relax();
        scalar residual = daEqn.solve().initialResidual();
        Info<< "Max da " << gMax(mag(da_)()) << endl;

        mesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance_)
        {
            Info<< "\n***Reached adjoint eikonal convergence limit, iteration "
                << iter << "***\n\n";
            break;
        }
    }
    if (debug)
    {
        da_.write();
    }
}


void adjointEikonalSolver::reset()
{
    source_ == dimensionedScalar(source_.dimensions(), Zero);
    distanceSensPtr_() = vector::zero;
}


boundaryVectorField& adjointEikonalSolver::distanceSensitivities()
{
    Info<< "Calculating distance sensitivities " << endl;

    boundaryVectorField& distanceSens = distanceSensPtr_();

    const volScalarField& d = RASModelVars_().d();
    for (const label patchi : sensitivityPatchIDs_)
    {
        vectorField nf(mesh_.boundary()[patchi].nf());

        // No surface area included. Will be done by the actual sensitivity tool
        distanceSens[patchi] =
           -2.*da_.boundaryField()[patchi]
           *d.boundaryField()[patchi].snGrad()
           *d.boundaryField()[patchi].snGrad()*nf;
    }
    return distanceSens;
}


tmp<volTensorField> adjointEikonalSolver::getFISensitivityTerm()  const
{
    Info<< "Calculating distance sensitivities " << endl;

    const volScalarField& d = RASModelVars_().d();
    const volVectorField gradD(fvc::grad(d));

    volVectorField gradDDa
    (
        IOobject
        (
            "gradDDa",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(d.dimensions()*da_.dimensions()/dimLength, Zero),
        patchDistMethod::patchTypes<vector>(mesh_, wallPatchIDs_)
    );
    gradDDa = fvc::grad(d*da_);

    tmp<volTensorField> tdistanceSens
    (
        new volTensorField
        (
            IOobject
            (
                "distanceSensFI",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedTensor(da_.dimensions(), Zero)
        )
    );
    volTensorField& distanceSens = tdistanceSens.ref();

    distanceSens =
        - 2.*da_*gradD*gradD
        - epsilon_*gradD*gradDDa
        + epsilon_*da_*d*fvc::grad(gradD);

    return tdistanceSens;
}


const volScalarField& adjointEikonalSolver::da()
{
    return da_;
}


tmp<volVectorField> adjointEikonalSolver::gradEikonal()
{
    const volScalarField& d = RASModelVars_().d();
    volVectorField gradD(fvc::grad(d));
    return tmp<volVectorField>::New("gradEikonal", 2*gradD & fvc::grad(gradD));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
