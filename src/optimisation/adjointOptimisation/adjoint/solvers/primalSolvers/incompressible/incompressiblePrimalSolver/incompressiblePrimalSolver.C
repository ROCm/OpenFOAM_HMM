/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "incompressiblePrimalSolver.H"
#include "adjustPhi.H"
#include "adjointSolver.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressiblePrimalSolver, 0);
    defineRunTimeSelectionTable(incompressiblePrimalSolver, dictionary);
    addToRunTimeSelectionTable
    (
        primalSolver,
        incompressiblePrimalSolver,
        primalSolver
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressiblePrimalSolver::incompressiblePrimalSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict
)
:
    primalSolver(mesh, managerType, dict),
    vars_(nullptr),
    phiReconstructionTol_
    (
        dict.subOrEmptyDict("fieldReconstruction").
            lookupOrDefault<scalar>("tolerance", scalar(5.e-5))
    ),
    phiReconstructionIters_
    (
        dict.subOrEmptyDict("fieldReconstruction").
            lookupOrDefault<label>("iters", label(10))
    ),
    fvOptions_
    (
        mesh_,
        dict.subOrEmptyDict("fvOptions")
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::incompressiblePrimalSolver>
Foam::incompressiblePrimalSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict
)
{
    const word solverType(dict.get<word>("solver"));
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(solverType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown incompressiblePrimalSolver type " << solverType
            << nl << nl
            << "Valid incompressiblePrimalSolver types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        autoPtr<incompressiblePrimalSolver>
        (
            cstrIter()(mesh, managerType, dict)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::incompressiblePrimalSolver::readDict(const dictionary& dict)
{
    if (primalSolver::readDict(dict))
    {
        fvOptions_.read(dict.subOrEmptyDict("fvOptions"));

        return true;
    }

    return false;
}


Foam::List<Foam::objective*>
Foam::incompressiblePrimalSolver::getObjectiveFunctions() const
{
    DynamicList<objective*> objectives(10);

    auto adjointSolvers = mesh_.lookupClass<adjointSolver>();

    for (adjointSolver* adjointPtr : adjointSolvers)
    {
        adjointSolver& adjoint = *adjointPtr;

        if (adjoint.active() && adjoint.primalSolverName() == solverName_)
        {
            PtrList<objective>& managerObjectives =
                adjoint.getObjectiveManager().getObjectiveFunctions();

            for (objective& obj : managerObjectives)
            {
                objectives.append(&obj);
            }
        }
    }

    return objectives;
}


bool Foam::incompressiblePrimalSolver::useSolverNameForFields() const
{
    return vars_().useSolverNameForFields();
}


const Foam::incompressibleVars&
Foam::incompressiblePrimalSolver::getVars() const
{
    return vars_();
}


Foam::incompressibleVars&
Foam::incompressiblePrimalSolver::getVars()
{
    return vars_.ref();
}


void Foam::incompressiblePrimalSolver::correctBoundaryConditions()
{
    incompressibleVars& vars = vars_();
    // Update boundary conditions for all primal volFields,
    // including averaged ones, if present
    vars.correctBoundaryConditions();

    // phi cannot be updated through correctBoundayrConditions.
    // Re-compute based on the Rhie-Chow interpolation scheme.
    // This is a non-linear process
    // (phi depends on UEqn().A() which depends on phi)
    // so iterations are required

    volScalarField& p = vars.p();
    volVectorField& U = vars.U();
    surfaceScalarField& phi = vars.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence = vars.turbulence();

    scalar contError(GREAT), diff(GREAT);
    for (label iter = 0; iter < phiReconstructionIters_; ++iter)
    {
        Info<< "phi correction iteration " << iter << endl;

        // Form momentum equations to grab A
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, U)
          + turbulence->divDevReff(U)
          ==
            fvOptions_(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();
        UEqn.relax();
        fvOptions_.constrain(UEqn);

        // Pressure equation will give the Rhie-Chow correction
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA("HbyA", U);
        HbyA = rAU*UEqn.H();
        tUEqn.clear();

        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        adjustPhi(phiHbyA, U, p);

        //fvOptions_.makeRelative(phiHbyA);

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );
        phi = phiHbyA - pEqn.flux();

        // Check convergence
        // ~~~~~~~~~~~~~~~~~
        scalar contErrorNew =
            mesh_.time().deltaTValue()*
            mag(fvc::div(phi)())().weightedAverage(mesh_.V()).value();

        Info<< "sum local = " << contErrorNew << endl;
        diff = mag(contErrorNew - contError)/contError;
        contError = contErrorNew;

        if (diff < phiReconstructionTol_) break;

        Info<< endl;
     }
}


// ************************************************************************* //
