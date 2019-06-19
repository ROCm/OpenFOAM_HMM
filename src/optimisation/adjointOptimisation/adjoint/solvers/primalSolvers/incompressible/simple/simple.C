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

#include "simple.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simple, 0);
    addToRunTimeSelectionTable(incompressiblePrimalSolver, simple, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::simple::addExtraSchemes()
{
    if (vars_().useSolverNameForFields())
    {
        WarningInFunction
            << "useSolverNameForFields is set to true for primalSolver "
            << solverName() << nl << tab
            << "Appending variable names with the solver name" << nl << tab
            << "Please adjust the necessary entries in fvSchemes and fvSolution"
            << nl << endl;
    }
}


void Foam::simple::continuityErrors()
{
    surfaceScalarField& phi = vars_().phiInst();
    volScalarField contErr(fvc::div(phi));

    scalar sumLocalContErr = mesh_.time().deltaTValue()*
        mag(contErr)().weightedAverage(mesh_.V()).value();

    scalar globalContErr = mesh_.time().deltaTValue()*
        contErr.weightedAverage(mesh_.V()).value();
    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simple::simple
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict
)
:
    incompressiblePrimalSolver(mesh, managerType, dict),
    solverControl_(SIMPLEControl::New(mesh, managerType, *this)),
    MRF_(mesh),
    cumulativeContErr_(Zero),
    objectives_(0)
{
    vars_.reset(new incompressibleVars(mesh, solverControl_()));
    addExtraSchemes();
    setRefCell
    (
        vars_().pInst(),
        solverControl_().dict(),
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::simple::readDict(const dictionary& dict)
{
    if (incompressiblePrimalSolver::readDict(dict))
    {
        return true;
    }

    return false;
}

void Foam::simple::solveIter()
{
    const Time& time = mesh_.time();
    Info<< "Time = " << time.timeName() << "\n" << endl;

    // Grab references
    incompressibleVars& vars = vars_();
    volScalarField& p = vars.pInst();
    volVectorField& U = vars.UInst();
    surfaceScalarField& phi = vars.phiInst();
    autoPtr<incompressible::turbulenceModel>& turbulence = vars.turbulence();
    label&  pRefCell  = solverControl_().pRefCell();
    scalar& pRefValue = solverControl_().pRefValue();

    // Momentum predictor
    //~~~~~~~~~~~~~~~~~~~

    MRF_.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF_.DDt(U)
      + turbulence->divDevReff(U)
      ==
        fvOptions_(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    addOptimisationTypeSource(UEqn);

    UEqn.relax();

    fvOptions_.constrain(UEqn);

    if (solverControl_().momentumPredictor())
    {
        Foam::solve(UEqn == -fvc::grad(p));

        fvOptions_.correct(U);
    }

    // Pressure Eq
    //~~~~~~~~~~~~
    {
        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        MRF_.makeRelative(phiHbyA);
        adjustPhi(phiHbyA, U, p);

        tmp<volScalarField> rAtU(rAU);

        if (solverControl_().consistent())
        {
            rAtU = 1.0/(1.0/rAU - UEqn.H1());
            phiHbyA +=
                fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh_.magSf();
            HbyA -= (rAU - rAtU())*fvc::grad(p);
        }

        tUEqn.clear();

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAtU(), MRF_);

        // Non-orthogonal pressure corrector loop
        while (solverControl_().correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
            );

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve();

            if (solverControl_().finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
            }
        }

        continuityErrors();

        // Explicitly relax pressure for momentum corrector
        p.relax();

        // Momentum corrector
        U = HbyA - rAtU()*fvc::grad(p);
        U.correctBoundaryConditions();
        fvOptions_.correct(U);
    }

    vars_().laminarTransport().correct();
    turbulence->correct();

    solverControl_().write();

    // Print objective values to screen and compute mean value
    Info<< endl;
    for (objective* obj : objectives_)
    {
        Info<< obj->objectiveName() << " : " << obj->J() << endl;
        obj->accumulateJMean(solverControl_());
        obj->writeInstantaneousValue();
    }

    // Average fields if necessary
    vars.computeMeanFields();

    // Print execution time
    time.printExecutionTime(Info);
}


void Foam::simple::solve()
{
    // Iterate
    if (active_)
    {
        // Get the objectives for this solver
        if (objectives_.empty())
        {
            objectives_ = getObjectiveFunctions();
        }

        // Reset initial and mean fields before solving
        restoreInitValues();
        vars_().resetMeanFields();

        // Validate turbulence model fields
        vars_().turbulence()->validate();

        while (solverControl_().loop())
        {
            solveIter();
        }

        // Safety
        objectives_.clear();
    }
}


bool Foam::simple::loop()
{
    return solverControl_().loop();
}


void Foam::simple::restoreInitValues()
{
    vars_().restoreInitValues();
}


bool Foam::simple::writeData(Ostream& os) const
{
    os.writeEntry("averageIter", solverControl_().averageIter());

    return true;
}


// ************************************************************************* //
