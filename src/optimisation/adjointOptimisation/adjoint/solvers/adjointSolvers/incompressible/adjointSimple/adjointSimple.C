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

#include "adjointSimple.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointSimple, 0);
    addToRunTimeSelectionTable
    (
        incompressibleAdjointSolver,
        adjointSimple,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::incompressibleAdjointVars& Foam::adjointSimple::allocateVars()
{
    vars_.reset
    (
        new incompressibleAdjointVars
        (
            mesh_,
            solverControl_(),
            objectiveManagerPtr_(),
            primalVars_
        )
    );
    return getAdjointVars();
}


void Foam::adjointSimple::addExtraSchemes()
{
    if (adjointVars_.useSolverNameForFields())
    {
        WarningInFunction
            << "useSolverNameForFields is set to true for adjointSolver "
            << solverName() << nl << tab
            << "Appending variable names with the solver name" << nl << tab
            << "Please adjust the necessary entries in fvSchemes and fvSolution"
            << nl << endl;
    }
}


void Foam::adjointSimple::continuityErrors()
{
    const surfaceScalarField& phia = adjointVars_.phiaInst();
    volScalarField contErr(fvc::div(phia));

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

Foam::adjointSimple::adjointSimple
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName
)
:
    incompressibleAdjointSolver(mesh, managerType, dict, primalSolverName),
    solverControl_(SIMPLEControl::New(mesh, managerType, *this)),
    adjointVars_(allocateVars()),
    cumulativeContErr_(Zero),
    adjointSensitivity_(nullptr)
{
    ATCModel_.reset
    (
        ATCModel::New
        (
            mesh,
            primalVars_,
            adjointVars_,
            dict.subDict("ATCModel")
        ).ptr()
    );

    addExtraSchemes();
    setRefCell
    (
        adjointVars_.paInst(),
        solverControl_().dict(),
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );

    if (computeSensitivities_)
    {
        const IOdictionary& optDict =
            mesh.lookupObject<IOdictionary>("optimisationDict");

        adjointSensitivity_.reset
        (
            incompressible::adjointSensitivity::New
            (
                mesh,
                optDict.subDict("optimisation").subDict("sensitivities"),
                primalVars_,
                adjointVars_,
                objectiveManagerPtr_()
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adjointSimple::readDict(const dictionary& dict)
{
    if (incompressibleAdjointSolver::readDict(dict))
    {
        if (adjointSensitivity_)
        {
            const IOdictionary& optDict =
                mesh_.lookupObject<IOdictionary>("optimisationDict");

            adjointSensitivity_().readDict
            (
                optDict.subDict("optimisation").subDict("sensitivities")
            );
        }

        return true;
    }

    return false;
}


void Foam::adjointSimple::solveIter()
{
    preIter();
    mainIter();
    postIter();
}


void Foam::adjointSimple::preIter()
{
    Info<< "Time = " << mesh_.time().timeName() << "\n" << endl;
}


void Foam::adjointSimple::mainIter()
{
    // Grab primal references
    const surfaceScalarField& phi = primalVars_.phi();
    // Grab adjoint references
    volScalarField& pa = adjointVars_.paInst();
    volVectorField& Ua = adjointVars_.UaInst();
    surfaceScalarField& phia = adjointVars_.phiaInst();
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars_.adjointTurbulence();
    const label&  paRefCell  = solverControl_().pRefCell();
    const scalar& paRefValue = solverControl_().pRefValue();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Momentum predictor
    //~~~~~~~~~~~~~~~~~~~

    tmp<fvVectorMatrix> tUaEqn
    (
        fvm::div(-phi, Ua)
      + adjointTurbulence->divDevReff(Ua)
      + adjointTurbulence->adjointMeanFlowSource()
      ==
        fvOptions(Ua)
    );
    fvVectorMatrix& UaEqn = tUaEqn.ref();

    // Add sources from boundary conditions
    UaEqn.boundaryManipulate(Ua.boundaryFieldRef());

    // Add sources from volume-based objectives
    objectiveManagerPtr_().addUaEqnSource(UaEqn);

    // Add ATC term
    ATCModel_->addATC(UaEqn);

    // Additional source terms (e.g. energy equation)
    addMomentumSource(UaEqn);

    UaEqn.relax();

    fvOptions.constrain(UaEqn);

    if (solverControl_().momentumPredictor())
    {
        Foam::solve(UaEqn == -fvc::grad(pa));

        fvOptions.correct(Ua);
    }

    // Pressure Eq
    //~~~~~~~~~~~~
    {
        volScalarField rAUa(1.0/UaEqn.A());
        // 190402: Vag: to be updated.
        // Probably a constrainHabyA by class is needed?
        volVectorField HabyA(constrainHbyA(rAUa*UaEqn.H(), Ua, pa));
        surfaceScalarField phiaHbyA("phiaHbyA", fvc::flux(HabyA));
        adjustPhi(phiaHbyA, Ua, pa);

        tmp<volScalarField> rAtUa(rAUa);

        if (solverControl_().consistent())
        {
            rAtUa = 1.0/(1.0/rAUa - UaEqn.H1());
            phiaHbyA +=
                fvc::interpolate(rAtUa() - rAUa)*fvc::snGrad(pa)*mesh_.magSf();
            HabyA -= (rAUa - rAtUa())*fvc::grad(pa);
        }

        tUaEqn.clear();

        // Update the pressure BCs to ensure flux consistency
        // constrainPressure(p, U, phiHbyA, rAtU(), MRF_);

        // Non-orthogonal pressure corrector loop
        while (solverControl_().correctNonOrthogonal())
        {
            fvScalarMatrix paEqn
            (
                fvm::laplacian(rAtUa(), pa) == fvc::div(phiaHbyA)
            );

            paEqn.boundaryManipulate(pa.boundaryFieldRef());

            addPressureSource(paEqn);

            fvOptions.constrain(paEqn);
            paEqn.setReference(paRefCell, paRefValue);

            paEqn.solve();

            if (solverControl_().finalNonOrthogonalIter())
            {
                phia = phiaHbyA - paEqn.flux();
            }
        }

        continuityErrors();

        // Explicitly relax pressure for adjoint momentum corrector
        pa.relax();

        // Momentum corrector
        Ua = HabyA - rAtUa()*fvc::grad(pa);
        Ua.correctBoundaryConditions();
        fvOptions.correct(Ua);
        pa.correctBoundaryConditions();
    }

    adjointTurbulence->correct();

    if (solverControl_().printMaxMags())
    {
        dimensionedScalar maxUa = gMax(mag(Ua)());
        dimensionedScalar maxpa = gMax(mag(pa)());
        Info<< "Max mag of adjoint velocity = " << maxUa.value() << endl;
        Info<< "Max mag of adjoint pressure = " << maxpa.value() << endl;
    }
}


void Foam::adjointSimple::postIter()
{
    solverControl_().write();

    // Average fields if necessary
    adjointVars_.computeMeanFields();

    // Print execution time
    mesh_.time().printExecutionTime(Info);
}


void Foam::adjointSimple::solve()
{
    if (active_)
    {
        preLoop();
        while (solverControl_().loop())
        {
            solveIter();
        }
    }
}


bool Foam::adjointSimple::loop()
{
    return solverControl_().loop();
}


void Foam::adjointSimple::preLoop()
{
    // Reset mean fields before solving
    adjointVars_.resetMeanFields();
}


void Foam::adjointSimple::computeObjectiveSensitivities()
{
    if (computeSensitivities_)
    {
        adjointSensitivity_->accumulateIntegrand(scalar(1));
        const scalarField& sens = adjointSensitivity_->calculateSensitivities();
        if (!sensitivities_)
        {
            sensitivities_.reset(new scalarField(sens.size(), Zero));
        }
        *sensitivities_ = sens;
    }
    else
    {
        sensitivities_.reset(new scalarField());
    }
}


const Foam::scalarField& Foam::adjointSimple::getObjectiveSensitivities()
{
    if (!sensitivities_)
    {
        computeObjectiveSensitivities();
    }

    return sensitivities_();
}


void Foam::adjointSimple::clearSensitivities()
{
    if (computeSensitivities_)
    {
        adjointSensitivity_->clearSensitivities();
        adjointSolver::clearSensitivities();
    }
}


Foam::sensitivity& Foam::adjointSimple::getSensitivityBase()
{
    if (!adjointSensitivity_.valid())
    {
        FatalErrorInFunction
            << "Sensitivity object not allocated" << nl
            << "Turn computeSensitivities on in "
            << solverName_
            << nl << nl
            << exit(FatalError);
    }

    return adjointSensitivity_();
}


void Foam::adjointSimple::addMomentumSource(fvVectorMatrix& matrix)
{
    // Does nothing
}


void Foam::adjointSimple::addPressureSource(fvScalarMatrix& matrix)
{
    // Does nothing
}


void Foam::adjointSimple::updatePrimalBasedQuantities()
{
    incompressibleAdjointSolver::updatePrimalBasedQuantities();

    // Update objective function related quantities
    objectiveManagerPtr_->updateAndWrite();
}


bool Foam::adjointSimple::writeData(Ostream& os) const
{
    os.writeEntry("averageIter", solverControl_().averageIter());

    return true;
}


// ************************************************************************* //
