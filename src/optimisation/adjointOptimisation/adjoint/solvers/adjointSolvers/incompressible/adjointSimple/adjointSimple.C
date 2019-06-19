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

#include "adjointSimple.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
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

void Foam::adjointSimple::addExtraSchemes()
{
    if (getAdjointVars().useSolverNameForFields())
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
    const surfaceScalarField& phia = getAdjointVars().phiaInst();
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
    cumulativeContErr_(Zero),
    adjointSensitivity_(nullptr)
{
    adjointVars_.reset
    (
        new incompressibleAdjointVars
        (
            mesh,
            solverControl_(),
            objectiveManagerPtr_(),
            primalVars_
        )
    ),
    ATCModel_.reset
    (
        ATCModel::New
        (
            mesh,
            primalVars_,
            getAdjointVars(),
            dict.subDict("ATCModel")
        ).ptr()
    );

    addExtraSchemes();
    setRefCell
    (
        getAdjointVars().paInst(),
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
                getAdjointVars(),
                objectiveManagerPtr_(),
                fvOptionsAdjoint_
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adjointSimple::readDict(const dictionary& dict)
{
    if (incompressibleAdjointSolver::readDict(dict))
    {
        if (adjointSensitivity_.valid())
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
    const Time& time = mesh_.time();
    Info<< "Time = " << time.timeName() << "\n" << endl;

    // Grab primal references
    const surfaceScalarField& phi = primalVars_.phi();
    // Grab adjoint references
    incompressibleAdjointVars& adjointVars = getAdjointVars();
    volScalarField& pa = adjointVars.paInst();
    volVectorField& Ua = adjointVars.UaInst();
    surfaceScalarField& phia = adjointVars.phiaInst();
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars.adjointTurbulence();
    const label&  paRefCell  = solverControl_().pRefCell();
    const scalar& paRefValue = solverControl_().pRefValue();

    // Momentum predictor
    //~~~~~~~~~~~~~~~~~~~

    tmp<fvVectorMatrix> tUaEqn
    (
        fvm::div(-phi, Ua)
      + adjointTurbulence->divDevReff(Ua)
      + adjointTurbulence->adjointMeanFlowSource()
      ==
        fvOptionsAdjoint_(Ua)
    );
    fvVectorMatrix& UaEqn = tUaEqn.ref();

    // Add sources from boundary conditions
    UaEqn.boundaryManipulate(Ua.boundaryFieldRef());

    // Add sources from volume-based objectives
    objectiveManagerPtr_().addUaEqnSource(UaEqn);

    // Add ATC term
    ATCModel_->addATC(UaEqn);

    // Add source from optimisationType (e.g. topology)
    addOptimisationTypeSource(UaEqn);

    UaEqn.relax();

    fvOptionsAdjoint_.constrain(UaEqn);

    if (solverControl_().momentumPredictor())
    {
        Foam::solve(UaEqn == -fvc::grad(pa));

        fvOptionsAdjoint_.correct(Ua);
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
        fvOptionsAdjoint_.correct(Ua);
        pa.correctBoundaryConditions();
    }

    adjointTurbulence->correct();

    if (solverControl_().printMaxMags())
    {
        dimensionedScalar maxUa = max(mag(Ua));
        dimensionedScalar maxpa = max(mag(pa));
        Info<< "Max mag of adjoint velocity = " << maxUa.value() << endl;
        Info<< "Max mag of adjoint pressure = " << maxpa.value() << endl;
    }

    solverControl_().write();

    // Average fields if necessary
    adjointVars.computeMeanFields();

    // Print execution time
    time.printExecutionTime(Info);
}


void Foam::adjointSimple::solve()
{
    if (active_)
    {
        // Reset mean fields before solving
        getAdjointVars().resetMeanFields();

        // Iterate
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


void Foam::adjointSimple::computeObjectiveSensitivities()
{
    if (computeSensitivities_)
    {
        sensitivities_.reset
        (
            new scalarField(adjointSensitivity_().calculateSensitivities())
        );
    }
    else
    {
        sensitivities_.reset(new scalarField(0));
    }
}


const Foam::scalarField& Foam::adjointSimple::getObjectiveSensitivities()
{
    if (!sensitivities_.valid())
    {
        computeObjectiveSensitivities();
    }

    return sensitivities_();
}


Foam::sensitivity& Foam::adjointSimple::getSensitivityBase()
{
    if (adjointSensitivity_.valid())
    {
        return adjointSensitivity_();
    }
    else
    {
        FatalErrorInFunction
            << "Sensitivity object not allocated \n"
            << "Turn computeSensitivities on in "
            << solverName_
            << nl << nl
            << exit(FatalError);

        return adjointSensitivity_();
    }
}


bool Foam::adjointSimple::writeData(Ostream& os) const
{
    os.writeEntry("averageIter", solverControl_().averageIter());

    return true;
}


// ************************************************************************* //
