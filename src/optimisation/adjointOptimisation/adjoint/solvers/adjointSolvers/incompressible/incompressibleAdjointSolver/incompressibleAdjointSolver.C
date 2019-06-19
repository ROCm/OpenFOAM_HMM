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

#include "incompressibleAdjointSolver.H"
#include "incompressiblePrimalSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleAdjointSolver, 0);
    defineRunTimeSelectionTable(incompressibleAdjointSolver, dictionary);
    addToRunTimeSelectionTable
    (
        adjointSolver,
        incompressibleAdjointSolver,
        adjointSolver
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleAdjointSolver::incompressibleAdjointSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName
)
:
    adjointSolver(mesh, managerType, dict, primalSolverName),
    primalVars_
    (
        mesh.lookupObjectRef<incompressiblePrimalSolver>(primalSolverName).
            getVars()
    ),
    adjointVars_(nullptr),
    ATCModel_(nullptr),
    fvOptionsAdjoint_
    (
        mesh_,
        dict.subOrEmptyDict("fvOptions")
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::incompressibleAdjointSolver>
Foam::incompressibleAdjointSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName
)
{
    const word solverType(dict.get<word>("solver"));
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(solverType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown incompressibleAdjointSolver type "
            << solverType << nl << nl
            << "Valid incompressibleAdjointSolver types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        autoPtr<incompressibleAdjointSolver>
        (
            cstrIter()(mesh, managerType, dict, primalSolverName)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::incompressibleAdjointSolver::readDict(const dictionary& dict)
{
    if (adjointSolver::readDict(dict))
    {
        fvOptionsAdjoint_.read(dict.subOrEmptyDict("fvOptions"));

        return true;
    }

    return false;
}

bool Foam::incompressibleAdjointSolver::useSolverNameForFields() const
{
    return getAdjointVars().useSolverNameForFields();
}


const Foam::incompressibleVars&
Foam::incompressibleAdjointSolver::getPrimalVars() const
{
    return primalVars_;
}



const Foam::incompressibleAdjointVars&
Foam::incompressibleAdjointSolver::getAdjointVars() const
{
    return adjointVars_();
}


Foam::incompressibleAdjointVars&
Foam::incompressibleAdjointSolver::getAdjointVars()
{
    return adjointVars_.ref();
}



const Foam::autoPtr<Foam::ATCModel>&
Foam::incompressibleAdjointSolver::getATCModel() const
{
    return ATCModel_;
}


Foam::autoPtr<Foam::ATCModel>& Foam::incompressibleAdjointSolver::getATCModel()
{
    return ATCModel_;
}


Foam::fv::optionAdjointList&
Foam::incompressibleAdjointSolver::getFvOptionsAdjoint()
{
    return fvOptionsAdjoint_;
}


void Foam::incompressibleAdjointSolver::updatePrimalBasedQuantities()
{
    getAdjointVars().adjointTurbulence()->setChangedPrimalSolution();
}


// ************************************************************************* //
