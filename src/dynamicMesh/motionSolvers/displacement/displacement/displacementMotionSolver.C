/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementMotionSolver, 0);
    defineRunTimeSelectionTable(displacementMotionSolver, displacement);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementMotionSolver::displacementMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const word& type
)
:
    points0MotionSolver(mesh, dict, type),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement",
            time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    )
{}


Foam::displacementMotionSolver::displacementMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0,
    const word& type
)
:
    points0MotionSolver(mesh, dict, points0, type),
    pointDisplacement_
    (
        IOobject(pointDisplacement, "pointDisplacement"),
        pointDisplacement
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::displacementMotionSolver>
Foam::displacementMotionSolver::New
(
    const word& solverTypeName,
    const polyMesh& mesh,
    const IOdictionary& solverDict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
{
    //const word solverTypeName(solverDict.lookup("solver"));

    Info<< "Selecting motion solver: " << solverTypeName << endl;

    const_cast<Time&>(mesh.time()).libs().open
    (
        solverDict,
        "motionSolverLibs",
        displacementConstructorTablePtr_
    );

    if (!displacementConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "solver table is empty"
            << exit(FatalError);
    }

    auto cstrIter = displacementConstructorTablePtr_->cfind(solverTypeName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown solver type "
            << solverTypeName << nl << nl
            << "Valid solver types :" << endl
            << displacementConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<displacementMotionSolver>
    (
        cstrIter()
        (
            mesh,
            solverDict,
            pointDisplacement,
            points0
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementMotionSolver::~displacementMotionSolver()
{}


// ************************************************************************* //
