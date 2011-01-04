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

#include "injectionModel.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(injectionModel, 0);
        defineRunTimeSelectionTable(injectionModel, dictionary);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::surfaceFilmModels::injectionModel::correctDetachedFilm
(
    scalarField& mass
) const
{
    mass = 0.0;

    const kinematicSingleLayer& film =
        refCast<const kinematicSingleLayer>(owner_);

    const scalarField gNorm(film.gNorm());
    const scalarField& delta = film.delta();
    const scalarField& rho = film.rho();
    const scalarField& magSf = film.magSf();
    const scalarField& massPhaseChange = film.massPhaseChangeForPrimary();

    forAll(gNorm, i)
    {
        if (gNorm[i] > SMALL)
        {
            const scalar ddelta = max(0.0, delta[i] - deltaStable_);
            mass[i] = max(0.0, ddelta*rho[i]*magSf[i] - massPhaseChange[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::injectionModel::injectionModel
(
    const surfaceFilmModel& owner
)
:
    owner_(owner),
    coeffs_(dictionary::null),
    injectedMass_(0.0)
{}


Foam::surfaceFilmModels::injectionModel::injectionModel
(
    const word& type,
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    owner_(owner),
    coeffs_(dict.subDict(type + "Coeffs")),
    injectedMass_(0.0),
    deltaStable_(readScalar(coeffs_.lookup("deltaStable")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::injectionModel::~injectionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::injectionModel::correct
(
    volScalarField& massToInject,
    volScalarField& diameterToInject
)
{
    inject(massToInject, diameterToInject);
    massToInject.correctBoundaryConditions();
    diameterToInject.correctBoundaryConditions();

    injectedMass_ += sum(massToInject.internalField());
}


void Foam::surfaceFilmModels::injectionModel::info() const
{
    Info<< indent << "injected mass      = "
        << returnReduce<scalar>(injectedMass_, sumOp<scalar>()) << nl;
}


// ************************************************************************* //
