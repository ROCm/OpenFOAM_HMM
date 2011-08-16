/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "surfaceShearForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceShearForce, 0);
addToRunTimeSelectionTable(force, surfaceShearForce, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceShearForce::surfaceShearForce
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    force(typeName, owner, dict),
    Cf_(readScalar(coeffs_.lookup("Cf")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

surfaceShearForce::~surfaceShearForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> surfaceShearForce::correct(volVectorField& U)
{
    const kinematicSingleLayer& film =
        static_cast<const kinematicSingleLayer&>(owner_);

    const volScalarField& rho = film.rho();
    const volScalarField& mu = film.mu();
    const volVectorField& Us = film.Us();
    const volVectorField& Uw = film.Uw();
    const volScalarField& delta = film.delta();

    // Calculate shear stress
    volScalarField Cs("Cs", rho*Cf_*mag(Us - U));
    volScalarField Cw
    (
        "Cw",
        mu/(0.3333*(delta + dimensionedScalar("SMALL", dimLength, SMALL)))
    );
    Cw.min(1.0e+06);

    return
    (
       - fvm::Sp(Cs, U) + Cs*Us // surface contribution
       - fvm::Sp(Cw, U) + Cw*Uw // wall contribution
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
