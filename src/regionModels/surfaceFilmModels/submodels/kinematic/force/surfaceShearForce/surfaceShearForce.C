/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
#include "turbulenceModel.H"

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
    // local reference to film model
    const kinematicSingleLayer& film =
        static_cast<const kinematicSingleLayer&>(owner_);

    // local references to film fields
    const volScalarField& mu = film.mu();
    const volVectorField& Uw = film.Uw();
    const volScalarField& delta = film.delta();
    const volVectorField& Up = film.UPrimary();

    // film surface linear coeff to apply to velocity
    tmp<volScalarField> tCs;

    typedef compressible::turbulenceModel turbModel;
    if (film.primaryMesh().foundObject<turbModel>("turbulenceModel"))
    {
        // local reference to turbulence model
        const turbModel& turb =
            film.primaryMesh().lookupObject<turbModel>("turbulenceModel");

        // calculate and store the stress on the primary region
        const volSymmTensorField primaryReff(turb.devRhoReff());

        // create stress field on film
        // - note boundary condition types (mapped)
        // - to map, the field name must be the same as the field on the
        //   primary region
        volSymmTensorField Reff
        (
            IOobject
            (
                primaryReff.name(),
                film.regionMesh().time().timeName(),
                film.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film.regionMesh(),
            dimensionedSymmTensor
            (
                "zero",
                primaryReff.dimensions(),
                symmTensor::zero
            ),
            film.mappedFieldAndInternalPatchTypes<symmTensor>()
        );

        // map stress from primary region to film region
        Reff.correctBoundaryConditions();

        dimensionedScalar U0("SMALL", U.dimensions(), SMALL);
        volVectorField UHat("UHat", (Up - U)/(mag(Up - U) + U0));

        // shear stress tangential to the film
        volVectorField tauTan
        (
            "tauTan",
            UHat & (Reff + film.nHat()*(-film.nHat() & Reff))
        );

        // note: Cf_ 'should' be 1 in this case
        tCs = Cf_*mag(tauTan)/(mag(Up - U) + U0);
    }
    else
    {
        // laminar case - employ simple coeff-based model
        const volScalarField& rhop = film.rhoPrimary();
        tCs = Cf_*rhop*mag(Up - U);
    }

    dimensionedScalar d0("SMALL", delta.dimensions(), SMALL);

    // linear coeffs to apply to velocity
    const volScalarField& Cs = tCs();
    volScalarField Cw("Cw", mu/(0.3333*(delta + d0)));
    Cw.min(5000.0);

    return
    (
       - fvm::Sp(Cs, U) + Cs*Up // surface contribution
       - fvm::Sp(Cw, U) + Cw*Uw // wall contribution
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
