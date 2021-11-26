/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "perturbedTemperatureDependentContactAngleForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(perturbedTemperatureDependentContactAngleForce, 0);
addToRunTimeSelectionTable
(
    force,
    perturbedTemperatureDependentContactAngleForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

perturbedTemperatureDependentContactAngleForce::
perturbedTemperatureDependentContactAngleForce
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    thetaPtr_(Function1<scalar>::New("theta", coeffDict_, &film.primaryMesh())),
    rndGen_(label(0)),
    distribution_
    (
        distributionModel::New
        (
            coeffDict_.subDict("distribution"),
            rndGen_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

perturbedTemperatureDependentContactAngleForce::
~perturbedTemperatureDependentContactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<areaScalarField> perturbedTemperatureDependentContactAngleForce::
theta() const
{
    tmp<areaScalarField> ttheta
    (
        new areaScalarField
        (
            IOobject
            (
                typeName + ":theta",
                film().primaryMesh().time().timeName(),
                film().primaryMesh()
            ),
            film().regionMesh(),
            dimensionedScalar(dimless, Zero)
        )
    );

    areaScalarField& theta = ttheta.ref();
    scalarField& thetai = theta.ref();

    const areaScalarField& T = film().Tf();

    // Initialize with the function of temperature
    thetai = thetaPtr_->value(T());

    // Add the stochastic perturbation
    forAll(thetai, facei)
    {
        thetai[facei] += distribution_->sample();
    }

    return ttheta;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
