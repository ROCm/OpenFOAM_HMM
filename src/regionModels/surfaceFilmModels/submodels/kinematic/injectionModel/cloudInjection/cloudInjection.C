/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "cloudInjection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cloudInjection, 0);
addToRunTimeSelectionTable(injectionModel, cloudInjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cloudInjection::cloudInjection
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    injectionModel(type(), owner, dict),
    particlesPerParcel_(readScalar(coeffs_.lookup("particlesPerParcel"))),
    rndGen_(label(0), -1),
    parcelDistribution_
    (
        distributionModels::distributionModel::New
        (
            coeffs_.subDict("parcelDistribution"),
            rndGen_
        )
    ),
    diameter_(owner.regionMesh().nCells(), 0.0)
{
    forAll(diameter_, faceI)
    {
        diameter_[faceI] = parcelDistribution_->sample();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cloudInjection::~cloudInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void cloudInjection::correct
(
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const scalar pi = constant::mathematical::pi;
    const scalarField& rhoFilm = owner().rho();

    // Collect the data to be transferred
    forAll(massToInject, cellI)
    {
        scalar rho = rhoFilm[cellI];
        scalar diam = diameter_[cellI];
        scalar minMass = particlesPerParcel_*rho*pi/6*pow3(diam);

        if (massToInject[cellI] > minMass)
        {
            // All mass can be injected - set particle diameter
            diameterToInject[cellI] = diameter_[cellI];

            // Retrieve new particle diameter sample
            diameter_[cellI] = parcelDistribution_->sample();
        }
        else
        {
            // Mass below minimum threshold - cannot be injected
            massToInject[cellI] = 0.0;
            diameterToInject[cellI] = -1.0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
