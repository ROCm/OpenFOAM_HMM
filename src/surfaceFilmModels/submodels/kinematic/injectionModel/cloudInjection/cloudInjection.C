/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "cloudInjection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "mathematicalConstants.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(cloudInjection, 0);
        addToRunTimeSelectionTable(injectionModel, cloudInjection, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::cloudInjection::cloudInjection
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    injectionModel(type(), owner, dict),
    particlesPerParcel_(readScalar(coeffs_.lookup("particlesPerParcel"))),
    rndGen_(label(0), -1),
    parcelPDF_(pdfs::pdf::New(coeffs_.subDict("parcelPDF"), rndGen_)),
    diameter_(owner.film().nCells(), 0.0)
{
    forAll(diameter_, faceI)
    {
        diameter_[faceI] = parcelPDF_->sample();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::cloudInjection::~cloudInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::cloudInjection::inject
(
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    correctDetachedFilm(massToInject);

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
            diameter_[cellI] = parcelPDF_->sample();
        }
        else
        {
            // Mass below minimum threshold - cannot be injected
            massToInject[cellI] = 0.0;
            diameterToInject[cellI] = -1.0;
        }
    }
}


// ************************************************************************* //
