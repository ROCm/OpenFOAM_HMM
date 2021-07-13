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

#include "liquidFilmModel.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "gravityMeshObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liquidFilmModel, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool liquidFilmModel::read(const dictionary& dict)
{
    liquidFilmBase::read(dict);
    return true;
}

void liquidFilmModel::correctThermoFields()
{
    scalarField X(thermo_.size(), 1);

    forAll (rho_, faceI)
    {
        rho_[faceI] = thermo_.rho(pRef_, Tf_[faceI], X);
        mu_[faceI] = thermo_.mu(pRef_, Tf_[faceI], X);
        sigma_[faceI] = thermo_.sigma(pRef_, Tf_[faceI], X);
    }

    forAll (regionMesh().boundary(), patchI)
    {
        const scalarField& patchTf = Tf_.boundaryFieldRef()[patchI];

        scalarField& patchRho = rho_.boundaryFieldRef()[patchI];
        scalarField& patchmu = mu_.boundaryFieldRef()[patchI];
        scalarField& patchsigma = sigma_.boundaryFieldRef()[patchI];

        forAll(patchRho, edgeI)
        {
            patchRho[edgeI] = thermo_.rho(pRef_, patchTf[edgeI], X);
            patchmu[edgeI] = thermo_.mu(pRef_, patchTf[edgeI], X);
            patchsigma[edgeI] = thermo_.sigma(pRef_, patchTf[edgeI], X);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

liquidFilmModel::liquidFilmModel
(
    const word& modelType,
    const fvPatch& patch,
    const dictionary& dict
)
:
    liquidFilmBase(modelType, patch, dict),

    thermo_(dict.subDict("thermo")),

    rho_
    (
        IOobject
        (
            "rhof",
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity, Zero)
    ),
    mu_
    (
        IOobject
        (
            "muf",
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimViscosity, Zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigmaf",
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/sqr(dimTime), Zero)
    ),
    hRho_
    (
        IOobject
        (
            h_.name() + "*" + rho_.name(),
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(h_.dimensions()*rho_.dimensions(), Zero)
    ),
    rhoSp_
    (
        IOobject
        (
            "rhoSp",
            primaryMesh().time().timeName(),
            primaryMesh()
        ),
        regionMesh(),
        dimensionedScalar(dimVelocity, Zero)
    ),
    USp_
    (
        IOobject
        (
            "USp",
            primaryMesh().time().timeName(),
            primaryMesh()
        ),
        regionMesh(),
        dimensionedVector(sqr(dimVelocity), Zero)
    ),
    pnSp_
    (
        IOobject
        (
            "pnSp",
            primaryMesh().time().timeName(),
            primaryMesh()
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    turbulence_(filmTurbulenceModel::New(*this, dict))
{
    correctThermoFields();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liquidFilmModel::~liquidFilmModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const areaScalarField& liquidFilmModel::mu() const
{
     return mu_;
}


const areaScalarField& liquidFilmModel::rho() const
{
     return rho_;
}


const areaScalarField& liquidFilmModel::sigma() const
{
     return sigma_;
}

void liquidFilmModel::preEvolveRegion()
{
    const scalar deltaT = primaryMesh().time().deltaTValue();
    const scalarField rAreaDeltaT = 1/deltaT/regionMesh().S().field();

    // Map the total mass, mom and pnSource from particles
    rhoSp_.primitiveFieldRef() =
        vsm().mapToSurface(massSource_.boundaryField()[patchID()]);
    // [kg.m/s]
    USp_.primitiveFieldRef() =
        vsm().mapToSurface(momentumSource_.boundaryField()[patchID()]);
    pnSp_.primitiveFieldRef() =
        vsm().mapToSurface(pnSource_.boundaryField()[patchID()]);


    // Calculate rate per area
    rhoSp_.primitiveFieldRef() *= rAreaDeltaT/rho_;
    USp_.primitiveFieldRef() *= rAreaDeltaT/rho_;
    pnSp_.primitiveFieldRef() *= rAreaDeltaT/rho_;
}

void liquidFilmModel::postEvolveRegion()
{
    massSource_ == dimensionedScalar(massSource_.dimensions(), Zero);
    momentumSource_ == dimensionedVector(momentumSource_.dimensions(), Zero);
    pnSource_ == dimensionedScalar(pnSource_.dimensions(), Zero);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
