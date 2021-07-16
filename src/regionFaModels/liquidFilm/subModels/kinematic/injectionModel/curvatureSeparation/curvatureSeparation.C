/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "curvatureSeparation.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "stringListOps.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(curvatureSeparation, 0);
addToRunTimeSelectionTable
(
    injectionModel,
    curvatureSeparation,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<areaScalarField> curvatureSeparation::calcInvR1
(
    const areaVectorField& U,
    const scalarField& calcCosAngle
) const
{
    const dimensionedScalar smallU(dimVelocity, ROOTVSMALL);
    const areaVectorField UHat(U/(mag(U) + smallU));
    tmp<areaScalarField> tinvR1
    (
        new areaScalarField("invR1", UHat & (UHat & -gradNHat_))
    );

    scalarField& invR1 = tinvR1.ref().primitiveFieldRef();

    // apply defined patch radii
    const scalar rMin = 1e-6;
    const scalar definedInvR1 = 1.0/max(rMin, definedPatchRadii_);

    if (definedPatchRadii_ > 0)
    {
        invR1 = definedInvR1;
    }

    // filter out large radii
    const scalar rMax = 1e6;
    forAll(invR1, i)
    {
        if ((mag(invR1[i]) < 1/rMax))
        {
            invR1[i] = -1.0;
        }
    }

    return tinvR1;
}


tmp<scalarField> curvatureSeparation::calcCosAngle
(
    const edgeScalarField& phi
) const
{
    const areaVectorField& U = film().Uf();
    const dimensionedScalar smallU(dimVelocity, ROOTVSMALL);
    const areaVectorField UHat(U/(mag(U) + smallU));

    const faMesh& mesh = film().regionMesh();
    const labelUList& own = mesh.edgeOwner();
    const labelUList& nbr = mesh.edgeNeighbour();

    scalarField phiMax(mesh.nFaces(), -GREAT);
    scalarField cosAngle(UHat.size(), Zero);

    const scalarField invR1(calcInvR1(U, cosAngle));

    forAll(nbr, edgei)
    {
        const label cellO = own[edgei];
        const label cellN = nbr[edgei];

        if (phi[edgei] > phiMax[cellO])
        {
            phiMax[cellO] = phi[edgei];
            cosAngle[cellO] = -gHat_ & UHat[cellN];
        }
        if (-phi[edgei] > phiMax[cellN])
        {
            phiMax[cellN] = -phi[edgei];
            cosAngle[cellN] = -gHat_ & UHat[cellO];
        }
    }

    cosAngle *= pos(invR1);

    // checks
    if (debug && mesh.time().writeTime())
    {
        areaScalarField volCosAngle
        (
            IOobject
            (
                "cosAngle",
                film().primaryMesh().time().timeName(),
                film().primaryMesh(),
                IOobject::NO_READ
            ),
            film().regionMesh(),
            dimensionedScalar(dimless, Zero)
        );
        volCosAngle.primitiveFieldRef() = cosAngle;
        volCosAngle.correctBoundaryConditions();
        volCosAngle.write();
    }

    return max(min(cosAngle, scalar(1)), scalar(-1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

curvatureSeparation::curvatureSeparation
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    gradNHat_(fac::grad(film.regionMesh().faceAreaNormals())),
    deltaByR1Min_(coeffDict_.getOrDefault<scalar>("deltaByR1Min", 0)),
    definedPatchRadii_
    (
        coeffDict_.getOrDefault<scalar>("definedPatchRadii", 0)
    ),
    magG_(mag(film.g().value())),
    gHat_(Zero),
    fThreshold_
    (
        coeffDict_.getOrDefault<scalar>("fThreshold", 1e-8)
    ),
    minInvR1_
    (
        coeffDict_.getOrDefault<scalar>("minInvR1", 5)
    )
{
    if (magG_ < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Acceleration due to gravity must be non-zero"
            << exit(FatalError);
    }

    gHat_ = film.g().value()/magG_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void curvatureSeparation::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const faMesh& mesh = film().regionMesh();

    const areaScalarField& delta = film().h();
    const areaVectorField& U = film().Uf();
    const edgeScalarField& phi = film().phi2s();
    const areaScalarField& rho = film().rho();
    const scalarField magSqrU(magSqr(film().Uf()));
    const areaScalarField& sigma = film().sigma();

    const scalarField cosAngle(calcCosAngle(phi));
    const scalarField invR1(calcInvR1(U, cosAngle));


    // calculate force balance
    scalarField Fnet(mesh.nFaces(), Zero);
    scalarField separated(mesh.nFaces(), Zero);

    forAll(invR1, i)
    {
        if ((invR1[i] > minInvR1_) && (delta[i]*invR1[i] > deltaByR1Min_))
        {
            const scalar R1 = 1.0/(invR1[i] + ROOTVSMALL);
            const scalar R2 = R1 + delta[i];

            // inertial force
            const scalar Fi = -delta[i]*rho[i]*magSqrU[i]*72.0/60.0*invR1[i];

            // body force
            const scalar Fb =
               - 0.5*rho[i]*magG_*invR1[i]*(sqr(R1) - sqr(R2))*cosAngle[i];

            // surface force
            const scalar Fs = sigma[i]/R2;

            Fnet[i] = Fi + Fb + Fs;

            if (Fnet[i] + fThreshold_ < 0)
            {
                separated[i] = 1.0;
            }
        }
    }

    // inject all available mass
    massToInject = separated*availableMass;
    diameterToInject = separated*delta;
    availableMass -= separated*availableMass;

    addToInjectedMass(sum(massToInject));

    if (debug && mesh.time().writeTime())
    {
        areaScalarField volFnet
        (
            IOobject
            (
                "Fnet",
                film().primaryMesh().time().timeName(),
                film().primaryMesh(),
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar(dimForce, Zero)
        );
        volFnet.primitiveFieldRef() = Fnet;
        volFnet.write();

        areaScalarField areaSeparated
        (
            IOobject
            (
                "separated",
                film().primaryMesh().time().timeName(),
                film().primaryMesh(),
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar(dimMass, Zero)
        );
        areaSeparated.primitiveFieldRef() = separated;
        areaSeparated.write();

        areaScalarField areaMassToInject
        (
            IOobject
            (
                "massToInject",
                film().primaryMesh().time().timeName(),
                film().primaryMesh(),
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar(dimMass, Zero)
        );
        areaMassToInject.primitiveFieldRef() = massToInject;
        areaMassToInject.write();

        areaScalarField areaInvR1
        (
            IOobject
            (
                "InvR1",
                film().primaryMesh().time().timeName(),
                film().primaryMesh(),
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar(inv(dimLength), Zero)
        );
        areaInvR1.primitiveFieldRef() = invR1;
        areaInvR1.write();
    }

    injectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
