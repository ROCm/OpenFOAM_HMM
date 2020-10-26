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

#include "contactAngleForce.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "unitConversion.H"
#include "meshWavePatchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(contactAngleForce, 0);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void contactAngleForce::initialise()
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactAngleForce::contactAngleForce
(
    const word& typeName,
    liquidFilmBase& film,
    const dictionary& dict
)
:
    force(typeName, film, dict),
    Ccf_(coeffDict_.get<scalar>("Ccf")),
    mask_
    (
        IOobject
        (
            typeName + ":fContactForceMask",
            film.primaryMesh().time().timeName(),
            film.primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("mask", dimless, 1.0)
    )
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

contactAngleForce::~contactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<faVectorMatrix> contactAngleForce::correct(areaVectorField& U)
{
    tmp<areaVectorField> tForce
    (
        new areaVectorField
        (
            IOobject
            (
                typeName + ":contactForce",
                film().primaryMesh().time().timeName(),
                film().primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedVector(dimForce/dimDensity/dimArea, Zero)
        )
    );

    vectorField& force = tForce.ref().primitiveFieldRef();

    const labelUList& own = film().regionMesh().owner();
    const labelUList& nbr = film().regionMesh().neighbour();

    const DimensionedField<scalar, areaMesh>& magSf = film().regionMesh().S();

    tmp<areaScalarField> talpha = film().alpha();
    const areaScalarField& sigma = film().sigma();

    const areaScalarField& rhof = film().rho();

    tmp<areaScalarField> ttheta = theta();
    const areaScalarField& theta = ttheta();

    const areaVectorField gradAlpha(fac::grad(talpha()));

    forAll(nbr, edgei)
    {
        const label faceO = own[edgei];
        const label faceN = nbr[edgei];

        label facei = -1;
        if ((talpha()[faceO] > 0.5) && (talpha()[faceN] < 0.5))
        {
            facei = faceO;
        }
        else if ((talpha()[faceO] < 0.5) && (talpha()[faceN] > 0.5))
        {
            facei = faceN;
        }

        if (facei != -1 && mask_[facei] > 0.5)
        {
            const scalar invDx = film().regionMesh().deltaCoeffs()[edgei];
            const vector n
            (
                gradAlpha[facei]/(mag(gradAlpha[facei]) + ROOTVSMALL)
            );
            const scalar cosTheta = cos(degToRad(theta[facei]));

            force[facei] +=
                Ccf_*n*sigma[facei]*(1 - cosTheta)/invDx/rhof[facei];
        }
    }

    forAll(sigma.boundaryField(), patchi)
    {
        const faPatchField<scalar>& alphaPf = sigma.boundaryField()[patchi];
        const faPatchField<scalar>& sigmaPf = sigma.boundaryField()[patchi];
        const labelUList& faces = alphaPf.patch().edgeFaces();

        forAll(sigmaPf, edgei)
        {
            label face0 = faces[edgei];
            force[face0] = vector::zero;
        }
    }


    force /= magSf.field();

    if (film().regionMesh().time().writeTime())
    {
        tForce().write();
        gradAlpha.write();
    }

    tmp<faVectorMatrix> tfvm
    (
        new faVectorMatrix(U, dimForce/dimDensity)
    );

    tfvm.ref() += tForce;

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
