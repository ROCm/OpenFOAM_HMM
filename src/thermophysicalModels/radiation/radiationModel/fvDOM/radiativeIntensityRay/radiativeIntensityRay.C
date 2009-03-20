/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "radiativeIntensityRay.H"
#include "fvm.H"
#include "fvDOM.H"

#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"
#include "radiationModel.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::radiation::radiativeIntensityRay::rayId = 0;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiativeIntensityRay::radiativeIntensityRay
(
    scalar& phii,
    scalar& thetai,
    scalar& deltaPhi,
    scalar& deltaTheta,
    label& lambdaj,
    const fvMesh& mesh,
    const absorptionEmissionModel& absEmmModel,
    const blackBodyEmission& blackBody
)
:
    absEmmModel_(absEmmModel),
    mesh_(mesh),
    blackBody_(blackBody),
    I_
    (
        IOobject
        (
            "I" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("I", dimMass/pow3(dimTime), 0.0)
    ),
    Qri_
    (
        IOobject
        (
            "Qr" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    )
{
    init(phii,thetai,deltaPhi,deltaTheta,lambdaj);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::radiativeIntensityRay::~radiativeIntensityRay()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::radiativeIntensityRay::init
(
    const scalar phii,
    const scalar thetai,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const scalar lambdaj
)
{
    phii_ = phii;
    thetai_ = thetai;
    nLambdaj_ = lambdaj;

    scalar sinTheta = Foam::sin(thetai);
    scalar cosTheta = Foam::cos(thetai);
    scalar sinPhi = Foam::sin(phii);
    scalar cosPhi = Foam::cos(phii);
    Si_ = vector(sinTheta*sinPhi, sinTheta*cosPhi, cosTheta);
    omegai_ = 2.0*Foam::sin(thetai)*Foam::sin(deltaTheta/2.0)*deltaPhi;
    Ilambdaj_.setSize(nLambdaj_);
    Di_ = vector
    (
        sinPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*thetai)
       *Foam::sin(deltaTheta)),
        cosPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*thetai)
       *Foam::sin(deltaTheta)),
        0.5*deltaPhi*Foam::sin(2.0*thetai)*Foam::sin(deltaTheta)
    );

    forAll(Ilambdaj_, i)
    {
        IOobject header
        (
            "Ilambda_" + name(rayId) + "_" + name(i),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            Ilambdaj_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Ilambda_" + name(rayId) + "_" + name(i),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_
                )
            );
        }
        else
        {
            volScalarField Idefault
            (
                IOobject
                (
                    "Idefault",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

            Ilambdaj_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Ilambda_" + name(rayId) + "_"+ name(i),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Idefault
                )
            );
        }
    }
    rayId++;
}


Foam::scalar Foam::radiation::radiativeIntensityRay::correct
(
    fvDOM* DomPtr
)
{
    Qri_ =  dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);

    scalar maxResidual = 0.0;

    for(label i = 0; i < nLambdaj_; i++)
    {
        volScalarField k = DomPtr->aj(i);

        volScalarField E = absEmmModel_.ECont(i)/Foam::mathematicalConstant::pi;

        surfaceScalarField Ji = Di_ & mesh_.Sf();

        volScalarField Ib = blackBody_.bj(i)/Foam::mathematicalConstant::pi;

        fvScalarMatrix IiEq
        (
            fvm::div(Ji, Ilambdaj_[i], " div(Ji,Ii_h)")
          + fvm::Sp(k*omegai_, Ilambdaj_[i])
         ==
            k*omegai_*Ib + E
        );

        IiEq.relax();

        scalar eqnResidual = solve
        (
            IiEq,
            mesh_.solver("Ii")
        ).initialResidual();

        maxResidual = max(eqnResidual, maxResidual);

    }
    return maxResidual;
}


void Foam::radiation::radiativeIntensityRay::addIntensity()
{
    I_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);

    for (label i = 0; i < nLambdaj_; i++)
    {
        I_ += absEmmModel_.addRadInt(i, Ilambdaj_[i]);
    }
}


void Foam::radiation::radiativeIntensityRay::add
(
    const scalarField& qr,
    const label patchI
) const
{
    Qri_.boundaryField()[patchI] += qr;
}


// ************************************************************************* //
