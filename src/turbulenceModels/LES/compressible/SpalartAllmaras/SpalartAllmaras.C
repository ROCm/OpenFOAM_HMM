/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "SpalartAllmaras.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LES
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmaras, 0);
addToRunTimeSelectionTable(LESmodel, SpalartAllmaras, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmaras::fv1() const
{
    volScalarField chi3 = pow3(nuTilda_/(mu()/rho()));
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> SpalartAllmaras::fv2() const
{
    volScalarField chi = nuTilda_/(mu()/rho());
    //return scalar(1) - chi/(scalar(1) + chi*fv1());
    return 1.0/pow3(scalar(1) + chi/Cv2_);
}


tmp<volScalarField> SpalartAllmaras::fv3() const
{
    volScalarField chi = nuTilda_/(mu()/rho());
    volScalarField chiByCv2 = (1/Cv2_)*chi;

    return
        (scalar(1) + chi*fv1())
       *(1/Cv2_)
       *(3*(scalar(1) + chiByCv2) + sqr(chiByCv2))
       /pow3(scalar(1) + chiByCv2);
}


tmp<volScalarField> SpalartAllmaras::fw(const volScalarField& Stilda) const
{
    volScalarField r = min
    (
        nuTilda_
       /(
           max(Stilda, dimensionedScalar("SMALL", Stilda.dimensions(), SMALL))
          *sqr(kappa_*dTilda_)
        ),
        scalar(10.0)
    );
    r.boundaryField() == 0.0;

    volScalarField g = r + Cw2_*(pow6(r) - r);

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmaras::SpalartAllmaras
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermoPhysicalModel
)
:
    LESmodel(typeName, rho, U, phi, thermoPhysicalModel),

    alphaNut_(LESmodelProperties().lookupOrAddDefault<scalar>("alphaNut", 1.5)),
    Cb1_(LESmodelProperties().lookupOrAddDefault<scalar>("Cb1", 0.1355)),
    Cb2_(LESmodelProperties().lookupOrAddDefault<scalar>("Cb2", 0.622)),
    Cv1_(LESmodelProperties().lookupOrAddDefault<scalar>("Cv1", 7.1)),
    Cv2_(LESmodelProperties().lookupOrAddDefault<scalar>("Cv2", 5.0)),
    CDES_(LESmodelProperties().lookupOrAddDefault<scalar>("CDES", 0.65)),
    ck_(LESmodelProperties().lookupOrAddDefault<scalar>("ck", 0.07)),
    kappa_(lookupOrAddDefault<scalar>("kappa", 0.4187)),
    Cw1_(Cb1_/sqr(kappa_) + alphaNut_*(1.0 + Cb2_)),
    Cw2_(LESmodelProperties().lookupOrAddDefault<scalar>("Cw2", 0.3)),
    Cw3_(LESmodelProperties().lookupOrAddDefault<scalar>("Cw3", 2.0)),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    dTilda_(min(CDES_*delta(), wallDist(mesh_).y())),
    muSgs_
    (
        IOobject
        (
            "muSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )

{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> SpalartAllmaras::B() const
{
    return ((2.0/3.0)*I)*k() - (muSgs_/rho())*dev(twoSymm(fvc::grad(U())));
}


tmp<volSymmTensorField> SpalartAllmaras::devRhoBeff() const
{
    return -muEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<volScalarField> SpalartAllmaras::epsilon() const
{
    return 2*muEff()/rho()*magSqr(symm(fvc::grad(U())));
}


tmp<fvVectorMatrix> SpalartAllmaras::divDevRhoBeff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U) - fvc::div(muEff()*dev2(fvc::grad(U)().T()))
    );
}


void SpalartAllmaras::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();
    LESmodel::correct(gradU);

    if (mesh_.changing())
    {
        dTilda_ = min(CDES_*delta(), wallDist(mesh_).y());
    }

    volScalarField Stilda =
        fv3()*::sqrt(2.0)*mag(skew(gradU)) + fv2()*nuTilda_/sqr(kappa_*dTilda_);

    solve
    (
        fvm::ddt(rho(), nuTilda_)
      + fvm::div(phi(), nuTilda_)
      - fvm::laplacian
        (
            alphaNut_*(nuTilda_*rho() + mu()),
            nuTilda_,
            "laplacian(DnuTildaEff,nuTilda)"
        )
      - alphaNut_*rho()*Cb2_*magSqr(fvc::grad(nuTilda_))
     ==
        rho()*Cb1_*Stilda*nuTilda_
      - fvm::Sp(rho()*Cw1_*fw(Stilda)*nuTilda_/sqr(dTilda_), nuTilda_)
    );

    bound(nuTilda_, dimensionedScalar("zero", nuTilda_.dimensions(), 0.0));

    nuTilda_.correctBoundaryConditions();
    muSgs_.internalField() = rho()*fv1()*nuTilda_.internalField();
    muSgs_.correctBoundaryConditions();
}


bool SpalartAllmaras::read()
{
    if (LESmodel::read())
    {
        LESmodelProperties().readIfPresent<scalar>("alphaNut", alphaNut_);
        LESmodelProperties().readIfPresent<scalar>("Cb1", Cb1_);
        LESmodelProperties().readIfPresent<scalar>("Cb2", Cb2_);
        Cw1_ = Cb1_/sqr(kappa_) + alphaNut_*(1.0 + Cb2_);
        LESmodelProperties().readIfPresent<scalar>("Cw2", Cw2_);
        LESmodelProperties().readIfPresent<scalar>("Cw3", Cw3_);
        LESmodelProperties().readIfPresent<scalar>("Cv1", Cv1_);
        LESmodelProperties().readIfPresent<scalar>("Cv2", Cv2_);
        LESmodelProperties().readIfPresent<scalar>("CDES", CDES_);
        LESmodelProperties().readIfPresent<scalar>("ck", ck_);
        readIfPresent<scalar>("kappa", kappa_);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LES
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
