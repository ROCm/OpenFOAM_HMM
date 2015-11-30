/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "kOmegaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSSTBase<DESModel<BasicTurbulenceModel> >::correctNut(S2);

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;

    return min(CDES*this->delta(), sqrt(k)/(this->betaStar_*omega));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTDES<BasicTurbulenceModel>::kOmegaSSTDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTBase<DESModel<BasicTurbulenceModel> >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CDESkom_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDESkom",
            this->coeffDict_,
            0.78
        )
    ),
    CDESkeps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDESkeps",
            this->coeffDict_,
            0.61
        )
    )
{
    correctNut();

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTBase<DESModel<BasicTurbulenceModel> >::read())
    {
        kappa_.readIfPresent(this->coeffDict());
        CDESkom_.readIfPresent(this->coeffDict());
        CDESkeps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& k = this->k_;
    volScalarField& omega = this->omega_;
    volScalarField& nut = this->nut_;

    DESModel<BasicTurbulenceModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField magGradU(mag(tgradU()));
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField GbyNu((tgradU() && dev(twoSymm(tgradU()))));
    volScalarField G(this->GName(), nut*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega.boundaryField().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
    );

    volScalarField F1(this->F1(CDkOmega));

    {
        volScalarField gamma(this->gamma(F1));
        volScalarField beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega)
          + fvm::div(alphaRhoPhi, omega)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), omega)
         ==
            alpha*rho*gamma*GbyNu // Using unlimited GybNu
          - fvm::SuSp((2.0/3.0)*alpha*rho*gamma*divU, omega)
          - fvm::Sp(alpha*rho*beta*omega, omega)
          - fvm::SuSp(alpha*rho*(F1 - scalar(1))*CDkOmega/omega, omega)
          + this->omegaSource()
        );

        omegaEqn().relax();

        omegaEqn().boundaryManipulate(omega.boundaryField());

        solve(omegaEqn);
        bound(omega, this->omegaMin_);
    }

    {
        volScalarField CDES(this->CDES(F1));
        volScalarField dTilda(this->dTilda(magGradU, CDES));

        // Turbulent kinetic energy equation
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, k)
          + fvm::div(alphaRhoPhi, k)
          - fvm::laplacian(alpha*rho*this->DkEff(F1), k)
         ==
            min(alpha*rho*G, (this->c1_*this->betaStar_)*alpha*rho*k*omega)
          - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k)
          - fvm::Sp(alpha*rho*sqrt(k)/dTilda, k)  // modified for DES
          + this->kSource()
        );

        kEqn().relax();
        solve(kEqn);
        bound(k, this->kMin_);
    }

    this->correctNut(S2);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicTurbulenceModel>::LESRegion() const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volVectorField& U = this->U_;

    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
    );

    const volScalarField F1(this->F1(CDkOmega));

    tmp<volScalarField> tLESRegion
    (
        new volScalarField
        (
            IOobject
            (
                "DES::LESRegion",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neg
            (
                dTilda
                (
                    mag(fvc::grad(U)),
                    F1*CDESkom_ + (1 - F1)*CDESkeps_
                )
              - sqrt(k)/(this->betaStar_*omega)
            )
        )
    );

    return tLESRegion;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
