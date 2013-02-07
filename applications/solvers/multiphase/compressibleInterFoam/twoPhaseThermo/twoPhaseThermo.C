/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "twoPhaseThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseThermo, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseThermo::twoPhaseThermo
(
    const twoPhaseMixture& twoPhaseProperties
)
:
    rhoThermo(twoPhaseProperties.alpha1().mesh(), word::null),
    tpm_(twoPhaseProperties),
    thermo1_(rhoThermo::New(tpm_.alpha1().mesh(), tpm_.phase1Name())),
    thermo2_(rhoThermo::New(tpm_.alpha1().mesh(), tpm_.phase2Name()))
{
    thermo1_->validate(tpm_.phase1Name(), "e");
    thermo2_->validate(tpm_.phase2Name(), "e");

    rho_ = tpm_.alpha1()*thermo1_->rho() + tpm_.alpha2()*thermo2_->rho();

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseThermo::~twoPhaseThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::twoPhaseThermo::correct()
{
    thermo1_->he() = thermo1_->he(p_, T_);
    thermo1_->correct();

    thermo2_->he() = thermo2_->he(p_, T_);
    thermo2_->correct();

    psi_ = tpm_.alpha1()*thermo1_->psi() + tpm_.alpha2()*thermo2_->psi();
    mu_ = tpm_.alpha1()*thermo1_->mu() + tpm_.alpha2()*thermo2_->mu();
    alpha_ = tpm_.alpha1()*thermo1_->alpha() + tpm_.alpha2()*thermo2_->alpha();
}


bool Foam::twoPhaseThermo::incompressible() const
{
    return thermo1_->incompressible() && thermo2_->incompressible();
}


bool Foam::twoPhaseThermo::isochoric() const
{
    return thermo1_->isochoric() && thermo2_->isochoric();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return tpm_.alpha1()*thermo1_->he(p, T) + tpm_.alpha2()*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(tpm_.alpha1(), cells)*thermo1_->he(p, T, cells)
      + scalarField(tpm_.alpha2(), cells)*thermo2_->he(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->he(p, T, patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->he(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::hc() const
{
    return tpm_.alpha1()*thermo1_->hc() + tpm_.alpha2()*thermo2_->hc();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    notImplemented("twoPhaseThermo::THE(...)");
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    notImplemented("twoPhaseThermo::THE(...)");
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::Cp() const
{
    return tpm_.alpha1()*thermo1_->Cp() + tpm_.alpha2()*thermo2_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->Cp(p, T, patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->Cp(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::Cv() const
{
    return tpm_.alpha1()*thermo1_->Cv() + tpm_.alpha2()*thermo2_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->Cv(p, T, patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->Cv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::gamma() const
{
    return tpm_.alpha1()*thermo1_->gamma() + tpm_.alpha2()*thermo2_->gamma();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->gamma(p, T, patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->gamma(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::Cpv() const
{
    return tpm_.alpha1()*thermo1_->Cpv() + tpm_.alpha2()*thermo2_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->Cpv(p, T, patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->Cpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::CpByCpv() const
{
    return
        tpm_.alpha1()*thermo1_->CpByCpv()
      + tpm_.alpha2()*thermo2_->CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->CpByCpv(p, T, patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->CpByCpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::kappa() const
{
    return tpm_.alpha1()*thermo1_->kappa() + tpm_.alpha2()*thermo2_->kappa();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::kappa
(
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->kappa(patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        tpm_.alpha1()*thermo1_->kappaEff(alphat)
      + tpm_.alpha2()*thermo2_->kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->kappaEff(alphat, patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->kappaEff(alphat, patchi)
    ;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        tpm_.alpha1()*thermo1_->alphaEff(alphat)
      + tpm_.alpha2()*thermo2_->alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        tpm_.alpha1().boundaryField()[patchi]*thermo1_->alphaEff(alphat, patchi)
      + tpm_.alpha2().boundaryField()[patchi]*thermo2_->alphaEff(alphat, patchi)
    ;
}


// ************************************************************************* //
