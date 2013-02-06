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
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const word& phaseName
)
:
    rhoThermo(mesh, phaseName),
    phaseName1_("1"),
    phaseName2_("2"),
    alpha1_(alpha1),
    alpha2_(alpha2),
    thermo1_(rhoThermo::New(mesh, phaseName1_)),
    thermo2_(rhoThermo::New(mesh, phaseName2_)),
    he_
    (
        IOobject
        (
            phasePropertyName
            (
                "he"
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        heBoundaryTypes(),
        heBoundaryBaseTypes()
    )
{
    thermo1_->validate("phaseModel 1", "e");
    thermo2_->validate("phaseModel 2", "e");

    rho_ = alpha1_*thermo1_->rho() + alpha2_*thermo2_->rho();

    he_ =
    (
        alpha1_*thermo1_->rho()*thermo1_->he()
      + alpha2_*thermo2_->rho()*thermo2_->he()
    )/rho_;

    volScalarField::GeometricBoundaryField& hbf = he_.boundaryField();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }

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

    psi_ = alpha1_*thermo1_->psi() + alpha2_*thermo2_->psi();
    mu_ = alpha1_*thermo1_->mu() + alpha2_*thermo2_->mu();
    alpha_ = alpha1_*thermo1_->alpha() + alpha2_*thermo2_->alpha();
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
    return alpha1_*thermo1_->he(p, T) + alpha2_*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1_, cells)*thermo1_->he(p, T, cells)
      + scalarField(alpha2_, cells)*thermo2_->he(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->he(p, T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->he(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::hc() const
{
    return alpha1_*thermo1_->hc() + alpha2_*thermo2_->hc();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    notImplemented("Foam::twoPhaseThermo::THE");
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
    notImplemented("Foam::twoPhaseThermo::THE");
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::Cp() const
{
    return alpha1_*thermo1_->Cp() + alpha2_*thermo2_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->Cp(p, T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->Cp(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::Cv() const
{
    return alpha1_*thermo1_->Cv() + alpha2_*thermo2_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->Cv(p, T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->Cv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::gamma() const
{
    return alpha1_*thermo1_->gamma() + alpha2_*thermo2_->gamma();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->gamma(p, T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->gamma(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::Cpv() const
{
    return alpha1_*thermo1_->Cpv() + alpha2_*thermo2_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->Cpv(p, T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->Cpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::CpByCpv() const
{
    return alpha1_*thermo1_->CpByCpv() + alpha2_*thermo2_->CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->CpByCpv(p, T, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->CpByCpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::kappa() const
{
    return alpha1_*thermo1_->kappa() + alpha2_*thermo2_->kappa();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::kappa
(
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->kappa(patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1_*thermo1_->kappaEff(alphat)
      + alpha2_*thermo2_->kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->kappaEff(alphat, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->kappaEff(alphat, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1_*thermo1_->alphaEff(alphat)
      + alpha2_*thermo2_->alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1_.boundaryField()[patchi]*thermo1_->alphaEff(alphat, patchi)
      + alpha2_.boundaryField()[patchi]*thermo2_->alphaEff(alphat, patchi);
}


// ************************************************************************* //
