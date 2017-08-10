/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "MultiComponentPhaseModel.H"

#include "phaseSystem.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
MultiComponentPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseName),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.mesh().solverDict("Yi").lookup("residualAlpha")
    ),
    species_(),
    inertIndex_(-1)
{
    thermoPtr_.set
    (
        phaseThermo::New
        (
            fluid.mesh(),
            phaseName,
            basicThermo::phasePropertyName(basicThermo::dictName, phaseName)
        ).ptr()
    );

    if (thermoPtr_->composition().species().size() == 0)
    {
        FatalErrorInFunction
            << " The selected thermo is pure. Use a multicomponent thermo."
            << exit(FatalError);
    }

    species_ = thermoPtr_->composition().species();

    if (thermoPtr_().found("inertSpecie"))
    {
        inertIndex_ =
            species_
            [
                thermoPtr_().lookup("inertSpecie")
            ];
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
~MultiComponentPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
const phaseThermo&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::thermo() const
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
phaseThermo&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::thermo()
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::YiEqn
(
    volScalarField& Yi
)
{
    if
    (
        (inertIndex_ != -1)
     && (
            Yi.name()
         ==
            IOobject::groupName
            (
                thermoPtr_->composition().species()[inertIndex_],
                this->name()
            )
        )
    )
    {
        return tmp<fvScalarMatrix>();
    }

    const volScalarField& alpha = *this;
    const surfaceScalarField& alphaPhi = this->alphaPhi();

//     surfaceScalarField alphaRhoPhi
//     (
//         fvc::interpolate(alpha)*this->fluid().phi()
//     );

    return
        tmp<fvScalarMatrix>
        (
            fvm::ddt(alpha, Yi)
          + fvm::div(alphaPhi, Yi, "div(" + alphaPhi.name() + ",Yi)")
//          - fvm::Sp(fvc::div(alphaPhi), Yi)
          ==
            fvc::ddt(residualAlpha_, Yi)
          - fvm::ddt(residualAlpha_, Yi)
        );
}


template<class BasePhaseModel, class phaseThermo>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::Y() const
{
    return thermoPtr_->composition().Y();
}


template<class BasePhaseModel, class phaseThermo>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::Y()
{
    return thermoPtr_->composition().Y();
}


template<class BasePhaseModel, class phaseThermo>
Foam::dimensionedScalar
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::hf
(
    label index
) const
{
    return dimensionedScalar
    (
        "Hc",
        sqr(dimLength)/sqr(dimTime),
        thermoPtr_->composition().Hc(index)
    );
}


/*
template<class BasePhaseModel, class phaseThermo>
void Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
correctMassFractions()
{
    volScalarField Yt
    (
        IOobject
        (
            IOobject::groupName("Yt", this->name()),
            this->fluid().mesh().time().timeName(),
            this->fluid().mesh()
        ),
        this->fluid().mesh(),
        dimensionedScalar("zero", dimless, 0)
    );

    PtrList<volScalarField>& Yi = thermoPtr_->composition().Y();

    forAll(Yi, i)
    {
        if (i != inertIndex_)
        {
            Yt += Yi[i];
        }
    }

    if (inertIndex_ != -1)
    {
        Yi[inertIndex_] = scalar(1) - Yt;
        Yi[inertIndex_].max(0);
    }
    else
    {
        forAll(Yi, i)
        {
            Yi[i] /= Yt;
            Yi[i].max(0);
        }
    }
}
*/

template<class BasePhaseModel, class phaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::mu() const
{
    return thermoPtr_->mu();
}


template<class BasePhaseModel, class phaseThermo>
Foam::tmp<Foam::scalarField>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::mu
(
    const label patchI
) const
{
    return thermoPtr_->mu()().boundaryField()[patchI];
}


template<class BasePhaseModel, class phaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::nu() const
{
    return thermoPtr_->nu();
}


template<class BasePhaseModel, class phaseThermo>
Foam::tmp<Foam::scalarField>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::nu
(
    const label patchI
) const
{
    return  thermoPtr_->nu(patchI);
}


// ************************************************************************* //
