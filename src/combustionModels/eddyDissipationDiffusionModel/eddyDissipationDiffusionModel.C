/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

#include "eddyDissipationDiffusionModel.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationDiffusionModel<CombThermoType, ThermoType>::
eddyDissipationDiffusionModel
(
    const word& modelType,
    const fvMesh& mesh,
    const word& combustionProperties,
    const word& phaseName
)
:
    eddyDissipationModelBase<CombThermoType, ThermoType>
    (
        modelType,
        mesh,
        combustionProperties,
        phaseName
    ),
    Cd_(readScalar(this->coeffs().lookup("Cd")))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationDiffusionModel<CombThermoType, ThermoType>::
~eddyDissipationDiffusionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationDiffusionModel<CombThermoType, ThermoType>::timeScale()
{
    return (max(this->rtTurb(), this->rtDiff()));
}


template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationDiffusionModel<CombThermoType, ThermoType>::rtDiff() const
{
    tmp<volScalarField> tdelta
    (
        new volScalarField
        (
            IOobject
            (
                "tdelta",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("delta", dimLength, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& delta = tdelta.ref();
    delta.ref() = pow(this->mesh().V(), 1.0/3.0);
    delta.correctBoundaryConditions();

    // NOTE: Assume Prt = 1
    return Cd_*this->turbulence().nuEff()/sqr(delta);
}


template<class CombThermoType, class ThermoType>
bool eddyDissipationDiffusionModel<CombThermoType, ThermoType>::read()
{
    if (eddyDissipationModelBase<CombThermoType, ThermoType>::read())
    {
        this->coeffs().lookup("Cd") >> Cd_;
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
