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

#include "ClausiusClapeyron.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationPressureModels
{
    defineTypeNameAndDebug(ClausiusClapeyron, 0);
    addToRunTimeSelectionTable
    (
        saturationPressureModel,
        ClausiusClapeyron,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationPressureModels::ClausiusClapeyron::ClausiusClapeyron
(
    const dictionary& dict
)
:
    saturationPressureModel(),
    p0_("p0", dimPressure, dict.lookup("p0")),
    Tb_("Tb", dimTemperature, dict.lookup("Tb")),
    dHvm_
    (
        "dHvm",
        dimMass*sqr(dimLength)/sqr(dimTime)/dimMoles ,
        dict.lookup("dHvm")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationPressureModels::ClausiusClapeyron::~ClausiusClapeyron()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::saturationPressureModels::ClausiusClapeyron::pSat
(
    const volScalarField& T
) const
{
    volScalarField invTc(1.0/Tb_ - 1.0/T);

    return (p0_*exp(dHvm_*invTc/constant::physicoChemical::R));
}


Foam::tmp<Foam::volScalarField>
Foam::saturationPressureModels::ClausiusClapeyron::pSatPrime
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pSatPrime",
                T.mesh().time().timeName(),
                T.mesh()
            ),
            T.mesh(),
            dimensionedScalar("0", dimless, 0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::saturationPressureModels::ClausiusClapeyron::lnPSat
(
    const volScalarField& T
) const
{
    volScalarField invTc(1.0/Tb_ - 1.0/T);

    return
        log(p0_.value()) + dHvm_*invTc/constant::physicoChemical::R.value();
}


// ************************************************************************* //
