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

#include "dissipationSelector.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* NamedEnum<dissipationSelector::treatment, 3>::names[] =
{
    "none",
    "cascade",
    "equilibrium"
};

const NamedEnum<dissipationSelector::treatment, 3>
    dissipationSelector::typeNames;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


dissipationSelector::dissipationSelector
(
    const turbulenceModel& turbModel,
    const treatment treatmentType
)
:
    turbModel_(turbModel),
    treatment_(treatmentType)
{}


dissipationSelector::dissipationSelector
(
    const turbulenceModel& turbModel,
    const dictionary& dict
)
:
    turbModel_(turbModel),
    treatment_(none)
{

    word modelType;

    if (dict.readIfPresent<word>("dissipation", modelType))
    {
        treatment_ = typeNames[modelType];
    }
}


tmp<volScalarField>
dissipationSelector::dissipation() const
{
    if (treatment_ == cascade)
    {
        return turbModel_.thermalDissipation();
    }
    else if (treatment_ == equilibrium)
    {
        return turbModel_.thermalDissipationEff();
    }
    else
    {
        // a bit wasteful, but we'll avoid it with 'enabled' query anyhow
        tmp<volScalarField> tField = turbModel_.thermalDissipation();
        tField() = 0.0;

        return tField;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
