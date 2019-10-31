/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "interfaceCompositionModel.H"
#include "phaseModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceCompositionModel, 0);
    defineRunTimeSelectionTable(interfaceCompositionModel, dictionary);
}


const Foam::Enum<Foam::interfaceCompositionModel::modelVariable>
Foam::interfaceCompositionModel::modelVariableNames
{
    { modelVariable::T, "temperature" },
    { modelVariable::P, "pressure" },
    { modelVariable::Y, "massFraction" },
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::interfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    modelVariable_
    (
        modelVariableNames.lookupOrDefault
        (
            "variable",
            dict,
            modelVariable::T
        )
    ),
    pair_(pair),
    speciesName_(dict.lookupOrDefault<word>("species", "none")),
    mesh_(pair_.from().mesh())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::interfaceCompositionModel::transferSpecie() const
{
    return speciesName_;
}


const Foam::phasePair& Foam::interfaceCompositionModel::pair() const
{
    return pair_;
}


const Foam::word Foam::interfaceCompositionModel::variable() const
{
    return modelVariableNames[modelVariable_];
}


// ************************************************************************* //
