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

#include "fixedReferenceTemperature.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{
    defineTypeNameAndDebug(fixedReferenceTemperature, 0);
    addToRunTimeSelectionTable
    (
        heatTransferCoeffModel,
        fixedReferenceTemperature,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferCoeffModels::fixedReferenceTemperature::fixedReferenceTemperature
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& TName
)
:
    heatTransferCoeffModel(dict, mesh, TName),
    TRef_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::heatTransferCoeffModels::fixedReferenceTemperature::read
(
    const dictionary& dict
)
{
    if (heatTransferCoeffModel::read(dict))
    {
        dict.lookup("TRef") >> TRef_;

        return true;
    }

    return false;
}


void Foam::heatTransferCoeffModels::fixedReferenceTemperature::htc
(
    volScalarField& htc
)
{
    const FieldField<Field, scalar> qBf(q());
    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);
    const volScalarField::Boundary& Tbf = T.boundaryField();
    const scalar eps = ROOTVSMALL;

    volScalarField::Boundary& htcBf = htc.boundaryFieldRef();
    forAllConstIters(patchSet_, iter)
    {
        label patchi = iter.key();
        htcBf[patchi] = qBf[patchi]/(TRef_ - Tbf[patchi] + eps);
    }
}


// ************************************************************************* //
