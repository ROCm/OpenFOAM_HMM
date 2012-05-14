/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "temperatureJumpFvPatchScalarField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureJumpFvPatchScalarField::temperatureJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    jumpTable_(0)
{}


Foam::temperatureJumpFvPatchScalarField::temperatureJumpFvPatchScalarField
(
    const temperatureJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper),
    jumpTable_(ptf.jumpTable_().clone().ptr())
{}


Foam::temperatureJumpFvPatchScalarField::temperatureJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    jumpTable_(new DataEntry<scalar>("jumpTable"))
{

    if (this->cyclicPatch().owner())
    {
         jumpTable_ = DataEntry<scalar>::New("jumpTable", dict);
    }

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
}


Foam::temperatureJumpFvPatchScalarField::temperatureJumpFvPatchScalarField
(
    const temperatureJumpFvPatchScalarField& ptf
)
:
    cyclicLduInterfaceField(),
    fixedJumpFvPatchField<scalar>(ptf),
    jumpTable_(ptf.jumpTable_().clone().ptr())
{}


Foam::temperatureJumpFvPatchScalarField::temperatureJumpFvPatchScalarField
(
    const temperatureJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF),
    jumpTable_(ptf.jumpTable_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::temperatureJumpFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    if (this->cyclicPatch().owner())
    {
        jumpTable_->writeData(os);
    }
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       temperatureJumpFvPatchScalarField
   );
}

// ************************************************************************* //
