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
#include "enthalpyJumpFvPatchScalarField.H"
#include "temperatureJumpFvPatchScalarField.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enthalpyJumpFvPatchScalarField::enthalpyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{}


Foam::enthalpyJumpFvPatchScalarField::enthalpyJumpFvPatchScalarField
(
    const enthalpyJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::enthalpyJumpFvPatchScalarField::enthalpyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
}


Foam::enthalpyJumpFvPatchScalarField::enthalpyJumpFvPatchScalarField
(
    const enthalpyJumpFvPatchScalarField& ptf
)
:
    cyclicLduInterfaceField(),
    fixedJumpFvPatchField<scalar>(ptf)
{}


Foam::enthalpyJumpFvPatchScalarField::enthalpyJumpFvPatchScalarField
(
    const enthalpyJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::enthalpyJumpFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicPatch().owner())
    {
        const basicThermo& thermo = db().lookupObject<basicThermo>
        (
            "thermophysicalProperties"
        );

        label patchID = patch().index();

        const temperatureJumpFvPatchScalarField& TbPatch =
            refCast<const temperatureJumpFvPatchScalarField>
            (
                thermo.T().boundaryField()[patchID]
            );

        const scalar time = this->db().time().value();
        const scalarField jumpTb
        (
            patch().size(), TbPatch.jumpTable().value(time)
        );

        const labelUList& faceCells = this->patch().faceCells();

        if (db().foundObject<volScalarField>("h"))
        {
            jump_ = thermo.h(jumpTb, faceCells)();
        }
        else if (db().foundObject<volScalarField>("hs"))
        {
            jump_ = thermo.hs(jumpTb, faceCells)();
        }
        else
        {
             FatalErrorIn("enthalpyJumpFvPatchScalarField::updateCoeffs()")
            << " hs or h are not found in db()"
            << exit(FatalError);
        }
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::enthalpyJumpFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       enthalpyJumpFvPatchScalarField
   );
}

// ************************************************************************* //
