/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "fanFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Tuple2.H"
#include "PolynomialEntry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::fanFvPatchField<Foam::scalar>::calcFanJump()
{
    if (this->cyclicPatch().owner())
    {
        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);

        const fvsPatchField<scalar>& phip =
            patch().patchField<surfaceScalarField, scalar>(phi);

        scalarField Un(max(phip/patch().magSf(), scalar(0)));
        if (uniformJump_)
        {
            scalar area = gSum(patch().magSf());
            Un = gSum(Un*patch().magSf())/area;

            if (nonDimensional_)
            {
                // Create an adimensional volumetric flow rate
                Un =
                    120.0*Un/pow3(constant::mathematical::pi)
                  * patch().magSf()
                  / pow3(dm_)/rpm_;
            }
        }

        if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            Un /= patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        }

        if (nonDimensional_)
        {
            scalarField deltap(this->jumpTable_->value(Un));

            // Convert adimensional deltap from curve into deltaP
            scalarField pdFan
            (
                deltap*pow4(constant::mathematical::pi)*sqr(dm_*rpm_)/1800.0
            );

            this->jump_ = max(pdFan, scalar(0));
        }
        else
        {
            this->jump_ = max(this->jumpTable_->value(Un), scalar(0));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<>
Foam::fanFvPatchField<Foam::scalar>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    uniformJumpFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    uniformJump_(dict.lookupOrDefault<bool>("uniformJump", false)),
    nonDimensional_(dict.lookupOrDefault<Switch>("nonDimensional", false)),
    rpm_(dict.lookupOrDefault<scalar>("rpm", 0.0)),
    dm_(dict.lookupOrDefault<scalar>("dm", 0.0))
{
    if (this->cyclicPatch().owner())
    {
        this->jumpTable_ = Function1<scalar>::New("jumpTable", dict);
    }

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTemplatePatchTypeField(scalar, fan);
}


// ************************************************************************* //
