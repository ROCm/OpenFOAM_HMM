/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "greyDiffusiveViewFactorFixedValueFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Qro_(),
    solarLoad_(false)
{}


Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Qro_(ptf.Qro_, mapper),
    solarLoad_(ptf.solarLoad_)
{}


Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    Qro_("Qro", dict, p.size()),
    solarLoad_(dict.lookupOrDefault<bool>("solarLoad", false))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );

    }
    else
    {
         fvPatchScalarField::operator=(0.0);
    }
}


Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    Qro_(ptf.Qro_),
    solarLoad_(ptf.solarLoad_)
{}


Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    Qro_(ptf.Qro_),
    solarLoad_(ptf.solarLoad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    if (debug)
    {
        scalar Q = gSum((*this)*patch().magSf());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << " heat transfer rate:" << Q
            << " wall radiative heat flux "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::radiation::
greyDiffusiveViewFactorFixedValueFvPatchScalarField::Qro() const
{
    tmp<scalarField> tQrt(new scalarField(Qro_));

    if (solarLoad_)
    {
        const radiationModel& radiation =
            db().lookupObject<radiationModel>("radiationProperties");

        tQrt.ref() += patch().lookupPatchField<volScalarField,scalar>
        (
            radiation.externalRadHeatFieldName_
        );
    }

    return tQrt;
}


void Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    Qro_.writeEntry("Qro", os);
    os.writeKeyword("solarLoad") << solarLoad_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
    );
}
}


// ************************************************************************* //
