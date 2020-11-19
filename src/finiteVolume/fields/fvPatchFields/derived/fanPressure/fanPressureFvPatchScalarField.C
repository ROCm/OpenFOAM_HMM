/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd
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

#include "fanPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "TableFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fanPressureFvPatchScalarField::fanFlowDirection
>
Foam::fanPressureFvPatchScalarField::fanFlowDirectionNames_
({
    { fanFlowDirection::ffdIn, "in" },
    { fanFlowDirection::ffdOut, "out" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(p, iF),
    fanCurve_(nullptr),
    direction_(ffdOut),
    nonDimensional_(false),
    rpm_(0),
    dm_(0)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(ptf, p, iF, mapper),
    fanCurve_(ptf.fanCurve_.clone()),
    direction_(ptf.direction_),
    nonDimensional_(ptf.nonDimensional_),
    rpm_(ptf.rpm_),
    dm_(ptf.dm_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF, dict),
    fanCurve_(nullptr),
    direction_(fanFlowDirectionNames_.get("direction", dict)),
    nonDimensional_(dict.getOrDefault("nonDimensional", false)),
    rpm_(0),
    dm_(0)
{
    // Backwards compatibility
    if (dict.found("file"))
    {
        fanCurve_.reset
        (
            new Function1Types::TableFile<scalar>("fanCurve", dict)
        );
    }
    else
    {
        fanCurve_.reset
        (
            Function1<scalar>::New("fanCurve", dict)
        );
    }

    if (nonDimensional_)
    {
        dict.readEntry("rpm", rpm_);
        dict.readEntry("dm", dm_);
    }
}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& fppsf
)
:
    totalPressureFvPatchScalarField(fppsf),
    fanCurve_(fppsf.fanCurve_.clone()),
    direction_(fppsf.direction_),
    nonDimensional_(fppsf.nonDimensional_),
    rpm_(fppsf.rpm_),
    dm_(fppsf.dm_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& fppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(fppsf, iF),
    fanCurve_(fppsf.fanCurve_.clone()),
    direction_(fppsf.direction_),
    nonDimensional_(fppsf.nonDimensional_),
    rpm_(fppsf.rpm_),
    dm_(fppsf.dm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve flux field
    const auto& phi = db().lookupObject<surfaceScalarField>(phiName());

    const auto& phip = patch().patchField<surfaceScalarField, scalar>(phi);

    const int dir = 2*direction_ - 1;

    // Average volumetric flow rate
    scalar volFlowRate = 0;

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        volFlowRate = dir*gSum(phip);
    }
    else if (phi.dimensions() == dimVelocity*dimArea*dimDensity)
    {
        const scalarField& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName());
        volFlowRate = dir*gSum(phip/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct\n"
            << "    on patch " << patch().name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << nl
            << exit(FatalError);
    }

    if (nonDimensional_)
    {
        // Create an non-dimensional flow rate
        volFlowRate =
            120.0*volFlowRate/pow3(constant::mathematical::pi)/pow3(dm_)/rpm_;
    }

    // Pressure drop for this flow rate
    scalar pdFan = fanCurve_->value(max(volFlowRate, scalar(0)));

    if (nonDimensional_)
    {
        // Convert the non-dimensional deltap from curve into deltaP
        pdFan = pdFan*pow4(constant::mathematical::pi)*sqr(dm_*rpm_)/1800;
    }

    totalPressureFvPatchScalarField::updateCoeffs
    (
        p0() - dir*pdFan,
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}


void Foam::fanPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);
    fanCurve_->writeData(os);
    os.writeEntry("direction", fanFlowDirectionNames_[direction_]);

    if (nonDimensional_)
    {
        os.writeEntry("nonDimensional", "true");
        os.writeEntry("rpm", rpm_);
        os.writeEntry("dm", dm_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fanPressureFvPatchScalarField
    );
};


// ************************************************************************* //
