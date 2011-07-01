/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "fanPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::fanPressureFvPatchScalarField::fanFlowDirection,
        2
    >::names[] =
    {
        "in",
        "out"
    };
}

const Foam::NamedEnum
<
    Foam::fanPressureFvPatchScalarField::fanFlowDirection,
    2
> Foam::fanPressureFvPatchScalarField::fanFlowDirectionNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    p0_(p.size(), 0.0),
    fanCurve_(),
    direction_(ffdOut)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    p0_(ptf.p0_, mapper),
    fanCurve_(ptf.fanCurve_),
    direction_(ptf.direction_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    p0_("p0", dict, p.size()),
    fanCurve_(dict),
    direction_(fanFlowDirectionNames_.read(dict.lookup("direction")))
{
    // Assign initial pressure by "value"
    fvPatchField<scalar>::operator==(scalarField("value", dict, p.size()));
}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& pfopsf
)
:
    fixedValueFvPatchScalarField(pfopsf),
    phiName_(pfopsf.phiName_),
    rhoName_(pfopsf.rhoName_),
    p0_(pfopsf.p0_),
    fanCurve_(pfopsf.fanCurve_),
    direction_(pfopsf.direction_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& pfopsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pfopsf, iF),
    phiName_(pfopsf.phiName_),
    rhoName_(pfopsf.rhoName_),
    p0_(pfopsf.p0_),
    fanCurve_(pfopsf.fanCurve_),
    direction_(pfopsf.direction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve flux field
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    int dir = 2*direction_ - 1;

    // Average volumetric flow rate
    scalar aveFlowRate = 0;

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        aveFlowRate = dir*gSum(phip)/gSum(patch().magSf());
    }
    else if (phi.dimensions() == dimVelocity*dimArea*dimDensity)
    {
        const scalarField& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        aveFlowRate = dir*gSum(phip/rhop)/gSum(patch().magSf());
    }
    else
    {
        FatalErrorIn("fanPressureFvPatchScalarField::updateCoeffs()")
            << "dimensions of phi are not correct"
                << "\n    on patch " << patch().name()
                << " of field " << dimensionedInternalField().name()
                << " in file " << dimensionedInternalField().objectPath() << nl
                << exit(FatalError);
    }

    // Normal flow through fan
    if (aveFlowRate >= 0.0)
    {
        // Pressure drop for this flow rate
        const scalar pdFan = fanCurve_(aveFlowRate);

        operator==(p0_ - dir*pdFan);
    }
    // Reverse flow
    else
    {
        // Assume that fan has stalled if flow reversed
        // i.e. apply dp for zero flow rate
        const scalar pdFan = fanCurve_(0);

        // Flow speed across patch
        scalarField Up = phip/(patch().magSf());

        // Pressure drop associated withback flow = dynamic pressure
        scalarField pdBackFlow = 0.5*magSqr(Up);

        if (phi.dimensions() == dimVelocity*dimArea*dimDensity)
        {
            const scalarField& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);
            pdBackFlow /= rhop;
        }

        operator==(p0_ - dir*(pdBackFlow + pdFan));
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fanPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    fanCurve_.write(os);
    os.writeKeyword("direction")
        << fanFlowDirectionNames_[direction_] << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
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
