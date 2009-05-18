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

#include "mutSpalartAllmarasStandardWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho"),
    muName_("mu")
{}


mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasStandardWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_),
    muName_(ptf.muName_)
{}


mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    muName_(dict.lookupOrDefault<word>("mu", "mu"))
{}


mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasStandardWallFunctionFvPatchScalarField& rwfpsf
)
:
    fixedValueFvPatchScalarField(rwfpsf),
    UName_(rwfpsf.UName_),
    rhoName_(rwfpsf.rhoName_),
    muName_(rwfpsf.muName_)
{}


mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasStandardWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(rwfpsf, iF),
    UName_(rwfpsf.UName_),
    rhoName_(rwfpsf.rhoName_),
    muName_(rwfpsf.muName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const scalar kappa = rasModel.kappa().value();
    const scalar E = rasModel.E().value();
    scalar yPlusLam = rasModel.yPlusLam();

    const scalarField& ry = patch().deltaCoeffs();

    const fvPatchVectorField& U =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    scalarField magUp = mag(U.patchInternalField() - U);

    const scalarField& rhow =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const scalarField& muw =
        patch().lookupPatchField<volScalarField, scalar>(muName_);
    scalarField& mutw = *this;

    scalarField magFaceGradU = mag(U.snGrad());

    forAll(mutw, faceI)
    {
        scalar magUpara = magUp[faceI];

        scalar kappaRe = kappa*magUpara/((muw[faceI]/rhow[faceI])*ry[faceI]);

        scalar yPlus = yPlusLam;
        scalar ryPlusLam = 1.0/yPlus;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yPlus;
            yPlus = (kappaRe + yPlus)/(1.0 + log(E*yPlus));

        } while (mag(ryPlusLam*(yPlus - yPlusLast)) > 0.01 && ++iter < 10 );

        if (yPlus > yPlusLam)
        {
            mutw[faceI] = muw[faceI]*(yPlus*kappa/log(E*yPlus) - 1);
        }
        else
        {
            mutw[faceI] = 0.0;
        }
    }
}


void mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "mu", "mu", muName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
