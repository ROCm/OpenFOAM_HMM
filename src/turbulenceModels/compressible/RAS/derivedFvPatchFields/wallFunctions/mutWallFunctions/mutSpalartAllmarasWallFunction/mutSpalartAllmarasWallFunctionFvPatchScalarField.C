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

#include "mutSpalartAllmarasWallFunctionFvPatchScalarField.H"
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

mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
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


mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasWallFunctionFvPatchScalarField& ptf,
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


mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
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


mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasWallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    UName_(wfpsf.UName_),
    rhoName_(wfpsf.rhoName_),
    muName_(wfpsf.muName_)
{}


mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    UName_(wfpsf.UName_),
    rhoName_(wfpsf.rhoName_),
    muName_(wfpsf.muName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mutSpalartAllmarasWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const scalar kappa = rasModel.kappa().value();
    const scalar E = rasModel.E().value();

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

        scalar utau =
            sqrt((mutw[faceI] + muw[faceI])*magFaceGradU[faceI]/rhow[faceI]);

        if (utau > VSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa*magUpara/utau, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - utau/(ry[faceI]*(muw[faceI]/rhow[faceI]))
                    + magUpara/utau
                    + 1/E*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    1.0/(ry[faceI]*(muw[faceI]/rhow[faceI]))
                  + magUpara/sqr(utau)
                  + 1/E*kUu*fkUu/utau;

                scalar utauNew = utau + f/df;
                err = mag((utau - utauNew)/utau);
                utau = utauNew;

            } while (utau > VSMALL && err > 0.01 && ++iter < 10);

            mutw[faceI] = max
            (
                rhow[faceI]*sqr(max(utau, 0))/magFaceGradU[faceI]- muw[faceI],
                0.0
            );
        }
        else
        {
            mutw[faceI] = 0;
        }
    }
}


void mutSpalartAllmarasWallFunctionFvPatchScalarField::write
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

makePatchTypeField(fvPatchScalarField, mutSpalartAllmarasWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
