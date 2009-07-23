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
    mutWallFunctionFvPatchScalarField(p, iF)
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
    mutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasStandardWallFunctionFvPatchScalarField& sawfpsf
)
:
    mutWallFunctionFvPatchScalarField(sawfpsf)
{}


mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasStandardWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(sawfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalar yPlusLam = rasModel.yPlusLam(kappa_, E_);
    const scalarField& y = rasModel.y()[patchI];

    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];

    const scalarField magUp = mag(Uw.patchInternalField() - Uw);

    const fvPatchScalarField& rhow = rasModel.rho().boundaryField()[patchI];

    const fvPatchScalarField& muw = rasModel.mu().boundaryField()[patchI];

    tmp<scalarField> tmutw(new scalarField(patch().size(), 0.0));
    scalarField& mutw = tmutw();

    forAll(mutw, faceI)
    {
        scalar magUpara = magUp[faceI];

        scalar kappaRe = kappa_*magUpara*y[faceI]/(muw[faceI]/rhow[faceI]);

        scalar yPlus = yPlusLam;
        scalar ryPlusLam = 1.0/yPlus;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yPlus;
            yPlus = (kappaRe + yPlus)/(1.0 + log(E_*yPlus));

        } while (mag(ryPlusLam*(yPlus - yPlusLast)) > 0.01 && ++iter < 10);

        if (yPlus > yPlusLam)
        {
            mutw[faceI] = muw[faceI]*(yPlus*kappa_/log(E_*yPlus) - 1.0);
        }
    }

    return tmutw;
}


tmp<scalarField>
mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::yPlus() const
{
    notImplemented
    (
        "mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::yPlus() "
        "const"
    );

    return tmp<scalarField>(NULL);
}


void mutSpalartAllmarasStandardWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
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
