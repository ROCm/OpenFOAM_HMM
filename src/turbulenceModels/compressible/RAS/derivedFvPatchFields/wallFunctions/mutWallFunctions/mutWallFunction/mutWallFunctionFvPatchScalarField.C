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

#include "mutWallFunctionFvPatchScalarField.H"
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

mutWallFunctionFvPatchScalarField::
mutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


mutWallFunctionFvPatchScalarField::
mutWallFunctionFvPatchScalarField
(
    const mutWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


mutWallFunctionFvPatchScalarField::
mutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


mutWallFunctionFvPatchScalarField::
mutWallFunctionFvPatchScalarField
(
    const mutWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


mutWallFunctionFvPatchScalarField::
mutWallFunctionFvPatchScalarField
(
    const mutWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mutWallFunctionFvPatchScalarField::updateCoeffs()
{
    const RASModel& ras = db().lookupObject<RASModel>("RASProperties");

    const scalar Cmu = ras.Cmu().value();
    const scalar Cmu25 = pow(Cmu, 0.25);
    const scalar kappa = ras.kappa().value();
    const scalar E = ras.E().value();
    const scalar yPlusLam = ras.yPlusLam();

    const scalarField& y = ras.y()[patch().index()];

    const scalarField& rhow =
        patch().lookupPatchField<volScalarField, scalar>("rho");

    const volScalarField& k = db().lookupObject<volScalarField>("k");

    const scalarField& muw =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    scalarField& mutw = *this;

    forAll(mutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar yPlus =
            Cmu25*y[faceI]*sqrt(k[faceCellI])
           /(muw[faceI]/rhow[faceI]);

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


void mutWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, mutWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
