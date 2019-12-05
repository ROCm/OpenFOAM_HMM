/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "v2WallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cv2_(0.193),
    Bv2_(-0.94)
{}


v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const v2WallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cv2_(ptf.Cv2_),
    Bv2_(ptf.Bv2_)
{}


v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Cv2_(dict.getOrDefault<scalar>("Cv2", 0.193)),
    Bv2_(dict.getOrDefault<scalar>("Bv2", -0.94))
{}


v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const v2WallFunctionFvPatchScalarField& v2wfpsf
)
:
    fixedValueFvPatchField<scalar>(v2wfpsf),
    Cv2_(v2wfpsf.Cv2_),
    Bv2_(v2wfpsf.Bv2_)
{}


v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const v2WallFunctionFvPatchScalarField& v2wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(v2wfpsf, iF),
    Cv2_(v2wfpsf.Cv2_),
    Bv2_(v2wfpsf.Bv2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void v2WallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const scalar Cmu25 = pow025(nutw.Cmu());

    scalarField& v2 = *this;

    // Set v2 wall values
    forAll(v2, facei)
    {
        const label celli = patch().faceCells()[facei];

        const scalar uTau = Cmu25*sqrt(k[celli]);

        const scalar yPlus = uTau*y[facei]/nuw[facei];

        if (nutw.yPlusLam() < yPlus)
        {
            v2[facei] = Cv2_/nutw.kappa()*log(yPlus) + Bv2_;
        }
        else
        {
            v2[facei] = Cv2_*pow4(yPlus);
        }

        v2[facei] *= sqr(uTau);
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


void v2WallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    os.writeEntry("Cv2", Cv2_);
    os.writeEntry("Bv2", Bv2_);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    v2WallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
