/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "atmNutkWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> atmNutkWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const auto& turbModel =
        db().lookupObject<turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );
    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    auto tnutw = tmp<scalarField>::New(*this);
    auto& nutw = tnutw.ref();

    const scalar Cmu25 = pow025(Cmu_);

    const scalar t = db().time().timeOutputValue();
    const scalarField z0(z0_->value(t));

    #ifdef FULLDEBUG
    for (const scalar z : z0)
    {
        if (z < VSMALL)
        {
            FatalErrorInFunction
                << "z0 field can only contain positive values. "
                << "Please check input field z0."
                << exit(FatalError);
        }
    }
    #endif

    const labelList& faceCells = patch().faceCells();

    // (HW:Eq. 5)
    forAll(nutw, facei)
    {
        const label celli = faceCells[facei];

        const scalar uStar = Cmu25*sqrt(k[celli]);
        const scalar yPlus = uStar*y[facei]/nuw[facei];
        const scalar Edash = (y[facei] + z0[facei])/z0[facei];

        nutw[facei] = nuw[facei]*(yPlus*kappa_/log(max(Edash, 1 + 1e-4)) - 1);
    }

    if (boundNut_)
    {
        nutw = max(nutw, scalar(0));
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmNutkWallFunctionFvPatchScalarField::atmNutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    boundNut_(true),
    z0_(nullptr)
{}


atmNutkWallFunctionFvPatchScalarField::atmNutkWallFunctionFvPatchScalarField
(
    const atmNutkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    boundNut_(ptf.boundNut_),
    z0_(ptf.z0_.clone(p.patch()))
{}


atmNutkWallFunctionFvPatchScalarField::atmNutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    boundNut_(dict.getOrDefault<bool>("boundNut", false)),
    z0_(PatchFunction1<scalar>::New(p.patch(), "z0", dict))
{}


atmNutkWallFunctionFvPatchScalarField::atmNutkWallFunctionFvPatchScalarField
(
    const atmNutkWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    boundNut_(rwfpsf.boundNut_),
    z0_(rwfpsf.z0_.clone(this->patch().patch()))
{}


atmNutkWallFunctionFvPatchScalarField::atmNutkWallFunctionFvPatchScalarField
(
    const atmNutkWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    boundNut_(rwfpsf.boundNut_),
    z0_(rwfpsf.z0_.clone(this->patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmNutkWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    z0_->autoMap(m);
}


void atmNutkWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const atmNutkWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const atmNutkWallFunctionFvPatchScalarField>(ptf);

    z0_->rmap(nrwfpsf.z0_(), addr);
}


void atmNutkWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    os.writeEntry("boundNut", boundNut_);
    z0_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmNutkWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
