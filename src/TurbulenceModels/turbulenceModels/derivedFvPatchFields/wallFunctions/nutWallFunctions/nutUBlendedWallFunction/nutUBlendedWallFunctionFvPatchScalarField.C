/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "nutUBlendedWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUBlendedWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}


Foam::tmp<Foam::scalarField>
Foam::nutUBlendedWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const vectorField n(patch().nf());
    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    vectorField Up(Uw.patchInternalField() - Uw);
    Up -= n*(n & Up);
    const scalarField magUp(mag(Up));

    tmp<scalarField> tuTaup(new scalarField(patch().size(), Zero));
    scalarField& uTaup = tuTaup.ref();

    const scalarField& nutw = *this;

    forAll(uTaup, facei)
    {
        scalar ut = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);
        if (mag(ut) > ROOTVSMALL)
        {
            scalar error = GREAT;
            label iter = 0;
            while (iter++ < 10 && error > 0.001)
            {
                const scalar yPlus = y[facei]*ut/nuw[facei];
                const scalar uTauVis = magUp[facei]/yPlus;
                const scalar uTauLog = kappa_*magUp[facei]/log(E_*yPlus);

                const scalar utNew =
                    pow(pow(uTauVis, n_) + pow(uTauLog, n_), 1.0/n_);
                error = mag(ut - utNew)/(ut + ROOTVSMALL);
                ut = 0.5*(ut + utNew);
            }
        }
        uTaup[facei] = ut;
    }

    return tuTaup;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutUBlendedWallFunctionFvPatchScalarField::
nutUBlendedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    n_(4)
{}


Foam::nutUBlendedWallFunctionFvPatchScalarField::
nutUBlendedWallFunctionFvPatchScalarField
(
    const nutUBlendedWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    n_(ptf.n_)
{}


Foam::nutUBlendedWallFunctionFvPatchScalarField::
nutUBlendedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    n_(dict.getOrDefault<scalar>("n", 4.0))
{}


Foam::nutUBlendedWallFunctionFvPatchScalarField::
nutUBlendedWallFunctionFvPatchScalarField
(
    const nutUBlendedWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    n_(wfpsf.n_)
{}


Foam::nutUBlendedWallFunctionFvPatchScalarField::
nutUBlendedWallFunctionFvPatchScalarField
(
    const nutUBlendedWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    n_(wfpsf.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUBlendedWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));

    return y*calcUTau(magGradU)/nuw;
}


void Foam::nutUBlendedWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    os.writeEntry("n", n_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutUBlendedWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
