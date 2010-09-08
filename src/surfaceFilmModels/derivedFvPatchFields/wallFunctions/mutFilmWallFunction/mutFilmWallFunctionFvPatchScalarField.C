/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mutFilmWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "RASModel.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFilmModel.H"
#include "directMappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> mutFilmWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();

    bool ok =
        db().objectRegistry::foundObject
        <surfaceFilmModels::surfaceFilmModel>("surfaceFilmProperties");

    if (!ok)
    {
        // do nothing on construction - film model doesn't exist yet
        return tuTau;
    }

    const label patchI = patch().index();

    // Retrieve phase change mass from surface film model
    const surfaceFilmModels::surfaceFilmModel& filmModel =
        db().objectRegistry::lookupObject
            <surfaceFilmModels::surfaceFilmModel>("surfaceFilmProperties");

    const directMappedWallPolyPatch& wpp =
         refCast<const directMappedWallPolyPatch>(patch().patch());
    const mapDistribute& distMap = wpp.map();
    label filmPatchI = -1;
    forAll(filmModel.primaryPatchIDs(), i)
    {
        if (filmModel.primaryPatchIDs()[i] == patchI)
        {
            filmPatchI = filmModel.filmBottomPatchIDs()[i];
            break;
        }
    }

    scalarField mDotFilm =
        filmModel.massPhaseChangeForPrimary().boundaryField()[filmPatchI];
    distMap.distribute(mDotFilm);


    // Retrieve RAS turbulence model
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patchI];
    const fvPatchScalarField& rhow = rasModel.rho().boundaryField()[patchI];
    const fvPatchScalarField& muw = rasModel.mu().boundaryField()[patchI];
    const tmp<volScalarField> tk = rasModel.k();
    const volScalarField& k = tk();

    const scalar Cmu25 = pow(Cmu_, 0.25);

    forAll(uTau, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar ut = Cmu25*sqrt(k[faceCellI]);

        scalar yPlus = y[faceI]*ut/(muw[faceI]/rhow[faceI]);

        scalar mStar = mDotFilm[faceI]/(y[faceI]*ut);

        scalar factor = 0.0;
        if (yPlus > yPlusCrit_)
        {
            scalar expTerm = exp(min(50.0, B_*mStar));
            scalar powTerm = pow(yPlus, mStar/kappa_);
            factor = mStar/(expTerm*powTerm - 1.0 + ROOTVSMALL);
        }
        else
        {
            scalar expTerm = exp(min(50.0, mStar));
            factor = mStar/(expTerm*yPlus - 1.0 + ROOTVSMALL);
        }

        uTau[faceI] = sqrt(max(0, magGradU[faceI]*ut*factor));
    }

    return tuTau;
}


tmp<scalarField> mutFilmWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField magGradU = mag(Uw.snGrad());
    const scalarField& rhow = rasModel.rho().boundaryField()[patchI];
    const scalarField& muw = rasModel.mu().boundaryField()[patchI];

    return max
    (
        scalar(0),
        rhow*sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - muw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutFilmWallFunctionFvPatchScalarField::mutFilmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(p, iF),
    B_(5.5),
    yPlusCrit_(11.05)
{}


mutFilmWallFunctionFvPatchScalarField::mutFilmWallFunctionFvPatchScalarField
(
    const mutFilmWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    B_(5.5),
    yPlusCrit_(11.05)
{}


mutFilmWallFunctionFvPatchScalarField::mutFilmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutkWallFunctionFvPatchScalarField(p, iF, dict),
    B_(dict.lookupOrDefault("B", 5.5)),
    yPlusCrit_(dict.lookupOrDefault("yPlusCrit", 11.05))
{}


mutFilmWallFunctionFvPatchScalarField::mutFilmWallFunctionFvPatchScalarField
(
    const mutFilmWallFunctionFvPatchScalarField& wfpsf
)
:
    mutkWallFunctionFvPatchScalarField(wfpsf),
    B_(wfpsf.B_),
    yPlusCrit_(wfpsf.yPlusCrit_)
{}


mutFilmWallFunctionFvPatchScalarField::mutFilmWallFunctionFvPatchScalarField
(
    const mutFilmWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(wfpsf, iF),
    B_(wfpsf.B_),
    yPlusCrit_(wfpsf.yPlusCrit_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> mutFilmWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patchI];
    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField& rhow = rasModel.rho().boundaryField()[patchI];
    const scalarField& muw = rasModel.mu().boundaryField()[patchI];

    return y*calcUTau(mag(Uw.snGrad()))/(muw/rhow);
}


void mutFilmWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    os.writeKeyword("B") << B_ << token::END_STATEMENT << nl;
    os.writeKeyword("yPlusCrit") << yPlusCrit_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, mutFilmWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
