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

#include "mutRoughWallFunctionFvPatchScalarField.H"
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar mutRoughWallFunctionFvPatchScalarField::fnRough
(
    const scalar KsPlus,
    const scalar Cs
) const
{
    // Return fn based on non-dimensional roughness height

    if (KsPlus < 90.0)
    {
        return pow
        (
            (KsPlus - 2.25)/87.75 + Cs*KsPlus,
            sin(0.4258*(log(KsPlus) - 0.811))
        );
    }
    else
    {
        return (1.0 + Cs*KsPlus);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutRoughWallFunctionFvPatchScalarField::
mutRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    rhoName_("rho"),
    muName_("mu"),
    kName_("k"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    Ks_(p.size(), 0.0),
    Cs_(p.size(), 0.0)
{}


mutRoughWallFunctionFvPatchScalarField::
mutRoughWallFunctionFvPatchScalarField
(
    const mutRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_),
    muName_(ptf.muName_),
    kName_(ptf.kName_),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    Ks_(ptf.Ks_, mapper),
    Cs_(ptf.Cs_, mapper)
{}


mutRoughWallFunctionFvPatchScalarField::
mutRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    muName_(dict.lookupOrDefault<word>("mu", "mu")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size())
{}


mutRoughWallFunctionFvPatchScalarField::
mutRoughWallFunctionFvPatchScalarField
(
    const mutRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    fixedValueFvPatchScalarField(rwfpsf),
    rhoName_(rwfpsf.rhoName_),
    muName_(rwfpsf.muName_),
    kName_(rwfpsf.kName_),
    Cmu_(rwfpsf.Cmu_),
    kappa_(rwfpsf.kappa_),
    E_(rwfpsf.E_),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


mutRoughWallFunctionFvPatchScalarField::
mutRoughWallFunctionFvPatchScalarField
(
    const mutRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(rwfpsf, iF),
    rhoName_(rwfpsf.rhoName_),
    muName_(rwfpsf.muName_),
    kName_(rwfpsf.kName_),
    Cmu_(rwfpsf.Cmu_),
    kappa_(rwfpsf.kappa_),
    E_(rwfpsf.E_),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mutRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    Ks_.autoMap(m);
    Cs_.autoMap(m);
}


void mutRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const mutRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const mutRoughWallFunctionFvPatchScalarField>(ptf);

    Cs_.rmap(nrwfpsf.Cs_, addr);
    Ks_.rmap(nrwfpsf.Ks_, addr);
}


void mutRoughWallFunctionFvPatchScalarField::updateCoeffs()
{
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    scalar yPlusLam = rasModel.yPlusLam(kappa_, E_);
    const scalarField& y = rasModel.y()[patch().index()];

    const scalar Cmu25 = pow(Cmu_, 0.25);

    const scalarField& rhow =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const scalarField& k = db().lookupObject<volScalarField>(kName_);

    const scalarField& muw =
        patch().lookupPatchField<volScalarField, scalar>(muName_);

    scalarField& mutw = *this;

    forAll(mutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uStar = Cmu25*sqrt(k[faceCellI]);

        scalar yPlus = uStar*y[faceI]/(muw[faceI]/rhow[faceI]);

        scalar KsPlus = uStar*Ks_[faceI]/(muw[faceI]/rhow[faceI]);

        scalar Edash = E_;
        scalar yPlusLamNew = yPlusLam;
        if (KsPlus > 2.25)
        {
            Edash /= fnRough(KsPlus, Cs_[faceI]);
            yPlusLamNew = rasModel.yPlusLam(kappa_, Edash);
        }

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", yPlusLam = " << yPlusLam
                << endl;
        }

        if (yPlus > yPlusLamNew)
        {
            mutw[faceI] =
                muw[faceI]*(yPlus*kappa_/log(max(Edash*yPlus, 1+1e-4)) - 1);
        }
        else
        {
            mutw[faceI] = 0.0;
        }
    }
}


void mutRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "mu", "mu", muName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    Cs_.writeEntry("Cs", os);
    Ks_.writeEntry("Ks", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, mutRoughWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
