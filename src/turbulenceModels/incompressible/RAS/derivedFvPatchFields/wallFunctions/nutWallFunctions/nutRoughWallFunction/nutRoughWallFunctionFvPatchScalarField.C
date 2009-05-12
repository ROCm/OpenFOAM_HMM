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

#include "nutRoughWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar nutRoughWallFunctionFvPatchScalarField::fnRough
(
    const scalar KsPlus,
    const scalar Cs,
    const scalar kappa
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

nutRoughWallFunctionFvPatchScalarField::
nutRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    kName_("k"),
    nuName_("nu"),
    Ks_(p.size(), 0.0),
    Cs_(p.size(), 0.0)
{}


nutRoughWallFunctionFvPatchScalarField::
nutRoughWallFunctionFvPatchScalarField
(
    const nutRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    kName_(ptf.kName_),
    nuName_(ptf.nuName_),
    Ks_(ptf.Ks_, mapper),
    Cs_(ptf.Cs_, mapper)
{}


nutRoughWallFunctionFvPatchScalarField::
nutRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size())
{}


nutRoughWallFunctionFvPatchScalarField::
nutRoughWallFunctionFvPatchScalarField
(
    const nutRoughWallFunctionFvPatchScalarField& nrwfpsf
)
:
    fixedValueFvPatchScalarField(nrwfpsf),
    kName_(nrwfpsf.kName_),
    nuName_(nrwfpsf.nuName_),
    Ks_(nrwfpsf.Ks_),
    Cs_(nrwfpsf.Cs_)
{}


nutRoughWallFunctionFvPatchScalarField::
nutRoughWallFunctionFvPatchScalarField
(
    const nutRoughWallFunctionFvPatchScalarField& nrwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(nrwfpsf, iF),
    kName_(nrwfpsf.kName_),
    nuName_(nrwfpsf.nuName_),
    Ks_(nrwfpsf.Ks_),
    Cs_(nrwfpsf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    Ks_.autoMap(m);
    Cs_.autoMap(m);
}


void nutRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const nutRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutRoughWallFunctionFvPatchScalarField>(ptf);

    Cs_.rmap(nrwfpsf.Cs_, addr);
    Ks_.rmap(nrwfpsf.Ks_, addr);
}


void nutRoughWallFunctionFvPatchScalarField::updateCoeffs()
{
    const RASModel& ras = db().lookupObject<RASModel>("RASProperties");

    const scalar Cmu = ras.Cmu().value();
    const scalar Cmu25 = pow(Cmu, 0.25);
    const scalar kappa = ras.kappa().value();
    const scalar E = ras.E().value();
    const scalar yPlusLam = ras.yPlusLam();

    const scalarField& y = ras.y()[patch().index()];

    const scalarField& k = db().lookupObject<volScalarField>(kName_);

    const scalarField& nuw =
        patch().lookupPatchField<volScalarField, scalar>(nuName_);

    scalarField& nutw = *this;

    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uStar = Cmu25*sqrt(k[faceCellI]);
        scalar yPlus = uStar*y[faceI]/nuw[faceI];
        scalar KsPlus = uStar*Ks_[faceI]/nuw[faceI];

        scalar Edash = E;

        if (KsPlus > 2.25)
        {
            Edash = E/fnRough(KsPlus, Cs_[faceI], kappa);
        }

        if (yPlus > yPlusLam)
        {
            scalar limitingNutw = max(nutw[faceI], nuw[faceI]);

            // To avoid oscillations limit the change in the wall viscosity
            // which is particularly important if it temporarily becomes zero
            nutw[faceI] =
                max
                (
                    min
                    (
                        nuw[faceI]*(yPlus*kappa/log(Edash*yPlus) - 1),
                        2*limitingNutw
                    ), 0.5*limitingNutw
                );
        }
        else
        {
            nutw[faceI] = 0.0;
        }

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", nutw = " << nutw[faceI]
                << endl;
        }
    }
}


void nutRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    Cs_.writeEntry("Cs", os);
    Ks_.writeEntry("Ks", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutRoughWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
