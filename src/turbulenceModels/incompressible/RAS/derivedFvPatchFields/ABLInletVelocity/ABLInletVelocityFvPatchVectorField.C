/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "ABLInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ABLInletVelocityFvPatchVectorField::ABLInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Ustar_(0),
    n_(pTraits<vector>::zero),
    z_(pTraits<vector>::zero),
    z0_(0),
    kappa_(0.41)
{}


ABLInletVelocityFvPatchVectorField::ABLInletVelocityFvPatchVectorField
(
    const ABLInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Ustar_(ptf.Ustar_),
    n_(ptf.n_),
    z_(ptf.z_),
    z0_(ptf.z0_),
    kappa_(ptf.kappa_)
{}


ABLInletVelocityFvPatchVectorField::ABLInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Ustar_(readScalar(dict.lookup("Ustar"))),
    n_(dict.lookup("n")),
    z_(dict.lookup("z")),
    z0_(readScalar(dict.lookup("z0"))),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41))
{
    if (mag(n_) < SMALL || mag(z_) < SMALL)
    {
        FatalErrorIn("ABLInletVelocityFvPatchVectorField(dict)")
            << "n or z given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    z_ /= mag(z_);

    evaluate();
}


ABLInletVelocityFvPatchVectorField::ABLInletVelocityFvPatchVectorField
(
    const ABLInletVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Ustar_(fcvpvf.Ustar_),
    n_(fcvpvf.n_),
    z_(fcvpvf.z_),
    z0_(fcvpvf.z0_),
    kappa_(fcvpvf.kappa_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ABLInletVelocityFvPatchVectorField::updateCoeffs()
{
    const vectorField& c = patch().Cf();
    scalarField coord = (c & z_);
    vectorField::operator=(n_*(Ustar_/kappa_)*log((coord + z0_)/z0_));
}


// Write
void ABLInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Ustar")
        << Ustar_ << token::END_STATEMENT << nl;
    os.writeKeyword("z0")
        << z0_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
     os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, ABLInletVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
