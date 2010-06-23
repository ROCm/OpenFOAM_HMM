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

#include "atmBoundaryLayerInletEpsilonFvPatchScalarField.H"
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

atmBoundaryLayerInletEpsilonFvPatchScalarField::
atmBoundaryLayerInletEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Ustar_(0),
    z_(pTraits<vector>::zero),
    z0_(0),
    kappa_(0.41),
    zGround_(0)
{}


atmBoundaryLayerInletEpsilonFvPatchScalarField::
atmBoundaryLayerInletEpsilonFvPatchScalarField
(
    const atmBoundaryLayerInletEpsilonFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Ustar_(ptf.Ustar_),
    z_(ptf.z_),
    z0_(ptf.z0_),
    kappa_(ptf.kappa_),
    zGround_(ptf.zGround_)
{}


atmBoundaryLayerInletEpsilonFvPatchScalarField::
atmBoundaryLayerInletEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    Ustar_(readScalar(dict.lookup("Ustar"))),
    z_(dict.lookup("z")),
    z0_(readScalar(dict.lookup("z0"))),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    zGround_(readScalar(dict.lookup("zGround")))
{
    if (mag(z_) < SMALL)
    {
        FatalErrorIn("atmBoundaryLayerInletEpsilonFvPatchScalarField(dict)")
            << "z is not correct"
            << abort(FatalError);
    }

    z_ /= mag(z_);

    evaluate();
}


atmBoundaryLayerInletEpsilonFvPatchScalarField::
atmBoundaryLayerInletEpsilonFvPatchScalarField
(
    const atmBoundaryLayerInletEpsilonFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    Ustar_(fcvpvf.Ustar_),
    z_(fcvpvf.z_),
    z0_(fcvpvf.z0_),
    kappa_(fcvpvf.kappa_),
    zGround_(fcvpvf.zGround_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerInletEpsilonFvPatchScalarField::updateCoeffs()
{
    const vectorField& c = patch().Cf();
    scalarField coord = (c & z_);
    scalarField::operator=(pow(Ustar_, 3.0)/(kappa_*(coord - zGround_ + z0_)));
}


// Write
void atmBoundaryLayerInletEpsilonFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("Ustar")
        << Ustar_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    os.writeKeyword("z0")
        << z0_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("zGround")
        << zGround_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmBoundaryLayerInletEpsilonFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
