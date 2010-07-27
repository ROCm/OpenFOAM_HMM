/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "multiphaseAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

multiphaseAlphaContactAngleFvPatchScalarField::interfaceThetaProps::
interfaceThetaProps
(
    Istream& is
)
:
    theta0_(readScalar(is)),
    uTheta_(readScalar(is)),
    thetaA_(readScalar(is)),
    thetaR_(readScalar(is))
{}


Istream& operator>>
(
    Istream& is,
    multiphaseAlphaContactAngleFvPatchScalarField::interfaceThetaProps& tp
)
{
    is  >> tp.theta0_ >> tp.uTheta_ >> tp.thetaA_ >> tp.thetaR_;
    return is;
}


Ostream& operator<<
(
    Ostream& os,
    const multiphaseAlphaContactAngleFvPatchScalarField::interfaceThetaProps& tp
)
{
    os  << tp.theta0_ << token::SPACE
        << tp.uTheta_ << token::SPACE
        << tp.thetaA_ << token::SPACE
        << tp.thetaR_;

    return os;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiphaseAlphaContactAngleFvPatchScalarField::
multiphaseAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF)
{}


multiphaseAlphaContactAngleFvPatchScalarField::
multiphaseAlphaContactAngleFvPatchScalarField
(
    const multiphaseAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(gcpsf, p, iF, mapper),
    thetaProps_(gcpsf.thetaProps_)
{}


multiphaseAlphaContactAngleFvPatchScalarField::
multiphaseAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF),
    thetaProps_(dict.lookup("thetaProperties"))
{
    evaluate();
}


multiphaseAlphaContactAngleFvPatchScalarField::
multiphaseAlphaContactAngleFvPatchScalarField
(
    const multiphaseAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(gcpsf, iF),
    thetaProps_(gcpsf.thetaProps_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void multiphaseAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("thetaProperties")
        << thetaProps_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    multiphaseAlphaContactAngleFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
