/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "uniformTemperatureJumpAMIFvPatchScalarField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformTemperatureJumpAMIFvPatchScalarField::
uniformTemperatureJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    uniformJumpAMIFvPatchField<scalar>(p, iF)
{}


Foam::uniformTemperatureJumpAMIFvPatchScalarField::
uniformTemperatureJumpAMIFvPatchScalarField
(
    const uniformTemperatureJumpAMIFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    uniformJumpAMIFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::uniformTemperatureJumpAMIFvPatchScalarField::
uniformTemperatureJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    uniformJumpAMIFvPatchField<scalar>(p, iF)
{}


Foam::uniformTemperatureJumpAMIFvPatchScalarField::
uniformTemperatureJumpAMIFvPatchScalarField
(
    const uniformTemperatureJumpAMIFvPatchScalarField& ptf
)
:
    uniformJumpAMIFvPatchField<scalar>(ptf)
{}


Foam::uniformTemperatureJumpAMIFvPatchScalarField::
uniformTemperatureJumpAMIFvPatchScalarField
(
    const uniformTemperatureJumpAMIFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    uniformJumpAMIFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       uniformTemperatureJumpAMIFvPatchScalarField
   );
}

// ************************************************************************* //
