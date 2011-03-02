/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "fixedValueFvPatchScalarFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    bool ${typeName}_${SHA1sum}()
    {
        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}FixedValueFvPatchScalarField::
${typeName}FixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/DimensionedField\n";
    }
}


${typeName}FixedValueFvPatchScalarField::
${typeName}FixedValueFvPatchScalarField
(
    const ${typeName}FixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/DimensionedField/mapper\n";
    }
}


${typeName}FixedValueFvPatchScalarField::
${typeName}FixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/dictionary\n";
    }
}


${typeName}FixedValueFvPatchScalarField::
${typeName}FixedValueFvPatchScalarField
(
    const ${typeName}FixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " as copy\n";
    }
}


${typeName}FixedValueFvPatchScalarField::
${typeName}FixedValueFvPatchScalarField
(
    const ${typeName}FixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum} "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}FixedValueFvPatchScalarField::
~${typeName}FixedValueFvPatchScalarField()
{
    if (${verbose:-false})
    {
        Info<<"destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}FixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (${verbose:-false})
    {
        Info<<"updateCoeffs ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${code};
//}}} end code

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    ${typeName}FixedValueFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
