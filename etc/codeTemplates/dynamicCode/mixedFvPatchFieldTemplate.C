/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) YEAR AUTHOR,AFFILIATION
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

#include "mixedFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = ${SHA1sum}
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void ${typeName}_${SHA1sum}(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatch${FieldType},
    ${typeName}MixedValueFvPatch${FieldType}
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF
)
:
    mixedFvPatchField<${TemplateType}>(p, iF)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} : patch/DimensionedField");
    }
}


${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const ${typeName}MixedValueFvPatch${FieldType}& ptf,
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<${TemplateType}>(ptf, p, iF, mapper)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} : patch/DimensionedField/mapper");
    }
}


${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<${TemplateType}>(p, iF, dict)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} : patch/dictionary");
    }
}


${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const ${typeName}MixedValueFvPatch${FieldType}& ptf
)
:
    mixedFvPatchField<${TemplateType}>(ptf)
{
    if (${verbose:-false})
    {
        printMessage("Copy construct ${typeName}");
    }
}


${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const ${typeName}MixedValueFvPatch${FieldType}& ptf,
    const DimensionedField<${TemplateType}, volMesh>& iF
)
:
    mixedFvPatchField<${TemplateType}>(ptf, iF)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}MixedValueFvPatch${FieldType}::
~${typeName}MixedValueFvPatch${FieldType}()
{
    if (${verbose:-false})
    {
        printMessage("Destroy ${typeName}");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}MixedValueFvPatch${FieldType}::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (${verbose:-false})
    {
        printMessage("updateCoeffs ${typeName}");
    }

//{{{ begin code
    ${code}
//}}} end code

    this->mixedFvPatchField<${TemplateType}>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
