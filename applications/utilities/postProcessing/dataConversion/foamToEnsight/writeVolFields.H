/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Read volume fields from disk and write with ensightMesh

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef ensight_writeVolFields_H
#define ensight_writeVolFields_H

#include "readFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
bool writeVolField
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tfield,
    const bool nearCellValue = false
)
{
    if (!tfield.valid())
    {
        return false;
    }
    else if (nearCellValue)
    {
        auto tzgrad = makeZeroGradientField<Type>(tfield);

        // Recursive call
        return writeVolField
        (
            ensCase,
            ensMesh,
            tzgrad,
            false  // No nearCellValue, we already have zero-gradient
        );
    }

    const auto& field = tfield();

    // Forced use of node values?
    const bool nodeValues = ensCase.nodeValues();

    autoPtr<ensightFile> os =
        ensCase.newData<Type>(field.name(), nodeValues);

    bool wrote = ensightOutput::writeVolField<Type>
    (
        os.ref(),
        field,
        ensMesh,
        nodeValues
    );

    tfield.clear();
    return wrote;
}


template<class Type>
label writeVolFields
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const IOobjectList& objects,
    const bool nearCellValue = false
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    const fvMesh& mesh = dynamicCast<const fvMesh>(ensMesh.mesh());

    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        if
        (
            writeVolField<Type>
            (
                ensCase,
                ensMesh,
                getField<GeoField>(objects.findObject(fieldName), mesh),
                nearCellValue
            )
        )
        {
            Info<< ' ' << fieldName;
            ++count;
        }
    }

    return count;
}


label writeAllVolFields
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const IOobjectList& objects,
    const bool nearCellValue = false
)
{
    #undef  ensight_WRITE_FIELD
    #define ensight_WRITE_FIELD(PrimitiveType)          \
        writeVolFields<PrimitiveType>                   \
        (                                               \
            ensCase,                                    \
            ensMesh,                                    \
            objects,                                    \
            nearCellValue                               \
        )

    label count = 0;
    count += ensight_WRITE_FIELD(scalar);
    count += ensight_WRITE_FIELD(vector);
    count += ensight_WRITE_FIELD(sphericalTensor);
    count += ensight_WRITE_FIELD(symmTensor);
    count += ensight_WRITE_FIELD(tensor);

    #undef ensight_WRITE_FIELD
    return count;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
