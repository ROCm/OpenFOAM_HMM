/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Read volume fields from disk and write with ensightParts

SourceFiles
    writeVolFields.H

\*---------------------------------------------------------------------------*/

#ifndef writeVolFields_H
#define writeVolFields_H

#include "readFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
bool writeVolField
(
    ensightCase& ensCase,
    const ensightParts& ensParts,
    const fvMesh& mesh,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tfield
)
{
    if (!tfield.valid())
    {
        return false;
    }

    const auto& field = tfield();

    autoPtr<ensightFile> os = ensCase.newData<Type>(field.name());

    // Currently serial only
    bool wrote = ensightOutput::Serial::writeVolField<Type>
    (
        field,
        ensParts,
        os.ref()
    );

    tfield.clear();
    return wrote;
}


template<class Type>
label writeVolFields
(
    ensightCase& ensCase,
    const ensightParts& ensParts,
    const fvMesh& mesh,
    const IOobjectList& objects
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        if
        (
            writeVolField<Type>
            (
                ensCase,
                ensParts,
                mesh,
                getField<GeoField>(objects.findObject(fieldName), mesh)
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
    const ensightParts& ensParts,
    const fvMesh& mesh,
    const IOobjectList& objects
)
{
    #undef  foamToEnsight_WRITE_FIELD
    #define foamToEnsight_WRITE_FIELD(PrimitiveType)    \
        writeVolFields<PrimitiveType>                   \
        (                                               \
            ensCase, ensParts,                          \
            mesh,                                       \
            objects                                     \
        )

    label count = 0;
    count += foamToEnsight_WRITE_FIELD(scalar);
    count += foamToEnsight_WRITE_FIELD(vector);
    count += foamToEnsight_WRITE_FIELD(sphericalTensor);
    count += foamToEnsight_WRITE_FIELD(symmTensor);
    count += foamToEnsight_WRITE_FIELD(tensor);

    #undef foamToEnsight_WRITE_FIELD
    return count;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
