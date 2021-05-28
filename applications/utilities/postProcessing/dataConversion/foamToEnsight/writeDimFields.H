/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

InNamespace
    Foam

Description
    Read dimensioned fields from disk and write with ensightMesh

\*---------------------------------------------------------------------------*/

#ifndef ensight_writeDimFields_H
#define ensight_writeDimFields_H

#include "writeVolFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
bool writeDimField
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const tmp<DimensionedField<Type, volMesh>>& tdf
)
{
    if (!tdf.valid())
    {
        return false;
    }

    auto tfield = makeZeroGradientField<Type>(tdf);

    // Now a volField with zero-gradient boundaries

    return writeVolField<Type>
    (
        ensCase,
        ensMesh,
        tfield,
        false  // No nearCellValue, we already have zero-gradient
    );
}


template<class Type>
label writeDimFields
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const IOobjectList& objects
)
{
    typedef typename
        GeometricField
        <
            Type, fvPatchField, volMesh
        >::Internal DimField;

    const fvMesh& mesh = dynamicCast<const fvMesh>(ensMesh.mesh());

    label count = 0;

    for (const word& fieldName : objects.sortedNames<DimField>())
    {
        if
        (
            writeDimField<Type>
            (
                ensCase,
                ensMesh,
                getField<DimField>(objects.findObject(fieldName), mesh)
            )
        )
        {
            Info<< ' ' << fieldName;
            ++count;
        }
    }

    return count;
}


label writeAllDimFields
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const IOobjectList& objects
)
{
    #undef  ensight_WRITE_FIELD
    #define ensight_WRITE_FIELD(PrimitiveType)          \
        writeDimFields<PrimitiveType>                   \
        (                                               \
            ensCase,                                    \
            ensMesh,                                    \
            objects                                     \
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
