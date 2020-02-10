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
    Read dimensioned fields from disk and write with ensightMesh

SourceFiles

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
