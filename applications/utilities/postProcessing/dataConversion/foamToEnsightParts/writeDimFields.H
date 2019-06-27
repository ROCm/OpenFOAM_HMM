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
    Read dimensioned fields from disk and write with ensightParts

SourceFiles
    writeDimFields.H

\*---------------------------------------------------------------------------*/

#ifndef ensightParts_writeDimFields_H
#define ensightParts_writeDimFields_H

#include "writeVolFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
label writeDimFields
(
    ensightCase& ensCase,
    const ensightParts& ensParts,
    const fvMesh& mesh,
    const IOobjectList& objects
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    typedef typename
        GeometricField
        <
            Type, fvPatchField, volMesh
        >::Internal DimField;


    label count = 0;

    for (const word& fieldName : objects.sortedNames<DimField>())
    {
        const IOobject* io = objects.findObject(fieldName);

        if
        (
            writeVolField<Type>
            (
                ensCase,
                ensParts,
                mesh,
                getZeroGradField<GeoField>(io, mesh)
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
    const ensightParts& ensParts,
    const fvMesh& mesh,
    const IOobjectList& objects
)
{
    #undef  foamToEnsight_WRITE_FIELD
    #define foamToEnsight_WRITE_FIELD(PrimitiveType)    \
        writeDimFields<PrimitiveType>                   \
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
