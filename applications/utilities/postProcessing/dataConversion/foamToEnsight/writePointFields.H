/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Read point fields from disk
    and write as ensight data

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef ensight_writePointFields_H
#define ensight_writePointFields_H

#include "readFields.H"
#include "ensightMesh.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
bool writePointField
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const tmp<GeometricField<Type, pointPatchField, pointMesh>>& tfield
)
{
    if (!tfield.valid())
    {
        return false;
    }

    const auto& field = tfield();

    // PointData = true
    autoPtr<ensightFile> os = ensCase.newData<Type>(field.name(), true);

    bool wrote = ensightOutput::writePointField<Type>
    (
        os.ref(),
        field,
        ensMesh
    );

    tfield.clear();
    return wrote;
}


template<class Type>
label writePointFields
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const IOobjectList& objects
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> GeoField;

    const pointMesh& ptMesh = pointMesh::New(ensMesh.mesh());

    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        if
        (
            writePointField<Type>
            (
                ensCase,
                ensMesh,
                getField<GeoField>(ptMesh, objects, fieldName)
            )
        )
        {
            Info<< ' ' << fieldName;
            ++count;
        }
    }

    return count;
}


label writeAllPointFields
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const IOobjectList& objects
)
{
    #undef  ensight_WRITE_FIELD
    #define ensight_WRITE_FIELD(PrimitiveType)          \
        writePointFields<PrimitiveType>                 \
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
