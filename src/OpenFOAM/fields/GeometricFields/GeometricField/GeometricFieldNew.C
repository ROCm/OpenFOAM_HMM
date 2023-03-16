/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    const word& patchFieldType
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        mesh,
        dims,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& iField,
    const word& patchFieldType
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        mesh,
        dims,
        iField,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& iField,
    const word& patchFieldType
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        mesh,
        dims,
        std::move(iField),
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims,
    const word& patchFieldType
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        mesh,
        value,
        dims,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        mesh,
        value,
        dims,
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        mesh,
        dt.value(),
        dt.dimensions(),
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        mesh,
        dt.value(),
        dt.dimensions(),
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
    const word& patchFieldType
)
{
    const bool caching = tgf().db().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            tgf().instance(),
            tgf().local(),
            tgf().db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        tgf,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    const bool caching = tgf().db().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            tgf().instance(),
            tgf().local(),
            tgf().db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        tgf
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    const bool caching = tgf().db().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            tgf().instance(),
            tgf().local(),
            tgf().db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        tgf,
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class AnyType>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const GeometricField<AnyType, PatchField, GeoMesh>& fld,
    const word& name,
    const dimensionSet& dims,
    const word& patchFieldType
)
{
    const bool caching = fld.db().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            fld.instance(),
            fld.db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        fld.mesh(),
        dims,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class AnyType>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const GeometricField<AnyType, PatchField, GeoMesh>& fld,
    const word& name,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
{
    const bool caching = fld.db().cacheTemporaryObject(name);

    return tmp<GeometricField<Type, PatchField, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            fld.instance(),
            fld.db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        fld.mesh(),
        dt.value(),
        dt.dimensions(),
        patchFieldType
    );
}


// ************************************************************************* //
