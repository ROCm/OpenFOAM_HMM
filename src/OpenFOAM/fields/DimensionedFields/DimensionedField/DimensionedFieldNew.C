/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& iField
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>::NewImmovable
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
        iField
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& iField
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>::NewImmovable
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
        std::move(iField)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>::NewImmovable
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
        false  // checkIOFlags off
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims
)
{
    const bool caching = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>::NewImmovable
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
        false  // checkIOFlags off
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt
)
{
    return DimensionedField<Type, GeoMesh>::New
    (
        name,
        mesh,
        dt.value(),
        dt.dimensions()
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const tmp<DimensionedField<Type, GeoMesh>>& tfld
)
{
    const bool caching = tfld().db().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>::NewImmovable
    (
        caching,        // (true: immovable, false: movable)
        IOobject
        (
            name,
            tfld().instance(),
            tfld().local(),
            tfld().db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            caching     // (true: REGISTER, false: NO_REGISTER)
        ),
        tfld
    );
}


template<class Type, class GeoMesh>
template<class AnyType>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const DimensionedField<AnyType, GeoMesh>& fld,
    const word& name,
    const dimensionSet& dims
)
{
    const bool caching = fld.db().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>::NewImmovable
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
        dims
    );
}


template<class Type, class GeoMesh>
template<class AnyType>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const DimensionedField<AnyType, GeoMesh>& fld,
    const word& name,
    const dimensioned<Type>& dt
)
{
    const bool caching = fld.db().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>::NewImmovable
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
        dt.dimensions()
    );
}


// ************************************************************************* //
