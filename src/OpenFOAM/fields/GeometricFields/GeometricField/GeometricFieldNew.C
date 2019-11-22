/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    const dimensionSet& ds,
    const word& patchFieldType
)
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        ds,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& iField,
    const word& patchFieldType
)
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        ds,
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
    const dimensionSet& ds,
    Field<Type>&& iField,
    const word& patchFieldType
)
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        ds,
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
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dt,
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
    return tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dt,
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& newName,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            newName,
            tgf().instance(),
            tgf().local(),
            tgf().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        tgf
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& newName,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            newName,
            tgf().instance(),
            tgf().local(),
            tgf().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        tgf,
        patchFieldTypes,
        actualPatchTypes
    );
}


// ************************************************************************* //
