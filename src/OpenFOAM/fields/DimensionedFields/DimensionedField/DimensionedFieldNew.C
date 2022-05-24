/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    const dimensionSet& ds,
    const Field<Type>& iField
)
{
    return tmp<DimensionedField<Type, GeoMesh>>::New
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
        iField
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& ds,
    Field<Type>&& iField
)
{
    return tmp<DimensionedField<Type, GeoMesh>>::New
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
        std::move(iField)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& ds
)
{
    return tmp<DimensionedField<Type, GeoMesh>>::New
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
        false  // checkIOFlags = true
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
    return tmp<DimensionedField<Type, GeoMesh>>::New
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
        false  // checkIOFlags = true
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& newName,
    const tmp<DimensionedField<Type, GeoMesh>>& tfld
)
{
    return tmp<DimensionedField<Type, GeoMesh>>::New
    (
        IOobject
        (
            newName,
            tfld().instance(),
            tfld().local(),
            tfld().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        tfld
    );
}


// ************************************************************************* //
