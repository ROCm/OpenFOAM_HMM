/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "MeshObject.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mesh, class Type>
Foam::MeshObject<Mesh, Type>::MeshObject(const Mesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            Type::typeName,
            mesh.db().instance(),
            mesh.db()
        )
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Mesh, class Type>
const Type& Foam::MeshObject<Mesh, Type>::New
(
    const Mesh& mesh
)
{
    if (!mesh.db().objectRegistry::foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh));
    }
    else
    {
        return mesh.db().objectRegistry::lookupObject<Type>(Type::typeName);
    }
}


template<class Mesh, class Type>
template<class Data1>
const Type& Foam::MeshObject<Mesh, Type>::New
(
    const Mesh& mesh,
    const Data1& d
)
{
    if (!mesh.db().objectRegistry::foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh, d));
    }
    else
    {
        return mesh.db().objectRegistry::lookupObject<Type>(Type::typeName);
    }
}


template<class Mesh, class Type>
template<class Data1, class Data2>
const Type& Foam::MeshObject<Mesh, Type>::New
(
    const Mesh& mesh,
    const Data1& d1,
    const Data2& d2
)
{
    if (!mesh.db().objectRegistry::foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh, d1, d2));
    }
    else
    {
        return mesh.db().objectRegistry::lookupObject<Type>(Type::typeName);
    }
}


template<class Mesh, class Type>
template<class Data1, class Data2, class Data3>
const Type& Foam::MeshObject<Mesh, Type>::New
(
    const Mesh& mesh,
    const Data1& d1,
    const Data2& d2,
    const Data3& d3
)
{
    if (!mesh.db().objectRegistry::foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh, d1, d2, d3));
    }
    else
    {
        return mesh.db().objectRegistry::lookupObject<Type>(Type::typeName);
    }
}


template<class Mesh, class Type>
template<class Data1, class Data2, class Data3, class Data4>
const Type& Foam::MeshObject<Mesh, Type>::New
(
    const Mesh& mesh,
    const Data1& d1,
    const Data2& d2,
    const Data3& d3,
    const Data4& d4
)
{
    if (!mesh.db().objectRegistry::foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh, d3, d4));
    }
    else
    {
        return mesh.db().objectRegistry::lookupObject<Type>(Type::typeName);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Mesh, class Type>
bool Foam::MeshObject<Mesh, Type>::Delete(const Mesh& mesh)
{
    if (mesh.db().objectRegistry::foundObject<Type>(Type::typeName))
    {
        return mesh.db().objectRegistry::checkOut
        (
            const_cast<Type&>
            (
                mesh.db().objectRegistry::lookupObject<Type>
                (
                    Type::typeName
                )
            )
        );
    }
    else
    {
        return false;
    }
}


template<class Mesh, class Type>
Foam::MeshObject<Mesh, Type>::~MeshObject()
{
    release();
}


// ************************************************************************* //
