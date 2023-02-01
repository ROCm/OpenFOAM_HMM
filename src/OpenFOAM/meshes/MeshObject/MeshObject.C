/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "MeshObject.H"
#include "objectRegistry.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::MeshObject<Mesh, MeshObjectType, Type>::MeshObject(const Mesh& mesh)
:
    MeshObjectType<Mesh>(Type::typeName, mesh.thisDb()),
    mesh_(mesh)
{}


template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::MeshObject<Mesh, MeshObjectType, Type>::MeshObject
(
    const word& objName,
    const Mesh& mesh
)
:
    MeshObjectType<Mesh>(objName, mesh.thisDb()),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
template<class... Args>
const Type& Foam::MeshObject<Mesh, MeshObjectType, Type>::New
(
    const Mesh& mesh,
    Args&&... args
)
{
    Type* ptr =
        mesh.thisDb().objectRegistry::template
        getObjectPtr<Type>(Type::typeName);

    if (ptr)
    {
        return *ptr;
    }

    if (meshObject::debug)
    {
        Pout<< "MeshObject::New(const " << Mesh::typeName
            << "&, ...) : constructing " << Type::typeName
            << " for region " << mesh.name() << endl;
    }

    ptr = new Type(mesh, std::forward<Args>(args)...);

    regIOobject::store(static_cast<MeshObjectType<Mesh>*>(ptr));

    return *ptr;
}


template<class Mesh, template<class> class MeshObjectType, class Type>
template<class... Args>
const Type& Foam::MeshObject<Mesh, MeshObjectType, Type>::New
(
    const word& objName,
    const Mesh& mesh,
    Args&&... args
)
{
    Type* ptr =
        mesh.thisDb().objectRegistry::template
        getObjectPtr<Type>(objName);

    if (ptr)
    {
        return *ptr;
    }

    if (meshObject::debug)
    {
        Pout<< "MeshObject::New('" << objName
            << "', const " << Mesh::typeName
            << "&, ...) : constructing " << objName
            << " of type " << Type::typeName
            << " for region " << mesh.name() << endl;
    }

    ptr = new Type(objName, mesh, std::forward<Args>(args)...);

    regIOobject::store(static_cast<MeshObjectType<Mesh>*>(ptr));

    return *ptr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
bool Foam::MeshObject<Mesh, MeshObjectType, Type>::Delete
(
    const word& objName,
    const Mesh& mesh
)
{
    Type* ptr =
        mesh.thisDb().objectRegistry::template
        getObjectPtr<Type>(objName);

    if (ptr)
    {
        if (meshObject::debug)
        {
            Pout<< "MeshObject::Delete(const Mesh&) : deleting "
                << objName << endl;
        }

        return mesh.thisDb().checkOut(static_cast<MeshObjectType<Mesh>*>(ptr));
    }

    return false;
}


template<class Mesh>
void Foam::meshObject::movePoints(objectRegistry& obr)
{
    HashTable<GeometricMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<GeometricMeshObject<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::movePoints(objectRegistry&) :"
            << " moving " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIters(meshObjects, iter)
    {
        // isA<MoveableMeshObject<Mesh>>
        auto* objectPtr = dynamic_cast<MoveableMeshObject<Mesh>*>(*iter);

        if (objectPtr)
        {
            if (meshObject::debug)
            {
                Pout<< "    Moving " << (*iter)->name() << endl;
            }
            objectPtr->movePoints();
        }
        else
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << (*iter)->name() << endl;
            }
            obr.checkOut(*iter);
        }
    }
}


template<class Mesh>
void Foam::meshObject::updateMesh(objectRegistry& obr, const mapPolyMesh& mpm)
{
    HashTable<GeometricMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<GeometricMeshObject<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::updateMesh(objectRegistry&, "
               "const mapPolyMesh& mpm) : updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIters(meshObjects, iter)
    {
        // isA<MoveableMeshObject<Mesh>>
        auto* objectPtr = dynamic_cast<UpdateableMeshObject<Mesh>*>(*iter);

        if (objectPtr)
        {
            if (meshObject::debug)
            {
                Pout<< "    Updating " << (*iter)->name() << endl;
            }
            objectPtr->updateMesh(mpm);
        }
        else
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << (*iter)->name() << endl;
            }
            obr.checkOut(*iter);
        }
    }
}


template<class Mesh, template<class> class MeshObjectType>
void Foam::meshObject::clear(objectRegistry& obr)
{
    HashTable<MeshObjectType<Mesh>*> meshObjects
    (
        obr.lookupClass<MeshObjectType<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::clear(objectRegistry&) :"
            << " clearing " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIters(meshObjects, iter)
    {
        if (meshObject::debug)
        {
            Pout<< "    Destroying " << (*iter)->name() << endl;
        }
        obr.checkOut(*iter);
    }
}


template
<
    class Mesh,
    template<class> class FromType,
    template<class> class ToType
>
void Foam::meshObject::clearUpto(objectRegistry& obr)
{
    HashTable<FromType<Mesh>*> meshObjects
    (
        obr.lookupClass<FromType<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::clearUpto(objectRegistry&) :"
            << " clearing " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIters(meshObjects, iter)
    {
        // isA<ToType<Mesh>
        auto* objectPtr = dynamic_cast<ToType<Mesh>*>(*iter);

        if (!objectPtr)
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << (*iter)->name() << endl;
            }
            obr.checkOut(*iter);
        }
    }
}


// ************************************************************************* //
